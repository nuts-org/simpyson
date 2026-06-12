from __future__ import annotations

import contextlib
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from simpyson.converter import ppm2hz
from simpyson.io import read_simp
from simpyson.templates import (
    CPMAS,
    CustomPulseSequence,
    PulseSequenceTemplate,
    get_template,
    pulseq_templates,
)
from simpyson.utils import get_larmor_freq, get_spin

logger = logging.getLogger("simpyson")


def _proton_freq_to_b0(proton_freq):
    """
    Convert a proton_frequency value to a b0 string like '800.0MHz'.

    Parameters
    ----------
    proton_freq : int, float, or str
        Proton frequency in Hz (numeric) or as a string with units.

    Returns
    -------
    str or None
        B0 string suitable for ``hz2ppm`` / ``ppm2hz``, or None if conversion
        is not possible.
    """
    unit_to_hz = {'hz': 1.0, 'khz': 1e3, 'mhz': 1e6, 'ghz': 1e9, 'thz': 1e12}

    if isinstance(proton_freq, (int, float)):
        if proton_freq > 1e6:
            return f"{proton_freq / 1e6:.1f}MHz"
        # Small numeric values are assumed to already be in MHz
        return f"{float(proton_freq)}MHz"
    if isinstance(proton_freq, str):
        match = re.match(
            r'(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([kMGT]?Hz)?\s*$',
            proton_freq.strip(),
            re.IGNORECASE,
        )
        if match:
            value, unit = match.groups()
            if not unit:
                # Bare numeric strings follow the numeric rule:
                # > 1e6 means Hz (e.g. '8e8' -> 800 MHz), otherwise MHz.
                return _proton_freq_to_b0(float(value))
            # Normalise any Hz-based unit to MHz, since get_larmor_freq()
            # only accepts 'T' or 'MHz'.
            factor = unit_to_hz.get(unit.lower())
            if factor is not None:
                return f"{float(value) * factor / 1e6}MHz"
    return None


def _extract_nucleus(spinsys_str):
    """
    Extract the observed nucleus from a SIMPSON spinsys string.

    Looks for the first nucleus on the ``channels`` line, falling back to the
    first nucleus on the ``nuclei`` line.

    Parameters
    ----------
    spinsys_str : str
        SIMPSON spinsys block as a string.

    Returns
    -------
    str or None
        Nucleus string (e.g. ``'1H'``, ``'13C'``), or None if not found.
    """
    # [^\n]+ instead of [\w\s]+: \s matches newlines, which made the old
    # pattern swallow the following spinsys lines as well.
    channels_match = re.search(r'channels[ \t]+([^\n]+)', spinsys_str)
    if channels_match:
        nuclei_list = channels_match.group(1).split()
        if nuclei_list:
            return nuclei_list[0]
    nuclei_match = re.search(r'nuclei\s+(\w+)', spinsys_str)
    if nuclei_match:
        return nuclei_match.group(1)
    return None


def _extract_turnoff_interactions(spinsys_str: str) -> list[str]:
    """
    Return the interaction names present in a spinsys block.

    Scans for ``dipole`` and ``jcoupling`` lines and returns the names in the
    form expected by SIMPSON's ``turnoff`` command (e.g. ``dipole_1_2``).
    """
    interactions = []
    pattern = re.compile(r'(dipole|jcoupling)\s+(\d+)\s+(\d+)', re.IGNORECASE)
    for line in spinsys_str.splitlines():
        m = pattern.match(line.strip())
        if m:
            interactions.append(f"{m.group(1).lower()}_{m.group(2)}_{m.group(3)}")
    return interactions


def _is_ct_operator(detect_op: str) -> bool:
    """Return True if detect_op selects only the central transition."""
    return bool(re.search(r'I(?:n|\d+)c', detect_op))


class SimpCalc:
    """
    Generator for SIMPSON simulation input files (``.tcl``).

    Assembles the four main sections of a SIMPSON input file (spinsys, par,
    pulseq, main) from a spin system definition, pulse sequence template,
    and simulation parameters.

    Parameters
    ----------
    spinsys : str or object
        Spin system definition. Can be a SIMPSON spinsys block string, the
        body of a spinsys block (starting with ``channels`` or ``nuclei``),
        or a Soprano SpinSystem object with a ``to_simpson()`` method.
    pulse_sequence : str, PulseSequenceTemplate, or None
        Pulse sequence to use. Can be a template name (``'no_pulse'``,
        ``'pulse_90'``, ``'cp_mas'``), a custom Tcl code string, a
        ``PulseSequenceTemplate`` instance, or None for no pulse sequence.
    **kwargs
        Simulation parameters. Required: ``proton_frequency``, ``spin_rate``,
        ``start_operator``, ``detect_operator``, ``np``, ``sw``, ``method``,
        ``crystal_file``, ``gamma_angles``, ``verbose``.

        Output options: ``out_name``, ``out_format``, ``lb``, ``gauss_lb``,
        ``zerofill``. Note that ``zerofill`` defaults to ``np`` (i.e. no
        zero-filling); pass e.g. ``zerofill=2*np`` to interpolate the
        spectrum.

    Raises
    ------
    ValueError
        If ``spinsys`` is None or ``pulse_sequence`` has an invalid type.

    Examples
    --------
    >>> calc = SimpCalc(
    ...     spinsys="channels 1H\\nnuclei 1H\\nshift 1 5p 0 0 0 0 0",
    ...     pulse_sequence='no_pulse',
    ...     proton_frequency=400e6, spin_rate=10000, sw=20000, np=1024,
    ...     start_operator='Inx', detect_operator='Inp', method='direct',
    ...     crystal_file='rep100', gamma_angles=10, verbose=0,
    ... )
    >>> calc.save("simulation.in")
    """

    @staticmethod
    def _normalize_spinsys(spinsys: str | object) -> str:
        """Return a canonical ``spinsys { … }`` block string."""
        if spinsys is None:
            raise ValueError("spinsys cannot be None")
        if hasattr(spinsys, 'to_simpson'):
            spinsys = spinsys.to_simpson()
        if not isinstance(spinsys, str):
            raise ValueError(
                f"spinsys must be a string or a Soprano SpinSystem object, got {type(spinsys)}"
            )
        stripped = spinsys.strip()
        if stripped.startswith("spinsys"):
            return stripped + "\n"
        return f"spinsys {{\n{stripped}\n}}\n"

    @property
    def spinsys(self) -> str:
        """The normalized ``spinsys { … }`` block string."""
        return self._spinsys

    @spinsys.setter
    def spinsys(self, value: str | object) -> None:
        self._spinsys = self._normalize_spinsys(value)
        if hasattr(self, 'pulse_sequence') and isinstance(self.pulse_sequence, CPMAS):
            self.pulse_sequence.turnoff_interactions = _extract_turnoff_interactions(self._spinsys)

    def __init__(self, spinsys: str | object, pulse_sequence: str | PulseSequenceTemplate | None = None, **kwargs) -> None:
        self.spinsys = spinsys  # calls setter, raises on bad input
        self.parameters = kwargs
        self.output_config = {}

        output_keys = ['out_name', 'out_format', 'lb', 'zerofill', 'gauss_lb']
        for key in output_keys:
            if key in self.parameters:
                self.output_config[key.replace('out_', '')] = self.parameters.pop(key)

        #   hz = f_SPE - ref = (d_Hz - OFFSET) - (-OFFSET) = d_Hz
        if 'variable_offset' in self.parameters and 'variable_ref' not in self.parameters:
            self.parameters['variable_ref'] = -self.parameters['variable_offset']

        self.pulse_sequence = self._setup_pulse_sequence(pulse_sequence)

    def __str__(self) -> str:
        """Generate the complete SIMPSON input file as a string."""
        sections = []
        sections.append(self.generate_spinsys())
        sections.append(self.generate_par())
        sections.append(self.generate_pulseq())
        sections.append(self.generate_main())
        return "\n".join(sections)


    def _setup_pulse_sequence(self, pulse_sequence: str | PulseSequenceTemplate | None) -> PulseSequenceTemplate | None:
        """Set up the pulse sequence based on user input."""
        if pulse_sequence is None:
             return None

        if isinstance(pulse_sequence, str):
            if pulse_sequence in pulseq_templates:
                template_class = pulseq_templates[pulse_sequence]
                template_instance = template_class()

                # Get variable parameters
                required_clean = {param.replace('variable_', '')
                                 for param in template_instance.get_required_parameters()}

                # Extract matching parameters from user input
                pulseq_params = {}
                for param in required_clean:
                    if param in self.parameters:
                        pulseq_params[param] = self.parameters[param]

                # Pass offset if it exists
                if 'variable_offset' in self.parameters:
                    pulseq_params['offset'] = self.parameters['variable_offset']

                # Count channels so the template emits the right number of offset args
                channels_match = re.search(r'channels\s+([^\n]+)', self.spinsys)
                if channels_match:
                    pulseq_params['num_channels'] = len(channels_match.group(1).split())

                # Create template with extracted parameters
                return get_template(pulse_sequence, **pulseq_params)
            else:
                # Custom string sequence
                return CustomPulseSequence(pulse_sequence, **self.parameters)

        elif isinstance(pulse_sequence, PulseSequenceTemplate):
            return pulse_sequence

        else:
            raise ValueError(
                f"Invalid pulse_sequence type: {type(pulse_sequence)}. "
                "Must be a template name (str), a custom string, or "
                "a PulseSequenceTemplate object"
            )

    def generate_spinsys(self) -> str:
        """
        Return the spinsys section of the SIMPSON input file.

        The spinsys is normalised once at construction (or when the
        ``spinsys`` property is set) and returned verbatim here.

        Returns
        -------
        str
            The ``spinsys { … }`` block as a string.
        """
        return self.spinsys

    def generate_par(self) -> str:
        """
        Generate the par section of the SIMPSON input file.

        Returns
        -------
        str
            The par section as a string.

        Raises
        ------
        ValueError
            If required simulation parameters are missing.
        """

        # Parameters required for every simulation
        required_params = {
            "proton_frequency", "spin_rate", "start_operator", "detect_operator",
            "np", "sw", "method", "crystal_file", "gamma_angles", "verbose"
        }


        missing_params = required_params - set(self.parameters.keys())
        if missing_params:
            raise ValueError(f"Missing required parameters: {', '.join(missing_params)}")

        par_block = "par {\n"


        for param in sorted(required_params):
            if param in self.parameters:
                value = self.parameters[param]
                par_block += f"   {param:<20} {value}\n"

        # Add pulse sequence parameters as variables
        pulse_sequence_vars = set()
        if self.pulse_sequence:
            for param, value in sorted(self.pulse_sequence.parameters.items()):
                if param.startswith('variable_'):
                    var_name = param.replace('variable_', '')
                    pulse_sequence_vars.add(var_name)
                    # Only add as variable if it's not already a standard parameter
                    if var_name not in required_params:
                        par_block += f"   variable {var_name:<15} {value}\n"

        # Add other variables from parameters
        for param, value in sorted(self.parameters.items()):
            if param.startswith('variable_'):
                var_name = param.replace('variable_', '')
                if var_name not in pulse_sequence_vars:
                    par_block += f"   variable {var_name:<15} {value}\n"

        # Add any remaining parameters
        processed_params = required_params | {k for k in self.parameters if k.startswith('variable_')} | pulse_sequence_vars
        remaining_params = {k: v for k, v in self.parameters.items()
                          if k not in processed_params}

        for param, value in sorted(remaining_params.items()):
            if not param.startswith('out_'):
                par_block += f"   {param:<20} {value}\n"

        par_block += "}\n"
        return par_block

    def generate_pulseq(self) -> str:
        """
        Generate the pulseq section of the SIMPSON input file.

        Returns
        -------
        str
            The pulseq section as a string.
        """
        if not self.pulse_sequence:
            # Return empty pulseq block if no sequence provided
            return "proc pulseq {} {}\n"

        # Lazily populate CPMAS turnoff list from the spinsys so this works
        # even when a CPMAS instance is passed directly to SimpCalc.
        if isinstance(self.pulse_sequence, CPMAS) and not self.pulse_sequence.turnoff_interactions:
            self.pulse_sequence.turnoff_interactions = _extract_turnoff_interactions(self.spinsys)

        return self.pulse_sequence.generate_code()

    def generate_main(self) -> str:
        """
        Generate the main section of the SIMPSON input file.

        Returns
        -------
        str
            The main section as a string.

        Raises
        ------
        ValueError
            If ``out_format`` is not one of ``'fid'``, ``'spe'``, ``'xreim'``.
        """

        # Output settings are popped from kwargs into output_config in __init__
        out_format = self.output_config.get('format', 'spe')
        out_name = self.output_config.get('name', '$par(name)')
        lb = self.output_config.get('lb', 20)
        gauss_lb = self.output_config.get('gauss_lb', 0)
        zerofill = self.output_config.get('zerofill', self.parameters.get('np', 0))


        indent = "    "

        gauss_line = f"{indent}faddlb $f {gauss_lb} 1\n" if gauss_lb else ""

        if out_format == "fid":
            return f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}faddlb $f {lb} 0
{gauss_line}{indent}fzerofill $f {zerofill}
{indent}fsave $f {out_name}.fid
}}
"""
        elif out_format == "spe":
            main_block = f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}faddlb $f {lb} 0
{gauss_line}{indent}fzerofill $f {zerofill}
{indent}fft $f
"""
            if 'variable_ref' in self.parameters or 'ref' in self.parameters:
                main_block += f"{indent}fset $f -ref $par(ref)\n"

            main_block += f"{indent}fsave $f {out_name}.spe\n}}\n"
            return main_block
        elif out_format == "xreim":
            return f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}fsave $f {out_name}.xreim -xreim
}}
"""
        else:
            raise ValueError(f"Unknown out_format '{out_format}'. Supported formats: 'fid', 'spe', 'xreim'")

    def save(self, filepath: str) -> None:
        """
        Save the SIMPSON input file to disk.

        Parameters
        ----------
        filepath : str
            Path to write the ``.in`` / ``.tcl`` file.
        """
        with Path(filepath).open('w') as file:
            file.write(str(self))

    def print(self) -> None:
        """Print the SIMPSON input file to the console (stdout)."""
        print(str(self))  # noqa: T201

    def run(
        self,
        filepath: str | None = None,
        timeout: int | None = None,
        read_output: bool = True,
        delete_files: bool = True,
        b0: str | None = None,
        nucleus: str | None = None,
        simpson_path: str | None = None,
        dry_run: bool = False,
    ) -> str | object:
        """
        Run the SIMPSON simulation and optionally read the results.

        Parameters
        ----------
        filepath : str or None
            Path to save the input file. If None, a temporary file is created.
        timeout : int or None
            Timeout in seconds for the SIMPSON process.
        read_output : bool
            If True (default), read the output file and return a Simpy object.
            If False, return SIMPSON's stdout as a string.
        delete_files : bool
            If True (default), delete input and output files after reading.
        b0 : str or None
            Magnetic field strength (e.g., ``'400MHz'``). Auto-derived from
            ``proton_frequency`` if not provided. Only used with
            ``read_output=True``.
        nucleus : str or None
            Nucleus type (e.g., ``'1H'``). Auto-extracted from the spinsys
            if not provided. Only used with ``read_output=True``.
        simpson_path : str or None
            Custom path to the SIMPSON executable.
        dry_run : bool
            If True, generate the input file and return the command string
            without running SIMPSON.

        Returns
        -------
        str or Simpy
            The SIMPSON command (if ``dry_run``), stdout output (if not
            ``read_output``), or a ``Simpy`` object with simulation results.

        Raises
        ------
        FileNotFoundError
            If SIMPSON is not found in PATH or at the specified path.
        subprocess.TimeoutExpired
            If the simulation exceeds the timeout.
        subprocess.CalledProcessError
            If SIMPSON returns a non-zero exit code.
        """
        # Find the SIMPSON executable
        simpson_executable = None

        if simpson_path:
            if Path(simpson_path).exists():
                simpson_executable = simpson_path
            else:
                raise FileNotFoundError(f"SIMPSON executable not found at specified path: {simpson_path}")

        if simpson_executable is None:
            simpson_executable = shutil.which("simpson")

        if simpson_executable is None and not dry_run:
            raise FileNotFoundError(
                "SIMPSON executable not found in PATH or at the specified path.\n"
                "Please ensure SIMPSON is installed and available in your PATH, or provide the correct path."
            )

        temp_file = None
        if filepath is None:
            temp_fd, temp_file = tempfile.mkstemp(suffix='.in')
            os.close(temp_fd)
            filepath = temp_file

        base_filepath = str(Path(filepath).with_suffix(''))

        self.save(filepath)

        if dry_run:
            cmd = [simpson_executable if simpson_executable else "simpson", filepath]
            logger.info("Dry run: Generated input file at %s", filepath)
            logger.info("Command: %s", ' '.join(cmd))
            return ' '.join(cmd)

        # Determine expected output filename/locations
        # (output settings are popped from kwargs into output_config in __init__)
        out_format = self.output_config.get('format', 'spe')

        out_name = self.output_config.get('name', base_filepath)
        if out_name == '$par(name)':
            out_name = base_filepath

        output_filename = f"{Path(out_name).name}.{out_format}"

        possible_locations = [
            str(Path(filepath).parent / output_filename),
            str(Path.cwd() / output_filename),
            f"{out_name}.{out_format}"
        ]

        # Remember which candidate output files already exist, so cleanup
        # never deletes a pre-existing user file that happens to share a name.
        preexisting = {loc for loc in possible_locations if Path(loc).exists()}

        try:
            cmd = [simpson_executable, filepath]

            try:
                result = subprocess.run(cmd,
                                     check=True,
                                     capture_output=True,
                                     text=True,
                                     timeout=timeout)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"SIMPSON failed with exit code {e.returncode}.\n"
                    f"SIMPSON stderr:\n{e.stderr}\n"
                    f"SIMPSON stdout:\n{e.stdout}"
                ) from e

            # Read output from run
            if read_output:
                output_file = None
                for location in possible_locations:
                    if Path(location).exists():
                        output_file = location
                        break

                if output_file is None:
                    raise FileNotFoundError(
                        f"SIMPSON output file not found. Looked in: {possible_locations}\n"
                        f"SIMPSON stdout: {result.stdout}\n"
                        f"SIMPSON stderr: {result.stderr}"
                    )

                # Auto-extract magnetic field if not provided
                if b0 is None and 'proton_frequency' in self.parameters:
                    b0 = _proton_freq_to_b0(self.parameters['proton_frequency'])

                # Auto-extract nucleus information if not provided.
                # Use detect_operator (e.g. 'I2p') to identify, nucleus compare with spinsys nuclei list.
                # Falls back to the first channel for global operators ('Inp', 'Inc').
                if nucleus is None:
                    spinsys_str = self.generate_spinsys()
                    detect_op = self.parameters.get('detect_operator', '')
                    indices = re.findall(r'I(\d+)', detect_op)
                    if indices:
                        nuclei_match = re.search(r'nuclei[ \t]+([^\n]+)', spinsys_str)
                        if nuclei_match:
                            nuclei_list = nuclei_match.group(1).split()
                            idx = int(indices[0]) - 1  # SIMPSON is 1-indexed
                            if 0 <= idx < len(nuclei_list):
                                nucleus = nuclei_list[idx]
                    if nucleus is None:
                        nucleus = _extract_nucleus(spinsys_str)

                # Read the output file
                return read_simp(output_file, format=out_format, b0=b0, nucleus=nucleus)
            else:
                return result.stdout

        finally:
            if delete_files:
                if temp_file and Path(temp_file).exists():
                    Path(temp_file).unlink()
                elif filepath and Path(filepath).exists():
                    Path(filepath).unlink()

                for location in possible_locations:
                    if location in preexisting:
                        continue
                    with contextlib.suppress(OSError):
                        Path(location).unlink(missing_ok=True)

def simulate_spectrum(
    spinsys: str | object,
    delete_files: bool = True,
    filepath: str | None = None,
    **kwargs,
) -> object:
    """
    Simulate a spectrum from a spin system with smart defaults.

    Automatically calculates the spectral width and center frequency from the
    chemical shifts in the spin system, then runs a ``no_pulse`` simulation
    via SIMPSON.

    Parameters
    ----------
    spinsys : str or object
        Spin system definition. Can be a SIMPSON spinsys string or a Soprano
        SpinSystem object with a ``to_simpson()`` method.
    delete_files : bool
        If True (default), delete input/output files after simulation.
    filepath : str or None
        Path to save the SIMPSON input file. If None, uses a temporary file.
    **kwargs
        Override any default simulation parameter. Common overrides:
        ``proton_frequency``, ``spin_rate``, ``lb``, ``zerofill``,
        ``detect_operator``.

    Returns
    -------
    Simpy
        Object containing the simulated spectrum.
    """
    # Defaults
    defaults = {
        'proton_frequency': 800e6,
        'spin_rate': 30e3,
        'start_operator': 'Inx',
        'detect_operator': 'Inp',
        'crystal_file': 'rep2000',
        'gamma_angles': 16,
        'np': 4096,
        'method': 'direct',
        'verbose': 0,
        'out_format': 'spe',
        'lb': 100,
        'zerofill': 4096,
        'pulse_sequence': 'no_pulse'
    }

    params = defaults.copy()
    # Update with user kwargs
    user_detect_op = 'detect_operator' in kwargs
    params.update(kwargs)

    # Static spectra require gamma_angles=1
    if params.get('spin_rate', 0) == 0 and 'gamma_angles' not in kwargs:
        params['gamma_angles'] = 1

    # Extract spinsys string
    spinsys_str = ""
    if hasattr(spinsys, 'to_simpson'):
        spinsys_str = spinsys.to_simpson()
    else:
        spinsys_str = str(spinsys)

    # Detect nucleus and spin early — needed for both detect_operator selection and SW estimation
    b0 = _proton_freq_to_b0(params['proton_frequency'])
    nucleus = _extract_nucleus(spinsys_str) or '1H'
    spin = None
    try:
        spin = get_spin(nucleus)
    except ValueError:
        logger.debug("Could not determine spin for nucleus %s", nucleus)

    # Switch to Inc (central transition) for quadrupolar nuclei unless the user explicitly set detect_operator.
    if not user_detect_op and spin is not None and spin > 0.5:
        params['detect_operator'] = 'Inc'

    # Parse chemical shifts and quadrupolar couplings in one pass
    shifts = []
    quadrupoles = []
    for line in spinsys_str.splitlines():
        stripped = line.strip()
        if stripped.startswith('shift'):
            parts = stripped.split()
            if len(parts) > 2:
                val_str = parts[2]
                if val_str.endswith('p'):
                    shifts.append(float(val_str[:-1]))
                else:
                    with contextlib.suppress(ValueError):
                        shifts.append(float(val_str))
        elif stripped.startswith('quadrupole'):
            parts = stripped.split()
            # quadrupole site order Cq eta alpha beta gamma
            if len(parts) >= 4:
                with contextlib.suppress(ValueError):
                    quadrupoles.append(abs(float(parts[3])))

    spin_rate = params['spin_rate']

    if 'sw' not in kwargs:
        # get_larmor_freq() is signed (follows gamma). Spectral widths must use
        # the magnitude, or every negative-gamma nucleus gets a collapsed /
        # negative SW. The sign itself is still needed for the carrier offset.
        nu_l_signed = get_larmor_freq(b0, nucleus)
        gamma_sign = 1.0 if nu_l_signed >= 0 else -1.0
        nu_l_hz = abs(nu_l_signed) * 1e6
        is_ct = _is_ct_operator(params['detect_operator'])

        if shifts:
            min_shift = min(shifts)
            max_shift = max(shifts)
            center_ppm = (min_shift + max_shift) / 2
            center_hz = ppm2hz(center_ppm, b0, nucleus)
            width_hz = abs(ppm2hz(max_shift, b0, nucleus) - ppm2hz(min_shift, b0, nucleus))

            if quadrupoles and is_ct:
                # CT width dominated by second-order quadrupolar broadening
                max_cq = max(quadrupoles)
                required_sw = max(max_cq ** 2 / nu_l_hz, nu_l_hz * 10e-6)
            else:
                # Floor at 10 ppm so a single peak gets a usable window
                required_sw = max(width_hz * 2, nu_l_hz * 10e-6)

            # The carrier offset is a rotating-frame frequency and follows the
            # sign of gamma; the REF correction acts on the (sign-free) absolute
            # axis, so it is always -center_hz. Verified empirically against
            # SIMPSON for both gamma signs (13C and 29Si).
            if 'variable_offset' not in kwargs:
                params['variable_offset'] = gamma_sign * center_hz
            if 'variable_ref' not in kwargs:
                params['variable_ref'] = -center_hz

        elif quadrupoles:
            max_cq = max(quadrupoles)
            raw_sw = max_cq ** 2 / nu_l_hz if is_ct else 2.5 * max_cq
            required_sw = max(raw_sw, nu_l_hz * 10e-6)

        else:
            raise ValueError(
                "Cannot auto-estimate spectral width: no 'shift' or 'quadrupole' "
                "interactions found in spinsys. Provide 'sw' explicitly."
            )

        # Round up to the nearest multiple of spin_rate for MAS; use directly for static
        if spin_rate > 0:
            n = max(1, math.ceil(required_sw / spin_rate))
            sw_hz = n * spin_rate
        else:
            sw_hz = required_sw

        params['sw'] = sw_hz

    # Create calculator and run
    calc = SimpCalc(spinsys, **params)
    return calc.run(read_output=True, filepath=filepath, delete_files=delete_files)
