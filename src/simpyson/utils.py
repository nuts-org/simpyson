from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from soprano.calculate.nmr.simpson import _header_template, _spinsys_template
from soprano.data.nmr import _get_isotope_list
from soprano.properties.nmr.dipolar import DipolarCoupling
from soprano.selection import AtomSelection


def _default_isotope_file() -> str:
    """Return the path to the bundled isotope data JSON file."""
    return str(Path(__file__).parent / 'isotope_data.json')


def _load_isotope_data(nucleus: str, isotope_file: str | None = None) -> dict:
    """
    Parse a nucleus string and load its isotope data row.

    Parameters
    ----------
    nucleus : str
        Nucleus type (e.g., ``'1H'``, ``'13C'``, ``'23Na'``).
    isotope_file : str or None
        Path to isotope data JSON file. If None, uses the bundled default.

    Returns
    -------
    dict
        The isotope data row (contains keys like ``'Gamma'``, ``'Spin'``).

    Raises
    ------
    ValueError
        If the nucleus is not found in the isotope data.
    """
    if isotope_file is None:
        isotope_file = _default_isotope_file()

    mass_number = int(''.join(filter(str.isdigit, nucleus)))
    element = ''.join(filter(str.isalpha, nucleus)).capitalize()

    with Path(isotope_file).open() as f:
        data = json.load(f)

    if element in data and str(mass_number) in data[element]:
        return data[element][str(mass_number)]

    raise ValueError(f'Nucleus {nucleus} not found in isotope data.')


def get_gamma(nucleus: str, isotope_file: str | None = None) -> float:
    """
    Get the gyromagnetic ratio for a given nucleus.

    Parameters
    ----------
    nucleus : str
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
    isotope_file : str or None
        Path to isotope data JSON file. If None, uses the bundled default.

    Returns
    -------
    float
        Gyromagnetic ratio in rad/(s*T) * 1e-7 (the value stored in the
        isotope data file).

    Raises
    ------
    ValueError
        If the nucleus is not found in the isotope data.
    """
    return _load_isotope_data(nucleus, isotope_file)['Gamma']


def get_spin(nucleus: str, isotope_file: str | None = None) -> float:
    """
    Get the spin quantum number for a given nucleus.

    Parameters
    ----------
    nucleus : str
        Nucleus type (e.g., ``'1H'``, ``'23Na'``).
    isotope_file : str or None
        Path to isotope data JSON file. If None, uses the bundled default.

    Returns
    -------
    float
        Spin quantum number (e.g., 0.5 for spin-1/2, 1.5 for spin-3/2).

    Raises
    ------
    ValueError
        If the nucleus is not found in the isotope data.
    """
    spin_str = _load_isotope_data(nucleus, isotope_file)['Spin']
    if '/' in str(spin_str):
        num, den = spin_str.split('/')
        return float(num) / float(den)
    return float(spin_str)


def get_larmor_freq(b0: str, nucleus: str, isotope_file: str | None = None) -> float:
    """
    Calculate the Larmor frequency for a given nucleus at a given field strength.

    Parameters
    ----------
    b0 : str
        Magnetic field strength (e.g., '400MHz' or '9.4T').
    nucleus : str
        Nucleus type (e.g., '1H' or '13C').
    isotope_file : str, optional
        Path to isotope data file. If None, uses the bundled default.

    Returns
    -------
    float
        Larmor frequency in MHz.

    Raises
    ------
    ValueError
        If B0 unit is not 'T' or 'MHz', or if the nucleus is not found.
    """
    gamma = get_gamma(nucleus, isotope_file=isotope_file)
    b0_unit = ''.join(filter(str.isalpha, b0)).lower()
    b0_value = float(''.join(filter(lambda x: x.isdigit() or x == '.', b0)))

    # gamma is stored as gamma / 1e7 in the isotope data file,
    # so gamma * 1e7 gives the true value in rad/(s*T).
    if b0_unit == 't':
        larmor_freq = gamma * 1e7 * b0_value / (2 * np.pi * 1e6)
    elif b0_unit == 'mhz':
        gamma_h = get_gamma('1H', isotope_file=isotope_file)
        b0_tesla = 2 * np.pi * b0_value * 1e6 / (gamma_h * 1e7)
        larmor_freq = gamma * 1e7 * b0_tesla / (2 * np.pi * 1e6)
    else:
        raise ValueError('B0 unit must be T or MHz.')

    return larmor_freq

def add_spectra(spectra_list: list, b0: str | None = None, nucleus: str | None = None):
    """
    Combine multiple Simpy objects into a single spectrum by summing.

    Parameters
    ----------
    spectra_list : list of Simpy
        Spectra to combine. All must have compatible frequency-domain data.
    b0 : str, optional
        Magnetic field override (e.g., '400MHz').
    nucleus : str, optional
        Nucleus override (e.g., '1H').

    Returns
    -------
    Simpy or None
        Combined spectrum, or None if the input list is empty.

    Raises
    ------
    ValueError
        If any spectrum in the list has no frequency-domain data.
    """
    if not spectra_list:
        return None

    # Validate spectra have frequency-domain data
    spe_data_list = []
    for i, spectrum in enumerate(spectra_list, start=1):
        spe = spectrum.spe
        if spe is None:
            raise ValueError(f"Spectrum #{i} has no frequency-domain data to combine.")
        spe_data_list.append(spe)

    result = spectra_list[0].copy()

    if b0:
        result.b0 = b0
    if nucleus:
        result.nucleus = nucleus

    if len(spectra_list) == 1:
        return result

    # Check whether all spectra share the same Hz axis
    ref_hz = spe_data_list[0]['hz']
    axes_match = all(
        len(spe['hz']) == len(ref_hz) and np.allclose(spe['hz'], ref_hz)
        for spe in spe_data_list[1:]
    )

    if axes_match:
        # Fast path: identical Hz axes, just sum element-wise
        for spe in spe_data_list[1:]:
            result.spe['real'] += spe['real']
            result.spe['imag'] += spe['imag']
    else:
        # Interpolate all spectra onto a common Hz grid
        hz_min = min(spe['hz'][0] for spe in spe_data_list)
        hz_max = max(spe['hz'][-1] for spe in spe_data_list)

        # Use the smallest Hz step among all spectra
        steps = [spe['sw'] / spe['np'] for spe in spe_data_list]
        finest_step = min(steps)

        n_common = int(np.ceil((hz_max - hz_min) / finest_step)) + 1
        common_hz = np.linspace(hz_min, hz_max, n_common)

        # Interpolate and sum all spectra on common grid
        sum_real = np.zeros(n_common)
        sum_imag = np.zeros(n_common)

        for spe in spe_data_list:
            sum_real += np.interp(common_hz, spe['hz'], spe['real'],
                                  left=0.0, right=0.0)
            sum_imag += np.interp(common_hz, spe['hz'], spe['imag'],
                                  left=0.0, right=0.0)

        # Full spectral width is N * step, not the coordinate span (N-1) * step
        common_sw = (common_hz[1] - common_hz[0]) * n_common
        result.from_spe(sum_real, sum_imag, n_common, common_sw, common_hz)

    if result.b0 and result.nucleus:
        result._compute_ppm()

    return result

def simple_spinsys(
    atoms: object,
    isotopes: dict,
    iso_ms: list | None = None,
    aniso_ms: list | None = None,
    eta_ms: list | None = None,
    euler_ms: list | None = None,
    cq: list | None = None,
    eta_q: list | None = None,
    q_order: list | None = None,
    euler_q: list | None = None,
    get_dipolar: bool = False,
    dip_sel: object | None = None,
    obs_nuc: str | None = None,
) -> str:
    """
    Build a SIMPSON spinsys string from manually provided NMR parameters.

    Use this when you don't have DFT-derived data and want to construct a
    spin system by hand (e.g., for quick tests or educational purposes).
    For DFT-derived spin systems, use Soprano instead.

    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object defining the molecular structure.
    isotopes : dict
        Mapping of element symbol to mass number (e.g., ``{'H': 1, 'C': 13}``).
    iso_ms : list of float or None
        Isotropic chemical shifts in ppm, one per atom.
    aniso_ms : list of float or None
        Shielding anisotropies in ppm. Defaults to zeros.
    eta_ms : list of float or None
        Shielding asymmetry parameters. Defaults to zeros.
    euler_ms : list of array_like or None
        Euler angles ``(alpha, beta, gamma)`` in degrees for the shielding
        tensor. Defaults to zeros.
    cq : list of float or None
        Quadrupolar coupling constants in Hz.
    eta_q : list of float or None
        Quadrupolar asymmetry parameters. Defaults to zeros.
    q_order : list of int or None
        Quadrupolar perturbation orders (must be <= 2). Defaults to 2.
    euler_q : list of array_like or None
        Euler angles ``(alpha, beta, gamma)`` in degrees for the EFG tensor.
        Defaults to zeros.
    get_dipolar : bool
        If True, calculate and include dipolar couplings.
    dip_sel : AtomSelection or None
        Selection of atoms for dipolar couplings. Defaults to all atoms.
    obs_nuc : str or None
        Observed nucleus (e.g., ``'13C'``). Placed first in the channels list.

    Returns
    -------
    str
        Complete SIMPSON spinsys block as a string.

    Raises
    ------
    ValueError
        If a quadrupolar order exceeds 2.
    """

    # Header Block
    symbols = atoms.get_chemical_symbols()
    isotope_list = _get_isotope_list(symbols, isotopes=isotopes)
    nuclei = [f"{iso}{el}" for el, iso in zip(symbols, isotope_list, strict=False)]

    if obs_nuc and obs_nuc in nuclei:
        channels = [obs_nuc] + [n for n in sorted(set(nuclei)) if n != obs_nuc]
    else:
        channels = sorted(set(nuclei))

    header = _header_template.format(
        channels=" ".join(channels), nuclei=" ".join(nuclei)
    )

    # Magnetic shielding Block
    ms_block = ""
    if iso_ms is not None:
        num_atoms = len(atoms)
        aniso_ms_vals = aniso_ms if aniso_ms is not None else np.zeros(num_atoms)
        eta_ms_vals = eta_ms if eta_ms is not None else np.zeros(num_atoms)
        euler_ms_vals = euler_ms if euler_ms is not None else np.zeros((num_atoms, 3))

        for i in range(num_atoms):
            ms_block += "shift {} {}p {}p {} {} {} {}\n".format(
                i + 1,
                iso_ms[i],
                aniso_ms_vals[i],
                eta_ms_vals[i],
                *euler_ms_vals[i],
            )

    # Quadrupolar Block
    efg_block = ""
    if cq is not None:
        num_atoms = len(atoms)
        eta_q_vals = eta_q if eta_q is not None else np.zeros(num_atoms)
        euler_q_vals = euler_q if euler_q is not None else np.zeros((num_atoms, 3))

        for i in range(num_atoms):
            if cq[i] != 0:
                if q_order is not None:
                    if q_order[i] > 2:
                        raise ValueError(
                            f"Quadrupolar order must be 2 or less, got {q_order[i]} for atom {i+1}"
                        )
                    efg_block += "quadrupole {} {} {} {} {} {} {}\n".format(
                        i + 1, q_order[i], cq[i], eta_q_vals[i], *euler_q_vals[i]
                    )
                else:
                    efg_block += "quadrupole {} 2 {} {} {} {} {}\n".format(
                        i + 1, cq[i], eta_q_vals[i], *euler_q_vals[i]
                    )

    # Dipolar Block
    dip_block = ""
    if get_dipolar:
        if dip_sel is None:
            dip_sel = AtomSelection.all(atoms)

        if len(dip_sel) > 1:
            dip_couplings = DipolarCoupling.get(
                atoms, sel_i=dip_sel, isotope_list=isotope_list
            )
            for (i, j), (d, v) in dip_couplings.items():
                # Convert units
                d_rad_s = d * 2 * np.pi

                if np.allclose(v, [0, 0, 1]):
                    beta, alpha = 0.0, 0.0
                elif np.allclose(v, [0, 0, -1]):
                    beta, alpha = 180.0, 0.0
                else:
                    beta = np.arccos(v[2]) * 180 / np.pi
                    alpha = np.arctan2(v[1], v[0]) * 180 / np.pi

                dip_block += f"dipole {i + 1} {j + 1} {d_rad_s} {beta} {alpha} 0.0\n"

    return _spinsys_template.format(
        header=header, ms=ms_block, efg=efg_block, dipolar=dip_block
    )
