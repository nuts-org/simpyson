from __future__ import annotations

from pathlib import Path

import ase.io
import numpy as np
import scipy.constants as const

from simpyson.utils import get_larmor_freq


def _find_line(lines: list[str], marker: str) -> int:
    """Return the index of the last line containing *marker*, or raise."""
    idx = None
    for i, line in enumerate(lines):
        if marker in line:
            idx = i
    if idx is None:
        raise ValueError(
            f"Expected OUTCAR section '{marker}' not found. "
            "Is this a VASP NMR calculation OUTCAR?"
        )
    return idx


def read_vasp(file: str, format: str) -> ase.Atoms:
    """
    Read NMR data from a VASP OUTCAR file.

    Parses electric field gradients (EFG) and magnetic shielding tensors from
    a VASP NMR calculation and attaches them to an ASE Atoms object as arrays
    named ``'efg'`` and ``'ms'``.

    Parameters
    ----------
    file : str
        Path to the VASP OUTCAR file.
    format : str
        ASE format string for reading the structure (e.g., ``'vasp-out'``).

    Returns
    -------
    ase.Atoms
        Atoms object with ``'efg'`` and ``'ms'`` arrays attached.

    Raises
    ------
    ValueError
        If the OUTCAR is missing expected NMR sections.

    Examples
    --------
    >>> atoms = read_vasp('OUTCAR', 'vasp-out')
    >>> atoms.get_array('ms')  # magnetic shielding tensors
    """
    atoms = ase.io.read(file, format=format)
    n_atoms = atoms.get_global_number_of_atoms()
    np.set_printoptions(suppress=True)

    with Path(file).open() as outcar:
        lines = outcar.readlines()

    # Locate required OUTCAR sections
    idx_efg = _find_line(lines, "Electric field gradients (V/A^2)")
    idx_sym = _find_line(lines, "SYMMETRIZED TENSORS")
    idx_g0 = _find_line(lines, "G=0 CONTRIBUTION TO CHEMICAL SHIFT")
    idx_core = _find_line(lines, "Core NMR properties")

    # Magnetic susceptibility
    idx_sus = _find_line(lines, "Core contribution to magnetic susceptibility:")
    sus_parts = lines[idx_sus].split()
    mag_sus = float(sus_parts[5]) * 10 ** int(sus_parts[6][-2:])

    # Cell volume (first occurrence)
    volume = None
    for line in lines:
        if "volume of cell" in line:
            volume = float(line.split()[4])
            break
    if volume is None:
        raise ValueError("Cell volume not found in OUTCAR.")

    # Avogadro-based conversion factor for magnetic susceptibility
    chi_fact = 3.0 / 8.0 / np.pi * volume * 6.022142e23 / 1e24

    # --- EFG tensors ---
    efg = []
    for i in range(n_atoms):
        # EFG data starts 4 lines after the header: V_xx, V_yy, V_zz, V_xy, V_xz, V_yz
        grad = lines[idx_efg + 4 + i].split()[1:]
        matrix = np.array(
            [
                [grad[0], grad[3], grad[4]],  # xx xy xz
                [grad[3], grad[1], grad[5]],  # xy yy yz
                [grad[4], grad[5], grad[2]],  # xz yz zz
            ]
        )
        efg.append(matrix)

    # --- Symmetrized shielding tensors (3 rows per atom, every 4th line is a separator) ---
    sym_tensor = []
    for i in range(n_atoms * 4):
        if i % 4 != 0:
            sym_tensor.append(lines[idx_sym + i + 1].split())

    # --- G=0 constant shielding (3x3 tensor) ---
    # VASP 6.4.1+ has an extra description line
    if lines[idx_g0 + 1].strip() == "using pGv susceptibility, excluding core contribution":
        start_idx = idx_g0 + 5
    else:
        start_idx = idx_g0 + 4

    const_shield = []
    for i in range(3):
        const_shield.append(lines[start_idx + i].split()[1:])

    # --- Core shielding per element type ---
    unique_elements = np.unique(atoms.get_chemical_symbols())
    core_shield_rows = []
    for i in range(len(unique_elements)):
        core_shield_rows.append(lines[idx_core + 4 + i].split()[1:])
    core_shield_map = {row[0]: float(row[1]) for row in core_shield_rows}
    core_shield = np.array([core_shield_map[el] for el in atoms.get_chemical_symbols()])

    # --- Assemble results ---
    efg = np.array(efg, dtype=float)
    efg = efg * 1e20 / const.physical_constants["atomic unit of electric field gradient"][0]
    sym_tensor = np.split(np.array(sym_tensor, dtype=float), n_atoms)
    const_shield = np.array(const_shield, dtype=float)

    ms = []
    for i in range(n_atoms):
        core_diag = np.diag(core_shield[i] * np.ones(3))
        ms_tensor = (
            sym_tensor[i]
            + const_shield
            + core_diag
            + mag_sus / chi_fact * 1e6 * np.eye(3)
        )
        ms.append(-ms_tensor)
    ms = np.array(np.array(ms).tolist(), dtype=float)

    atoms.set_array('efg', efg)
    atoms.set_array('ms', ms)

    return atoms


def hz2ppm(
    hz: np.ndarray | float,
    b0: str,
    nucleus: str,
    isotope_file: str | None = None,
) -> np.ndarray | float:
    """
    Convert frequency values from Hz to ppm.

    Parameters
    ----------
    hz : numpy.ndarray or float
        Frequency values in Hz.
    b0 : str
        Magnetic field strength (e.g., ``'400MHz'`` or ``'9.4T'``).
    nucleus : str
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
    isotope_file : str or None
        Path to isotope data JSON file. If None, uses the bundled default.

    Returns
    -------
    numpy.ndarray or float
        Chemical shift values in ppm.

    Raises
    ------
    ValueError
        If the B0 unit is invalid or the nucleus is not found.
    """
    # SIMPSON places +delta at +delta*|nu_L| regardless of the sign of gamma,
    # so the axis conversion must use the Larmor frequency magnitude.
    # Using the signed value mirrored every negative-gamma (29Si, 15N, ...) axis.
    larmor_freq = abs(get_larmor_freq(b0, nucleus, isotope_file))
    return hz / larmor_freq


def ppm2hz(
    ppm: np.ndarray | float,
    b0: str,
    nucleus: str,
    isotope_file: str | None = None,
) -> np.ndarray | float:
    """
    Convert chemical shift values from ppm to Hz.

    Parameters
    ----------
    ppm : numpy.ndarray or float
        Chemical shift values in ppm.
    b0 : str
        Magnetic field strength (e.g., ``'400MHz'`` or ``'9.4T'``).
    nucleus : str
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
    isotope_file : str or None
        Path to isotope data JSON file. If None, uses the bundled default.

    Returns
    -------
    numpy.ndarray or float
        Frequency values in Hz.

    Raises
    ------
    ValueError
        If the B0 unit is invalid or the nucleus is not found.
    """
    # See hz2ppm: SIMPSON's frequency axis convention uses |nu_L|.
    larmor_freq = abs(get_larmor_freq(b0, nucleus, isotope_file))
    return ppm * larmor_freq
