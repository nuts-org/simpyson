from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path

import numpy as np

from simpyson.simpy import Simpy

logger = logging.getLogger("simpyson")


def read_spe(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON SPE file.

    Parameters
    ----------
    filename : str
        Path to the ``.spe`` file.
    simpy_data : Simpy
        Object to populate with spectrum data.

    Raises
    ------
    ValueError
        If required header fields (NP, SW) are missing.
    """
    with Path(filename).open() as f:
        data_sec = False
        real: list[float] = []
        imag: list[float] = []
        ref = 0.0
        np_value: float | None = None
        sw: float | None = None
        for line in f:
            if line.startswith('NP'):
                np_value = float(line.split('=')[1])
            elif line.startswith('SW'):
                sw = float(line.split('=')[1])
            elif line.startswith('REF'):
                ref = float(line.split('=')[1])
            elif line.startswith('DATA'):
                data_sec = True
            elif data_sec and line.startswith('END'):
                break
            elif data_sec:
                a, b = map(float, line.split())
                real.append(a)
                imag.append(b)

        if np_value is None or sw is None:
            raise ValueError(
                f"Missing required header fields in {filename}: "
                f"{'NP' if np_value is None else ''}"
                f"{' and ' if np_value is None and sw is None else ''}"
                f"{'SW' if sw is None else ''} not found."
            )

        if ref != 0.0:
            logger.debug("Applying REF=%g Hz from SPE header: hz = f_SPE - REF", ref)
        indices = np.arange(np_value)
        hz = sw * (indices / np_value - 0.5) - ref

        simpy_data.from_spe(real, imag, np_value, sw, hz)


def read_fid(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON FID file.

    Parameters
    ----------
    filename : str
        Path to the ``.fid`` file.
    simpy_data : Simpy
        Object to populate with FID data.

    Raises
    ------
    ValueError
        If required header fields (NP, SW) are missing.
    """
    with Path(filename).open() as f:
        data_sec = False
        real: list[float] = []
        imag: list[float] = []
        np_value: float | None = None
        sw: float | None = None
        for line in f:
            if line.startswith('NP'):
                np_value = float(line.split('=')[1])
            elif line.startswith('SW'):
                sw = float(line.split('=')[1])
            elif line.startswith('DATA'):
                data_sec = True
            elif data_sec and line.startswith('END'):
                break
            elif data_sec:
                a, b = map(float, line.split())
                real.append(a)
                imag.append(b)

        if np_value is None or sw is None:
            raise ValueError(
                f"Missing required header fields in {filename}: "
                f"{'NP' if np_value is None else ''}"
                f"{' and ' if np_value is None and sw is None else ''}"
                f"{'SW' if sw is None else ''} not found."
            )

        # Let from_fid() compute the time axis (avoids duplicating the calculation)
        simpy_data.from_fid(np.array(real), np.array(imag), np_value, sw)


def read_xreim(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON file saved with the ``-xreim`` option.

    Parameters
    ----------
    filename : str
        Path to the ``.xreim`` file.
    simpy_data : Simpy
        Object to populate with xreim data.
    """
    with Path(filename).open() as f:
        time: list[float] = []
        real: list[float] = []
        imag: list[float] = []
        for line in f:
            parts = line.split()
            time.append(float(parts[0]))
            real.append(float(parts[1]))
            imag.append(float(parts[2]))

        simpy_data.from_xreim(np.array(time), np.array(real), np.array(imag))


def read_csdf(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a CSDF file (CSDM format).

    Requires the optional ``csdmpy`` dependency. The frequency axis is
    always converted to Hz; files storing coordinates in other units (e.g.
    kHz) are handled automatically via unit conversion.

    Parameters
    ----------
    filename : str
        Path to the ``.csdf`` file.
    simpy_data : Simpy
        Object to populate with spectrum data.
    """
    import csdmpy as csdm  # noqa: PLC0415

    data = csdm.load(filename)
    hz = data.dimensions[0].coordinates.to('Hz').value
    real = data.dependent_variables[0].components[0].real
    imag = data.dependent_variables[0].components[0].imag
    np_value = len(hz)
    sw = float(np.abs(hz[-1] - hz[0]))

    simpy_data.from_csdf(real, imag, hz, np_value, sw)


# Defined once, after the reader functions below
_EXT_TO_FMT: dict[str, str] = {
    '.spe':   'spe',
    '.fid':   'fid',
    '.xreim': 'xreim',
    '.csdf':  'csdf',
}

_READERS: dict[str, Callable[[str, Simpy], None]] = {
    'spe':   read_spe,
    'fid':   read_fid,
    'xreim': read_xreim,
    'csdf':  read_csdf,
}


def read_simp(
    filename: str,
    fmt: str | None = None,
    b0: str | None = None,
    nucleus: str | None = None,
) -> Simpy:
    """
    Read SIMPSON NMR data from a file into a unified Simpy object.

    The file format is determined from the extension if ``fmt`` is not
    given explicitly.

    Parameters
    ----------
    filename : str
        Path to the SIMPSON output file.
    fmt : str or None
        File format (``'spe'``, ``'fid'``, ``'xreim'``, ``'csdf'``).
        If None, guessed from the file extension.
    b0 : str or None
        Magnetic field strength (e.g., ``'9.4T'``, ``'400MHz'``).
        Needed for ppm conversion.
    nucleus : str or None
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
        Needed for ppm conversion.

    Returns
    -------
    Simpy
        Object containing the loaded data.

    Raises
    ------
    ValueError
        If the file format cannot be determined or is unsupported.
    OSError
        If the file cannot be read or parsed.
    """
    if fmt is not None:
        fmt = fmt.lower()
    else:
        ext = Path(filename).suffix.lower()
        fmt = _EXT_TO_FMT.get(ext)
        if fmt is None:
            raise ValueError(
                f"Cannot determine file format of {filename!r}. "
                f"Supported extensions: {sorted(_EXT_TO_FMT)}"
            )

    reader = _READERS.get(fmt)
    if reader is None:
        raise ValueError(
            f"Unsupported format {fmt!r}. Supported: {sorted(_READERS)}"
        )

    simpy_data = Simpy(b0=b0, nucleus=nucleus)
    try:
        reader(filename, simpy_data)
    except (ValueError, KeyError, IndexError, OSError) as e:
        raise OSError(f"Error reading {filename!r} as {fmt!r}: {e}") from e
    return simpy_data
