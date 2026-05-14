from __future__ import annotations

import copy as cp
import logging
import warnings
from pathlib import Path

import numpy as np

from simpyson.converter import hz2ppm

logger = logging.getLogger("simpyson")


class Simpy:
    """
    Unified container for SIMPSON NMR data with automatic format conversions.

    Stores and manages SIMPSON data in multiple representations (FID, spectrum
    in Hz, spectrum in ppm) and converts between them lazily on access.

    Parameters
    ----------
    b0 : str or None
        Magnetic field strength (e.g., ``'9.4T'``, ``'400MHz'``).
        Required for ppm conversion.
    nucleus : str or None
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
        Required for ppm conversion.

    Examples
    --------
    >>> from simpyson.io import read_simp
    >>> data = read_simp("spectrum.spe", b0="400MHz", nucleus="13C")
    >>> data.spe['hz']   # frequency axis in Hz
    >>> data.ppm['ppm']  # chemical shift axis in ppm
    """

    def __init__(self, b0: str | None = None, nucleus: str | None = None) -> None:
        self._b0 = b0
        self._nucleus = nucleus
        self._fid_data: dict | None = None
        self._spe_data: dict | None = None
        self._xreim_data: dict | None = None
        self._metadata: dict = {}

    def __repr__(self) -> str:
        meta = []
        if self._b0:
            meta.append(f"b0={self._b0!r}")
        if self._nucleus:
            meta.append(f"nucleus={self._nucleus!r}")

        data = []
        if self._fid_data is not None:
            n, sw = int(self._fid_data['np']), self._fid_data['sw']
            data.append(f"fid(np={n}, sw={sw:.0f} Hz)")
        if self._spe_data is not None:
            n, sw = int(self._spe_data['np']), self._spe_data['sw']
            ppm_flag = " +ppm" if 'ppm' in self._spe_data else ""
            data.append(f"spe(np={n}, sw={sw:.0f} Hz{ppm_flag})")
        if self._xreim_data is not None:
            data.append(f"xreim(n={len(self._xreim_data['time'])})")

        payload = ", ".join(data) if data else "empty"
        prefix = (", ".join(meta) + ", ") if meta else ""
        return f"Simpy({prefix}{payload})"

    @property
    def b0(self) -> str | None:
        """Magnetic field strength."""
        return self._b0

    @b0.setter
    def b0(self, value: str | None) -> None:
        self._b0 = value
        # Invalidate cached ppm data when b0 changes
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']

    @property
    def nucleus(self) -> str | None:
        """Nucleus type."""
        return self._nucleus

    @nucleus.setter
    def nucleus(self, value: str | None) -> None:
        self._nucleus = value
        # Invalidate cached ppm data if nucleus changes
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']

    @property
    def fid(self) -> dict | None:
        """Access time-domain data, converting from spectrum if needed.

        Returns
        -------
        dict or None
            Dictionary with keys ``'real'``, ``'imag'``, ``'np'``, ``'sw'``,
            ``'time'`` (in ms), or None if no data is available.
        """
        if self._fid_data is None and self._spe_data is not None:
            self._compute_fid()
        return self._fid_data

    @property
    def spe(self) -> dict | None:
        """Access frequency-domain data, converting from FID if needed.

        Returns
        -------
        dict or None
            Dictionary with keys ``'real'``, ``'imag'``, ``'np'``, ``'sw'``,
            ``'hz'`` (and optionally ``'ppm'``), or None if no data is
            available.
        """
        if self._spe_data is None and self._fid_data is not None:
            self._compute_spectrum()
        return self._spe_data

    @property
    def ppm(self) -> dict | None:
        """Access chemical shift data, computing from Hz if needed.

        Requires both ``b0`` and ``nucleus`` to be set.

        Returns
        -------
        dict or None
            Dictionary with keys ``'ppm'``, ``'real'``, ``'imag'``, or None
            if conversion is not possible.
        """
        if self.spe and 'ppm' not in self.spe and self._b0 and self._nucleus:
            self._compute_ppm()

        if self.spe and 'ppm' in self.spe:
            return {
                'ppm': self.spe['ppm'],
                'real': self.spe['real'],
                'imag': self.spe['imag']
            }
        return None

    @property
    def xreim(self) -> dict | None:
        """Access xreim data.

        Returns
        -------
        dict or None
            Dictionary with keys ``'time'``, ``'real'``, ``'imag'``, or None.
        """
        return self._xreim_data

    def _compute_spectrum(self) -> None:
        """Convert FID to spectrum via FFT."""
        if not self._fid_data:
            return

        npoints = self._fid_data['np']
        sw = self._fid_data['sw']
        signal = self._fid_data['real'] + 1j * self._fid_data['imag']
        spectrum = np.fft.fftshift(np.fft.fft(signal))

        self._spe_data = {
            'real': np.real(spectrum),
            'imag': np.imag(spectrum),
            'np': npoints,
            'sw': sw,
            'hz': sw * (np.arange(int(npoints)) / int(npoints) - 0.5)
        }

    def _compute_fid(self) -> None:
        """Convert spectrum to FID via inverse FFT."""
        if not self._spe_data:
            return

        npoints = self._spe_data['np']
        sw = self._spe_data['sw']
        signal = self._spe_data['real'] + 1j * self._spe_data['imag']
        time_signal = np.fft.ifft(np.fft.ifftshift(signal))
        dt = 1.0 / sw

        self._fid_data = {
            'real': np.real(time_signal),
            'imag': np.imag(time_signal),
            'np': npoints,
            'sw': sw,
            'time': np.arange(int(npoints)) * dt * 1e3 # seconds to milliseconds
        }

    def _compute_ppm(self) -> None:
        """Calculate ppm scale from Hz using b0 and nucleus."""
        if not (self._spe_data and 'hz' in self._spe_data):
            return

        if not (self._b0 and self._nucleus):
            return

        try:
            self._spe_data['ppm'] = hz2ppm(self._spe_data['hz'], self._b0, self._nucleus)
        except ValueError as e:
            warnings.warn(f"Could not convert to ppm: {e}", stacklevel=2)

    def from_fid(
        self,
        real: np.ndarray | list,
        imag: np.ndarray | list,
        np_value: int | float,
        sw: float,
        time: np.ndarray | list | None = None,
    ) -> Simpy:
        """
        Populate this object with FID (time-domain) data.

        Parameters
        ----------
        real : array_like
            Real part of the FID signal.
        imag : array_like
            Imaginary part of the FID signal.
        np_value : int or float
            Number of data points.
        sw : float
            Spectral width in Hz.
        time : array_like or None
            Time axis in milliseconds. If None, calculated from ``sw`` and
            ``np_value``.

        Returns
        -------
        Simpy
            Self, for method chaining.
        """
        if time is None:
            dt = 1.0 / sw
            time = np.arange(int(np_value)) * dt * 1e3  # seconds to milliseconds

        self._fid_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'time': np.array(time)
        }

        # Clear cached spectrum data
        self._spe_data = None
        return self

    def from_spe(
        self,
        real: np.ndarray | list,
        imag: np.ndarray | list,
        np_value: int | float,
        sw: float,
        hz: np.ndarray | list | None = None,
    ) -> Simpy:
        """
        Populate this object with spectrum (frequency-domain) data.

        Parameters
        ----------
        real : array_like
            Real part of the spectrum.
        imag : array_like
            Imaginary part of the spectrum.
        np_value : int or float
            Number of data points.
        sw : float
            Spectral width in Hz.
        hz : array_like or None
            Frequency axis in Hz. If None, calculated as a symmetric range
            around zero based on ``sw``.

        Returns
        -------
        Simpy
            Self, for method chaining.
        """
        if hz is None:
            hz = sw * (np.arange(int(np_value)) / int(np_value) - 0.5)

        self._spe_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'hz': np.array(hz)
        }

        # Calculate ppm if possible
        if self._b0 and self._nucleus:
            self._compute_ppm()

        # Clear cached FID data
        self._fid_data = None
        return self

    def from_xreim(
        self,
        time: np.ndarray | list,
        real: np.ndarray | list,
        imag: np.ndarray | list,
    ) -> Simpy:
        """
        Populate this object with xreim data.

        Parameters
        ----------
        time : array_like
            Time axis.
        real : array_like
            Real part of the signal.
        imag : array_like
            Imaginary part of the signal.

        Returns
        -------
        Simpy
            Self, for method chaining.
        """
        self._xreim_data = {
            'time': np.array(time),
            'real': np.array(real),
            'imag': np.array(imag)
        }

        # Clear cached FID and spectrum data
        self._fid_data = None
        self._spe_data = None
        return self

    def from_csdf(
        self,
        real: np.ndarray | list,
        imag: np.ndarray | list,
        hz: np.ndarray | list,
        np_value: int | float,
        sw: float,
    ) -> Simpy:
        """
        Populate this object from CSDF (csdmpy) data.

        Parameters
        ----------
        real : array_like
            Real part of the spectrum.
        imag : array_like
            Imaginary part of the spectrum.
        hz : array_like
            Frequency axis in Hz.
        np_value : int or float
            Number of data points.
        sw : float
            Spectral width in Hz.

        Returns
        -------
        Simpy
            Self, for method chaining.
        """
        self._spe_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'hz': np.array(hz)
        }

        if self._b0 and self._nucleus:
            self._compute_ppm()

        self._fid_data = None
        return self

    def copy(self) -> Simpy:
        """
        Create a deep copy of this Simpy object.

        Returns
        -------
        Simpy
            Independent copy with all data arrays duplicated.
        """
        new_obj = Simpy(b0=self._b0, nucleus=self._nucleus)

        if self._fid_data is not None:
            new_obj._fid_data = cp.deepcopy(self._fid_data)

        if self._spe_data is not None:
            new_obj._spe_data = cp.deepcopy(self._spe_data)

        if self._xreim_data is not None:
            new_obj._xreim_data = cp.deepcopy(self._xreim_data)

        new_obj._metadata = cp.deepcopy(self._metadata)

        return new_obj

    def write(self, filename: str, format: str = 'csv') -> Simpy:
        """
        Write data to file in the specified format.

        Parameters
        ----------
        filename : str
            Output file path.
        format : str
            Output format. One of ``'csv'``, ``'spe'``, ``'fid'``, ``'xreim'``.

        Returns
        -------
        Simpy
            Self, for method chaining.

        Raises
        ------
        ValueError
            If the format is unsupported or no data is available.
        """
        # CSV format
        if format == 'csv':
            if self.spe:
                if 'ppm' in self.spe and self._b0 and self._nucleus:
                    x_data, x_label = self.spe['ppm'], 'ppm'
                else:
                    x_data, x_label = self.spe['hz'], 'Hz'
                y_data = self.spe['real']
            elif self.fid:
                x_data, x_label = self.fid['time'], 'Time (ms)'
                y_data = self.fid['real']
            elif self.xreim:
                x_data, x_label = self.xreim['time'], 'Time (ms)'
                y_data = self.xreim['real']
            else:
                raise ValueError("No data to save")

            np.savetxt(
                filename,
                np.column_stack((x_data, y_data)),
                delimiter=",",
                header=f"{x_label},Real",
                comments=""
            )
            return self

        # SIMPSON formats
        data_dict = None
        if format == 'spe':
            data_dict = self.spe
            data_type = 'SPE'
        elif format == 'fid':
            data_dict = self.fid
            data_type = 'FID'
        elif format == 'xreim':
            data_dict = self.xreim
            data_type = 'XREIM'
        else:
            raise ValueError(f"Unsupported save format: {format}")

        if not data_dict:
            raise ValueError(f"No data available to save in {format} format")

        # Write SIMPSON format file
        with Path(filename).open('w') as f:
            f.write('SIMP\n')
            if 'np' in data_dict:
                f.write(f'NP={data_dict["np"]}\n')
            if 'sw' in data_dict:
                f.write(f'SW={data_dict["sw"]}\n')
            f.write(f'TYPE={data_type}\n')
            f.write('DATA\n')

            for re_val, im_val in zip(data_dict['real'], data_dict['imag'], strict=False):
                f.write(f'{re_val} {im_val}\n')

            f.write('END')

        return self
