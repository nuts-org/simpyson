import numpy as np
import os
import copy
from simpyson.converter import hz2ppm, ppm2hz
from simpyson.utils import get_larmor_freq

class Simpy:
    """
    A unified container Simpson NMR data that automatically handles conversions between formats.
    
    This class stores and manages Simpson's various formasts (FID, Spe in Hz and ppm) and 
    conversion between them.
    
    Attributes:
        b0: Magnetic field (e.g., '9.4T', '400MHz')
        nucleus: Nucleus type (e.g., '1H', '13C')
        
    Properties:
        fid: Time-domain data
        spe: Frequency-domain data (Hz)
        ppm: Chemical shift data (ppm)
    """
    def __init__(self, b0=None, nucleus=None):
        self._b0 = b0
        self._nucleus = nucleus
        self._fid_data = None  # Time-domain data
        self._spe_data = None  # Frequency-domain data
        self._xreim_data = None # xreim data
        self._metadata = {}    # Additional information
    
    @property
    def b0(self):
        return self._b0
    
    @b0.setter
    def b0(self, value):
        self._b0 = value
        # Invalidate cached ppm data when b0 change
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']
    
    @property
    def nucleus(self):
        return self._nucleus
    
    @nucleus.setter
    def nucleus(self, value):
        self._nucleus = value
        # Invalidate cached ppm data if nucleus changs
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']
    
    @property
    def fid(self):
        """Access time-domain data, converting from spectrum if needed."""
        if self._fid_data is None and self._spe_data is not None:
            self._compute_fid()
        return self._fid_data
    
    @property
    def spe(self):
        """Access frequency-domain data, converting from FID if needed."""
        if self._spe_data is None and self._fid_data is not None:
            self._compute_spectrum()
        return self._spe_data
    
    @property
    def ppm(self):
        """Access chemical shift data, computing from Hz if needed."""
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
    def xreim(self):
        """Access xreim data."""
        return self._xreim_data

    def _compute_spectrum(self):
        """Convert FID to spectrum."""
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
            'hz': np.linspace(-sw/2, sw/2, int(npoints))
        }
    
    def _compute_fid(self):
        """Convert spectrum to FID."""
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
            'time': np.linspace(0, npoints*dt, int(npoints)) * 10e3
        }
    
    def _compute_ppm(self):
        """Calculate ppm scale from Hz."""
        if not (self._spe_data and 'hz' in self._spe_data):
            return
        
        if not (self._b0 and self._nucleus):
            return
            
        try:
            self._spe_data['ppm'] = hz2ppm(self._spe_data['hz'], self._b0, self._nucleus)
        except ValueError as e:
            print(f"Error converting to ppm: {e}")
    
    def from_fid(self, real, imag, np_value, sw, time=None):
        """Set data from FID values."""
        if time is None:
            dt = 1.0 / sw
            time = np.linspace(0, np_value*dt, int(np_value)) * 10e3
            
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
    
    def from_spe(self, real, imag, np_value, sw, hz=None):
        """Set data from spectrum values."""
        if hz is None:
            hz = np.linspace(-int(sw) / 2, int(sw) / 2, int(np_value))
            
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
    
    def from_xreim(self, time, real, imag):
        """Set data from xreim values."""

        self._xreim_data = {
            'time': np.array(time),
            'real': np.array(real),
            'imag': np.array(imag)
        }

        # Clear cached FID and spectrum data
        self._fid_data = None
        self._spe_data = None
        return self
    
    def copy(self):
        """Create a copy of a Simpy object."""
        import copy as cp
    
        new_obj = Simpy(b0=self._b0, nucleus=self._nucleus)

        if self._fid_data is not None:
            new_obj._fid_data = cp.deepcopy(self._fid_data)

        if self._spe_data is not None:
            new_obj._spe_data = cp.deepcopy(self._spe_data)

        if self._xreim_data is not None:
            new_obj._xreim_data = cp.deepcopy(self._xreim_data)

        new_obj._metadata = cp.deepcopy(self._metadata)

        return new_obj
    
    def write(self, filename, format='csv'):
        """
        Write data to file in specified format.
        
        Args:
            filename: Output filename
            format: Format to save as ('csv', 'fid', 'spe')
            
        Returns:
            self for method chaining
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
        with open(filename, 'w') as f:
            f.write('SIMP\n')
            if 'np' in data_dict:
                f.write(f'NP={data_dict["np"]}\n')
            if 'sw' in data_dict:
                f.write(f'SW={data_dict["sw"]}\n')
            f.write(f'TYPE={data_type}\n')
            f.write('DATA\n')
            
            for re, im in zip(data_dict['real'], data_dict['imag']):
                f.write(f'{re} {im}\n')
                
            f.write('END')
            
        return self
    
# Simpson calculator
class Simpcalc:
    """
    Class to create a SIMPSON simulation input.

    Attributes:
        spinsys (str): Spin system from Soprano or a custom string.
        out_name (str): Output file name.
        out_format (str): Output format (fid, spe, xreim).
        spin_rate (float): Spin rate in Hz.
        np (int): Number of points.
        proton_freq (float): Proton frequency in Hz.
        start_op (str): Start operator.
        detect_op (str): Detect operator.
        crystal_file (str): Crystal file.
        gamma_angles (str): Gamma angles.
        sw (float): Spectral width.
        verbose (bool): Verbose output.
        lb (float): Line broadening.
        zerofill (int): Zero filling.
        method (str): Method of simulation (direct, indirect, ...).
        tsw (str, optional): Spectral width in time domain.
        pulse_sequence (str, optional): Pulse sequence from templates or a custom string.
        pH (float, optional): pulse for H in us.
        pX (float, optional): pulse for X in us.
        pY (float, optional): pulse for Y in us.
        plH (float, optional): power level for H in Hz.
        plX (float, optional): power level for X in Hz.
        plY (float, optional): power level for Y in Hz.
        phH (float, optional): phase for H pH.
        phX (float, optional): phase for X pH.
        phY (float, optional): phase for Y pH.

    Example:
        sim = SimpSim(spinsys=spinsys, out_name='output', out_format='spe', spin_rate=15e3, np=2048, proton_freq=400e6, start_op='Inx', detect_op='Inp', crystal_file='rep100', gamma_angles=4, sw=20e3, verbose=0, lb=20, zerofill=4096)
    """
    def __init__(self, spinsys, out_name, out_format, spin_rate, np, proton_freq, start_op, detect_op, crystal_file, gamma_angles, sw, verbose, lb, zerofill, method="direct", tsw=None, pulse_sequence=None, pH=None, pX=None, pY=None, plH=None, plX=None, plY=None, phH=None, phX=None, phY=None):
        self.spin_rate = spin_rate
        self.spinsys = spinsys
        self.out_name = out_name
        self.out_format = out_format
        self.np = np
        self.proton_freq = proton_freq
        self.start_op = start_op
        self.detect_op = detect_op
        self.crystal_file = crystal_file
        self.gamma_angles = gamma_angles
        self.sw = sw
        self.verbose = verbose
        self.tsw = tsw if tsw else f"1e6/{sw}"
        self.lb = lb
        self.zerofill = zerofill
        self.method = method
        self.pulse_sequence = pulse_sequence
        self.pH = pH
        self.pX = pX
        self.pY = pY
        self.plH = plH
        self.plX = plX
        self.plY = plY
        self.phH = phH
        self.phX = phX
        self.phY = phY

    def par_content(self):
        par_block = f"""
par {{
    spin_rate        {self.spin_rate}
    np               {self.np}
    proton_frequency {self.proton_freq}
    start_operator   {self.start_op}
    detect_operator  {self.detect_op}
    method           {self.method}
    crystal_file     {self.crystal_file}
    gamma_angles     {self.gamma_angles}
    variable sw      {self.sw}
    verbose          {self.verbose}
    variable tsw     {self.tsw}
"""
        if self.pH is not None:
            if self.plH is None or self.phH is None:
                raise ValueError("plH and phH must be defined if pH is defined")
            par_block += f"    pH               {self.pH}\n"
            par_block += f"    plH              {self.plH}\n"
            par_block += f"    phH              {self.phH}\n"

        if self.pX is not None:
            if self.plX is None or self.phX is None:
                raise ValueError("plX and phX must be defined if pX is defined")
            par_block += f"    pX               {self.pX}\n"
            par_block += f"    plX              {self.plX}\n"
            par_block += f"    phX              {self.phX}\n"

        if self.pY is not None:
            if self.plY is None or self.phY is None:
                raise ValueError("plY and phY must be defined if pY is defined")
            par_block += f"    pY               {self.pY}\n"
            par_block += f"    plY              {self.plY}\n"
            par_block += f"    phY              {self.phY}\n"

        par_block += "}\n"
        return par_block

    def pulseq_content(self):
        if self.pulse_sequence:
            return self.pulse_sequence
        else:
            return no_pulse

    def main_content(self):
        if self.out_format == "fid":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    faddlb $f {self.lb} 0
    fzerofill $f {self.zerofill}
    fsave $f {self.out_name}.fid
}}
"""
        elif self.out_format == "spe":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    faddlb $f {self.lb} 0
    fzerofill $f {self.zerofill}
    fft $f
    fsave $f {self.out_name}.spe
}}
"""
        elif self.out_format == "xreim":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    fsave $f {self.out_name}.xreim -xreim
}}
"""
        else:
            raise ValueError(f"Unknown out_format {self.out_format}")

    def __str__(self):
        return f"{self.spinsys}\n{self.par_content()}{self.pulseq_content()}{self.main_content()}"
    
    def save(self, filepath):
        with open(filepath, 'w') as file:
            file.write(str(self))
        return self