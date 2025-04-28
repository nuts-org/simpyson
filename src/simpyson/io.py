import numpy as np
import os
from simpyson.simpy import Simpy, Simpcalc

# Functions for reading Simpson output files
def read_simp(
        filename,
        format=None,
        b0=None,
        nucleus=None,
):
    """"
    Reads Simpson's NMR data from a file into a unified Simpy object.
    
    Args:
        filename: Path to the file
        format: File format ('spe', 'fid', 'xreim') or None to guess from extension
        b0: Magnetic field (e.g., '9.4T', '400MHz')
        nucleus: Nucleus (e.g., '1H', '13C')
        
    Returns:
        Simpy object with time- and frequency-domain data.
        """
    # Try to guess format
    ext = os.path.splitext(filename)[1].lower()
    if ext == '.spe':
        format = 'spe'
    elif ext == '.fid':
        format = 'fid'
    elif ext == '.xreim':
        format = 'xreim'
    elif ext == '.csdf':
        format = 'csdf'
    else:
        raise ValueError(f"Cannot determine file format of {filename}")
    
    simpy_data = Simpy(b0=b0, nucleus=nucleus)

    if format == 'spe':
        read_spe(filename, simpy_data)
    elif format == 'fid':
        read_fid(filename, simpy_data)
    elif format == 'xreim':
        read_xreim(filename, simpy_data)
    elif format == 'csdf':
        read_csdf(filename, simpy_data)
    else:
        raise ValueError(f"Unsupported format {format}")
    return simpy_data

# Define reading functions for each format
def read_spe(filename, simpy_data):
    """
    This method reads NMR data from a SIMPSON SPE file.
    """
    with open(filename) as f:
        data_sec = False
        real = []
        imag = []
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

        hz = np.linspace(-int(sw) / 2, int(sw) / 2, int(np_value))

        simpy_data.from_spe(real, imag, np_value, sw, hz)


def read_fid(filename, simpy_data):
    """
    This method reads NMR data from a SIMPSON FID file.
    """
    with open(filename) as f:
        data_sec = False
        real = []
        imag = []
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

        dt = 1.0 / sw
        time = np.linspace(0, np_value*dt, int(np_value))
        real = np.array(real)
        imag = np.array(imag)
        time = np.array(time)*10e3
        
        simpy_data.from_fid(real, imag, np_value, sw, time)

def read_xreim(filename, simpy_data):
    """
    This method reads NMR data from a SIMPSON saved with -xreim option.
    """
    with open(filename) as f:
        time = []
        real = []
        imag = []
        for line in f:
            time.append(float(line.split()[0]))
            real.append(float(line.split()[1]))
            imag.append(float(line.split()[2]))

        simpy_data.from_xreim(np.array(time), np.array(real), np.array(imag))

def read_csdf(filename, simpy_data):
    """
    This method reads NMR data from a SIMPSON CSDF file.
    """
    import csdmpy as cp
    data = cp.load(filename)
    hz = data.dimensions[0].coordinates.value
    real = data.dependent_variables[0].components[0].real
    imag = data.dependent_variables[0].components[0].imag
    np_value = np.array(len(hz))
    sw = np.abs(hz[-1] - hz[0])

    simpy_data.from_csdf(real, imag, hz, np_value, sw)

# Function to write Simpson simulations
def write_simp(spinsys, 
               out_name,
               out_format=".inp",
               spin_rate=10e3,
               np=1024,
               proton_freq=400e6,
               start_op="Inx",
               detect_op="Inp",
               crystal_file="rep100",
               gamma_angles=4,
               sw=20e3,
               verbose=0,
               lb=20,
               zerofill=4096,
               method="direct", 
               **kwargs
            ):
    """
    Create a SIMPSON input file with the specified parameters.
    
    Args:
        spinsys: Spin system 
        out_name: Output file name
        out_format: Output format
        spin_rate: Spin rate in Hz
        np: Number of points
        proton_freq: Proton frequency in Hz
        start_op: Start operator
        detect_op: Detect operator
        crystal_file: Crystal file
        gamma_angles: Gamma angles
        sw: Spectral width in Hz
        verbose: Verbose output (0 or 1)
        lb: Line broadening
        zerofill: Zero filling
        method: Simulation method ("direct", "reduced" etc.)
        **kwargs: Additional parameters for Simpcalc
        
    Returns:
        Simpcalc object that can be saved into a Simpson input file
    """
    
    # Create the Simpcalc object with all parameters
    sim = Simpcalc(
        spinsys=spinsys,
        out_name=out_name,
        out_format=out_format,
        spin_rate=spin_rate,
        np=np,
        proton_freq=proton_freq,
        start_op=start_op,
        detect_op=detect_op,
        crystal_file=crystal_file,
        gamma_angles=gamma_angles,
        sw=sw,
        verbose=verbose,
        lb=lb,
        zerofill=zerofill,
        method=method,
        **kwargs
    )
    
    return sim 