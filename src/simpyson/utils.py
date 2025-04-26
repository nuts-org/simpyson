import os
import json
import numpy as np

def get_larmor_freq(b0, nucleus, isotope_file=None):
    """
    Convert Hz values to ppm values.
    
    Args:
        b0 (str): Magnetic field strength (e.g., '400MHz' or '9.4T')
        nucleus (str): Nucleus type (e.g., '1H' or '13C')
        isotope_file (str, optional): Path to isotope data file. If None, uses default.
    
    Returns:
        float: Larmor Frequency
        
    Raises:
        ValueError: If B0 unit is invalid or nucleus not found
    """
    
    if isotope_file is None:
        dir = os.path.dirname(os.path.realpath(__file__))
        isotope_file = os.path.join(dir, 'isotope_data.json')
    
    isotope = int(''.join(filter(str.isdigit, nucleus)))
    element = ''.join(filter(str.isalpha, nucleus)).capitalize()
    b0_unit = ''.join(filter(str.isalpha, b0)).lower()
    
    with open(isotope_file) as f:
        data = json.load(f)
        if element in data and str(isotope) in data[element]:
            gamma = data[element][str(isotope)]['Gamma']
        else:
            raise ValueError(f'Nucleus {nucleus} not found in isotope data.')
    
    if b0_unit == 't':
        b0_value = float(''.join(filter(lambda x: x.isdigit() or x == '.', b0)))
        larmor_freq = gamma * 1e7 * b0_value / (2 * np.pi * 1e6)
        ppm = hz / np.abs(larmor_freq)
    elif b0_unit == 'mhz':
        b0_value = float(''.join(filter(lambda x: x.isdigit() or x == '.', b0))) 
        gamma_h = data['H']['1']['Gamma']
        b0_value_T = 2 * np.pi * b0_value * 1e6 / (gamma_h * 1e7)
        larmor_freq = gamma * 1e7 * b0_value_T / (2 * np.pi * 1e6)
    else:
        raise ValueError('B0 unit must be T or MHz.')
    
    return larmor_freq

def add_spectra(spectra_list, b0=None, nucleus=None):
    """Combine multiple Simpy objects into a single spectrum."""
    if not spectra_list:
        return None
        
    result = spectra_list[0].copy()
    
    if b0:
        result.b0 = b0
    if nucleus:
        result.nucleus = nucleus
        
    for spectrum in spectra_list[1:]:
        result._spe_data['real'] += spectrum.spe['real']
        result._spe_data['imag'] += spectrum.spe['imag']
        
    if result.b0 and result.nucleus:
        result._compute_ppm()
        
    return result