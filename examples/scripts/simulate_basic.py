"""
Example: Simulate a basic NMR spectrum using simulate_spectrum().

This script demonstrates the simplest way to run a SIMPSON simulation from
Python.  It creates a spin system by hand, runs a no-pulse simulation, and
prints basic info about the result.

NOTE: Requires a working SIMPSON installation in your PATH.
"""

from simpyson import simulate_spectrum

# Define a simple spin system with two 13C sites
spinsys = """
channels 13C
nuclei 13C 13C
shift 1 10p 0 0 0 0 0
shift 2 50p 0 0 0 0 0
"""

# simulate_spectrum() auto-calculates sw, offset, and ref from the shifts
result = simulate_spectrum(
    spinsys,
    proton_frequency=400e6,
    spin_rate=10000,
)

# The result is a Simpy object with spectrum data
spe = result.spe
print(f"Number of points: {spe['np']}")
print(f"Spectral width:   {spe['sw']} Hz")

# ppm is auto-calculated because proton_frequency and nucleus are known
ppm = result.ppm
if ppm is not None:
    print(f"ppm range: {ppm['ppm'][0]:.1f} to {ppm['ppm'][-1]:.1f}")
