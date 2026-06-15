"""
Example: Read and inspect a SIMPSON spectrum file.

This script shows how to load a .spe file, access its data, and optionally
convert the frequency axis to ppm.
"""

from simpyson import read_simp

# --- Read a spectrum file ---
data = read_simp("read/ethanol.spe")

# Access the frequency-domain data
spe = data.spe
print(f"Number of points: {spe['np']}")
print(f"Spectral width:   {spe['sw']} Hz")
print(f"Hz range:          {spe['hz'][0]:.1f} to {spe['hz'][-1]:.1f} Hz")

# --- Convert to ppm (requires b0 and nucleus) ---
data.b0 = "400MHz"
data.nucleus = "1H"

ppm = data.ppm
if ppm is not None:
    print(f"ppm range:         {ppm['ppm'][0]:.2f} to {ppm['ppm'][-1]:.2f} ppm")


# --- Read an FID file ---
fid_data = read_simp("read/ethanol.fid")

fid = fid_data.fid
print(f"\nFID points:  {fid['np']}")
print(f"Time range:  {fid['time'][0]:.3f} to {fid['time'][-1]:.3f} ms")

# The FID can be auto-converted to a spectrum:
spe_from_fid = fid_data.spe
print(f"FFT result:  {len(spe_from_fid['real'])} points")
