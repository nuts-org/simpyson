"""
Example: Build a SIMPSON input file step by step using SimpCalc.

This script shows how to use SimpCalc directly for full control over all
simulation parameters, pulse sequence selection, and output configuration.
The generated .in file can be run with SIMPSON or inspected manually.
"""

from simpyson import SimpCalc

# Define a spin system (two 1H sites with a shift difference)
spinsys = """
channels 1H
nuclei 1H 1H
shift 1 5p 0 0 0 0 0
shift 2 8p 0 0 0 0 0
"""

# Create a calculator with a no-pulse experiment
calc = SimpCalc(
    spinsys=spinsys,
    pulse_sequence="no_pulse",
    proton_frequency=400e6,
    spin_rate=10000,
    start_operator="Inx",
    detect_operator="Inp",
    np=4096,
    sw=10000,
    method="direct",
    crystal_file="rep100",
    gamma_angles=10,
    verbose=0,
    out_name="example_output",
    out_format="spe",
    lb=50,
    zerofill=4096,
)

# Preview the generated SIMPSON input file
print(calc)

# Save to disk
calc.save("example_simulation.in")
print("\nSaved to example_simulation.in")

# To actually run the simulation (requires SIMPSON in PATH):
# result = calc.run(read_output=True, delete_files=True)
# print(result.spe['hz'][:10])
