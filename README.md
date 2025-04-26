# SimPYson: A Pythonic Interface for SIMPSON

[![DOI](https://zenodo.org/badge/813117518.svg)](https://doi.org/10.5281/zenodo.14041918)

SimPYson is a Python package designed to simplify the use of [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson), a powerful code to simulate solid-state NMR experiments. Born out of the need to streamline the process, SimPYson makes it easier to prepare input files from DFT calculations and analyze results from SIMPSON simulations all within a Python.

## Features 🤌

- **Convert DFT data to SIMPSON input files**: Prepare SIMPSON input files from DFT data (CASTEP, Quantum Espresso, VASP).

- **Read SIMPSON output files**: Load and manipulate NMR data from SIMPSON `.spe`, `.fid`, and `.xreim` files directly in Python for further analysis and visualization.

- **Graphical User Interface**: You can type `simpyson gui` in your terminal to get Simpyson graphical user interface, and manipulate the data without any coding required.

- **Templates for common experiments**: Use ready-made templates for typical Simpson NMR simulations, currently 90-degree pulse and no-pulse. More soon.

## Documentation 📖

To learn more about simpyson, including some tutorials check the [documentation](https://nuts-org.github.io/simpyson/).

# Planned Features 🔜

- **Expand number of pulse sequences**: Additional templates for more complex NMR experiments.

- **Expand number of DFT codes**: Improve the support of other DFT codes, suggestions are welcomed.
