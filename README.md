# SimPYson: A Pythonic Interface for SIMPSON

[![DOI](https://zenodo.org/badge/813117518.svg)](https://doi.org/10.5281/zenodo.14041918)
![Python - Version](https://img.shields.io/pypi/pyversions/simpyson)
![PyPI - Version](https://img.shields.io/pypi/v/simpyson?color=blue)

SimPYson is a Python package that makes it easier to work with [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson), a code for simulating solid-state NMR experiments. It handles preparing input files from DFT calculations, running simulations, and reading results — all from Python.

## Features

- **Run SIMPSON simulations from Python**: Use `simulate_spectrum()` for smart defaults or `SimpCalc` for full control over spin systems, pulse sequences, and output.
- **Convert DFT data to SIMPSON input files**: Prepare spin systems from CASTEP, Quantum Espresso, and VASP calculations via ASE and Soprano.
- **Read SIMPSON output files**: Load `.spe`, `.fid`, `.xreim`, and `.csdf` files into a unified `Simpy` object with automatic FID↔spectrum conversion and ppm scaling.
- **Pulse sequence templates**: Built-in templates for no-pulse, 90° pulse, and CPMAS experiments. Custom Tcl sequences are also supported.
- **Graphical User Interface**: Launch with `simpyson gui` to inspect and process spectra without any coding.

## Quick Start

**Read a SIMPSON output file:**

```python
from simpyson import read_simp

data = read_simp("ethanol.spe", b0="400MHz", nucleus="1H")
print(data.ppm['ppm'])   # ppm axis, auto-calculated
print(data.spe['hz'])    # Hz axis
```

Accessing `.fid` on a spectrum file (or `.spe` on a FID) triggers automatic conversion via FFT — no manual processing needed.

**Simulate a spectrum directly from Python:**

```python
from simpyson import simulate_spectrum

spinsys = """
channels 13C
nuclei 13C 13C
shift 1 10p 0 0 0 0 0
shift 2 50p 0 0 0 0 0
"""

result = simulate_spectrum(spinsys, proton_frequency=400e6, spin_rate=10000)
print(result.ppm['ppm'])
```

`simulate_spectrum()` automatically estimates the spectral width and carrier offset from the chemical shifts. For full control over all parameters, use `SimpCalc` directly — see the [documentation](https://nuts-org.github.io/simpyson/).

## Installation

```bash
pip install git+https://github.com/nuts-org/simpyson.git
```

Requires Python ≥ 3.10 and a working [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) installation for running simulations.

## Documentation

Full documentation and tutorials are available at [nuts-org.github.io/simpyson](https://nuts-org.github.io/simpyson/).

## Planned Features

- Additional pulse sequence templates for more complex NMR experiments.
- Broader support for DFT codes — suggestions welcome.
