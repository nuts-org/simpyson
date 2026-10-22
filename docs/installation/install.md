# Installation

## Install SimPYson

SimPYson needs **Python 3.10 or newer**. Install it with pip:

```bash
pip install simpyson
```

All Python dependencies are installed automatically.

## When do you need SIMPSON?

The [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) program is only needed to *run* simulations. Everything else works without it.

| I want to... | SIMPSON needed? |
|---|---|
| Read, plot, convert or combine `.spe`, `.fid`, `.xreim`, `.csdf` files | No |
| Read DFT output (`.magres`, VASP `OUTCAR`) and build a spin system | No |
| Generate a SIMPSON input file with `SimpCalc` (`.save()`, `.print()`) | No |
| Run a simulation with `simulate_spectrum()` or `SimpCalc.run()` | **Yes** |

To run simulations, install SIMPSON (the [SIMPSON website](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) has downloads and build instructions) so that the `simpson` command is on your `PATH`. Running `simpson` without arguments prints its version, which is a quick way to check:

```bash
simpson
```

If the executable lives somewhere else, pass its location to `SimpCalc.run(simpson_path="/path/to/simpson")`.

## Verify the installation

Open a Python console and run:

```python
import simpyson

print(simpyson.__version__)
```

To check that the GUI works:

```bash
simpyson gui
```
