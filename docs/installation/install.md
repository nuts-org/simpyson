# Installation

## Prerequisites

- **Python >= 3.10**
- **SIMPSON** -- required only for *running* simulations (not needed for reading/converting files). Download from the [SIMPSON website](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) and make sure the `simpson` command is available in your PATH.

## Install SimPYson

If you have Python and pip installed, you can install SimPYson with:

```bash
pip install git+https://github.com/nuts-org/simpyson.git
```

All dependencies are installed automatically.

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
