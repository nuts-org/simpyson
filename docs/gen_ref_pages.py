"""
Generate the code reference pages and navigation.

You never have to run this manually!
"""

from __future__ import annotations

from pathlib import Path

import mkdocs_gen_files

# Modules that are internal (GUI and command-line entry points) get no page.
HIDDEN = {"gui", "cli"}

PUBLIC_API = """\
# Public API

Everything you need for normal use is importable from the top-level package:

```python
from simpyson import SimpCalc, Simpy, read_simp, simulate_spectrum
```

| Name | What it is | Guide |
|---|---|---|
| [`read_simp`][simpyson.io.read_simp] | Read a SIMPSON `.spe`, `.fid`, `.xreim` or `.csdf` file into a `Simpy` object. | [Analyse SIMPSON output](../user_guide/01_reading_files.ipynb) |
| [`Simpy`][simpyson.simpy.Simpy] | Container for FID, spectrum (Hz) and spectrum (ppm), with conversion and export. | [Analyse SIMPSON output](../user_guide/01_reading_files.ipynb) |
| [`simulate_spectrum`][simpyson.calculator.simulate_spectrum] | Simulate a spectrum from a spin system with automatic settings. Needs SIMPSON. | [From DFT to NMR spectrum](../user_guide/02_dft_to_simpson.ipynb) |
| [`SimpCalc`][simpyson.calculator.SimpCalc] | Build, save and run a SIMPSON input file. Needs SIMPSON to run. | [Custom simulations with SimpCalc](../user_guide/03_simulation_calculator.ipynb) |

## Pulse sequences

The pulse sequence templates are in `simpyson.templates` and are normally selected by name in `SimpCalc`:

| Name | Class |
|---|---|
| `'no_pulse'` | [`NoPulse`][simpyson.templates.NoPulse] |
| `'pulse_90'` | [`Pulse90`][simpyson.templates.Pulse90] |
| `'cp_mas'` | [`CPMAS`][simpyson.templates.CPMAS] |

Your own Tcl is wrapped by [`CustomPulseSequence`][simpyson.templates.CustomPulseSequence]; subclass [`PulseSequenceTemplate`][simpyson.templates.PulseSequenceTemplate] to write a reusable template.

## Helpers

Other modules contain helpers used in the guides: [`read_vasp`][simpyson.converter.read_vasp] (VASP `OUTCAR` to ASE `Atoms`), [`hz2ppm`][simpyson.converter.hz2ppm] and [`ppm2hz`][simpyson.converter.ppm2hz], and [`add_spectra`][simpyson.utils.add_spectra].

## Graphical interface

Simpyson provides a simple graphical user interface (GUI) for easy interaction with the software. Which is being tested, it can be
accessed via the command line wiht `simpyson gui`. Currently it is not documented.
"""

nav = mkdocs_gen_files.Nav()

nav[("Public API",)] = "public_api.md"
with mkdocs_gen_files.open("reference/public_api.md", "w") as fd:
    fd.write(PUBLIC_API)

for path in sorted(Path("src").rglob("*.py")):
    module_path = path.relative_to("src").with_suffix("")
    doc_path = path.relative_to("src").with_suffix(".md")
    full_doc_path = Path("reference", doc_path)

    parts = tuple(module_path.parts)

    if parts[-1] in ("__init__", "__main__") or parts[-1] in HIDDEN:
        continue

    nav[parts] = doc_path.as_posix()

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        ident = ".".join(parts)
        fd.write(f"::: {ident}")

    mkdocs_gen_files.set_edit_path(full_doc_path, Path("../") / path)

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
