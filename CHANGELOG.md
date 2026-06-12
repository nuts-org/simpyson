# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0]

### Added

- `SimpCalc` class for generating SIMPSON input files from Python.
- `simulate_spectrum()` convenience function with smart defaults (auto-calculates sw, offset, and ref from chemical shifts).
- CPMAS pulse sequence template.
- Support for custom pulse sequences via `CustomPulseSequence`.
- `Simpy` unified data container with lazy FID/spectrum conversion and automatic ppm calculation.
- `.csdf` (csdmpy) file format support, both reading and writing (`Simpy.write(format='csdf')`).
- GUI support for opening `.xreim` / `.csdf` files and saving in all supported formats.
- Comprehensive test suite.

### Changed

- Renamed `SimpSim` to `SimpCalc`.
- Bumped minimum Python version to 3.10.
- Replaced `print()` statements with `warnings.warn()` and `logging`.
- Replaced bare `except Exception:` catches with specific exception types.
- Added type hints and NumPy-style docstrings across all modules.
- Made VASP OUTCAR parser (`converter.py`) more robust with clear error messages for missing sections.
- DRY-ed up isotope data loading in `utils.py`.

### Fixed

- Time unit conversion bug (`* 10e3` should be `* 1e3` for seconds to milliseconds).
- `generate_spinsys()` double-wrapping mutation bug on repeated calls.
- `write_simp()` had wrong parameter names and caused circular import.
- `from_spe()` unnecessarily truncated spectral width with `int()`.
- `add_spectra()` accessed private `_spe_data` and crashed on FID-only input.
- `get_larmor_freq()` had copy-pasted docstring from `hz2ppm()`.
- Uninitialized variables in `read_spe()` / `read_fid()` gave confusing errors on malformed files.
- `read_simp()` was called with `format=` while its keyword was `fmt`, crashing `SimpCalc.run(read_output=True)`, `simulate_spectrum()`, and GUI file opening. The keyword is now `format` everywhere.
- `Simpy.write(format='xreim')` wrote a SIMP-header file without the time axis; it now writes the 3-column `time real imag` format that SIMPSON's `-xreim` flag produces and `read_simp()` expects.
- `Simpy.write(format='spe')` now preserves a shifted frequency axis via the `REF` header and writes `NP` as an integer.
- `read_csdf()` underestimated the spectral width by one bin (used the coordinate span instead of N x step).
- `_proton_freq_to_b0()` passed kHz/GHz strings through unconverted, causing a downstream `ValueError`; all Hz-based units are now normalised to MHz.
- GUI `open_files()` crashed with `UnboundLocalError` on uppercase or unknown file extensions; unsupported files are now skipped with a warning.
- `SimpCalc.run()` cleanup no longer deletes pre-existing files that share the output file name.

## [0.1.1]

- Added a GUI for Simpyson.
- Implemented FID to SPE conversion within the GUI.
- Implemented Hz to ppm conversion within the GUI.

## [0.1.0]

### Added

- Tutorial on converting DFT structures to SIMPSON simulations.
- Tutorial on reading and processing SIMPSON simulation results.
- Added isotope data to convert from Hz to ppm.
- Added `SimpSim` class to prepare SIMPSON input files.
- Added `read_vasp` to convert VASP NMR tensors into a format readable by Soprano.
- Templates for 90-degree pulse `pulse_90` and no-pulse `no_pulse` experiments.

## [0.0.1]

### Added

- The initial release!
