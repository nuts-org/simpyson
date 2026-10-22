# SimPYson

SimPYson is a Python interface to [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson), a software for numerical simulations of solid-state NMR experiments. SimPYson reads SIMPSON output, turns DFT calculations into simulated spectra, and writes and runs SIMPSON inputs within Python.

Pick the path that matches what you have in front of you:

<div class="grid cards" markdown>

-   **I have SIMPSON output files**

    ---

    You have `.spe`, `.fid`, `.xreim` or `.csdf` files and want to plot or process them in Python. SIMPSON is not needed.

    [:octicons-arrow-right-24: Analyse SIMPSON output](user_guide/01_reading_files.ipynb)

-   **I have a DFT calculation**

    ---

    You have a CASTEP `.magres` file or a VASP `OUTCAR` and want an simple NMR spectrum. `simulate_spectrum()` needs no further setup.

    [:octicons-arrow-right-24: From DFT to NMR spectrum](user_guide/02_dft_to_simpson.ipynb)

-   **I want to design a simulation**

    ---

    You want full control over the spin system, the pulse sequence and the acquisition: built-in templates or your own Tcl.

    [:octicons-arrow-right-24: Custom simulations with SimpCalc](user_guide/03_simulation_calculator.ipynb)

</div>

If you are new here, see the [Quickstart](quickstart.md) which shows each path in about ten lines. Installation instructions, including when you need the SIMPSON installation, are on the [Installation](installation/install.md) page.

## What you we provide

- **Read SIMPSON output**: `.spe`, `.fid`, `.xreim` and `.csdf` files load into one [`Simpy`][simpyson.simpy.Simpy] object that converts between FID, spectrum (Hz) and spectrum (ppm).
- **DFT to spectrum**: spin systems from CASTEP, Quantum ESPRESSO and VASP calculations, built with [ASE](https://ase-lib.org/) and [Soprano](https://ccp-nc.github.io/soprano/intro.html).
- **Simulations from Python**: [`simulate_spectrum()`][simpyson.calculator.simulate_spectrum] for one-call spectra, [`SimpCalc`][simpyson.calculator.SimpCalc] for full control, with ready-made `no_pulse`, `pulse_90` and `cp_mas` templates and support for your own pulse sequences.
- **A GUI** for quick inspection: `simpyson gui`.

## Why SimPYson?

SimPYson is built by [NUTS - **N**mr **U**nofficial **T**ool**S**et)](https://github.com/nuts-org/) for people who prefer to work within python. Everything is a script, so it is reproducible, easy to batch over many structures, and sits naturally in a DFT-to-NMR pipeline next to other Python tools.

If you would rather work through a graphical interface, these are good alternatives:

- [Simplot](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson): from the developers of SIMPSON, for analysing the results.
- [Simview](https://github.com/zdetos/Simpson-View): from [Zdenek Tosner](https://optimal-nmr.net/about.html) (SIMPSON developer), a GUI to run SIMPSON simulations and plot results.
- [EasyNMR](https://easynmr.pastis.dk/): a cloud-based platform to run and analyse SIMPSON simulations remotely.

We are open to suggestions: please open an issue on [GitHub](https://github.com/nuts-org/simpyson/issues).
