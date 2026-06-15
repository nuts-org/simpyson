`SimPYson` is a Python package developed to streamline the use of [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson), a simulation software for solid-state NMR. Born out of the frustration of trying to learn and use SIMPSON on our own, SimPYson aims to simplify the process of converting DFT calculation results into SIMPSON input files, running simulations, and reading SIMPSON output all within Python.

## Current features

- **Convert DFT data to SIMPSON input**: Prepare SIMPSON input files from DFT calculations (CASTEP, Quantum Espresso, VASP) using [ASE](https://ase-lib.org/) and [Soprano](https://ccp-nc.github.io/soprano/intro.html).
- **Run simulations from Python**: Use `SimpCalc` or `simulate_spectrum()` to generate, execute, and read SIMPSON simulations without writing Tcl manually.
- **Read and process output files**: Load `.spe`, `.fid`, `.xreim`, and `.csdf` files from SIMPSON into the unified `Simpy` data container.
- **Pulse sequence templates**: Ready-made templates for no-pulse, 90-degree pulse, and CPMAS experiments, plus support for custom pulse sequences.
- **Graphical User Interface**: Launch with `simpyson gui` to view, combine, and convert spectra interactively (PyQt5 + Plotly).

## Why choose SimPYson?

One possible advantage lies in its seamless integration with other Python packages, making it easy to incorporate into your existing Python workflow. However, depending on your needs/taste, there are a few alternatives that may better suit you, here are some examples:

- [Simplot](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson): From the developers of SIMPSON, offers a user-friendly interface to analyze the results.
- [Simview](https://github.com/zdetos/Simpson-View): From [Zdenek Tosner](https://optimal-nmr.net/about.html) (SIMPSON dev), provides a GUI to run SIMPSON simulations and plot results within a clean, intuitive interface.
- [EasyNMR](https://easynmr.pastis.dk/): A cloud-based platform that allows you to run and analyze SIMPSON simulations remotely.

## Future
We are open to suggestions.
