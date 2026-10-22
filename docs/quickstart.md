# Quickstart

Each tab is a complete script. The first two download a small example file from the SimPYson repository, so you can run them as written. You can also replace the file name with your own data when you are ready.

=== "Read SIMPSON output"

    Load a spectrum and plot it in ppm. No SIMPSON installation is needed.

    ```python
    from urllib.request import urlretrieve

    import matplotlib.pyplot as plt

    from simpyson import read_simp

    # Example file from the simpyson repository; use your own .spe/.fid instead
    urlretrieve(
        "https://raw.githubusercontent.com/nuts-org/simpyson/main/examples/read/ethanol.spe",
        "ethanol.spe",
    )

    data = read_simp("ethanol.spe", b0="400MHz", nucleus="1H")
    plt.plot(data.ppm["ppm"], data.ppm["real"])
    plt.xlim(8, 0)
    plt.xlabel("$^1$H (ppm)")
    plt.show()
    ```

    `b0` and `nucleus` are only needed for the ppm axis. The Hz axis is `data.spe["hz"]`. See [Analyse SIMPSON output](user_guide/01_reading_files.ipynb).

=== "DFT to spectrum"

    Turn a CASTEP `.magres` file into a simulated <sup>1</sup>H spectrum. **Requires SIMPSON on your `PATH`.**

    ```python
    from urllib.request import urlretrieve

    import matplotlib.pyplot as plt
    from ase.io import read
    from soprano.calculate.nmr.simpson import write_spinsys

    from simpyson import simulate_spectrum

    # Example CASTEP file from the simpyson repository; use your own .magres instead
    urlretrieve(
        "https://raw.githubusercontent.com/nuts-org/simpyson/main/examples/write/ethanol.magres",
        "ethanol.magres",
    )

    atoms = read("ethanol.magres")
    ch3 = atoms[[3, 4, 5]]  # the three CH3 hydrogens
    spinsys = write_spinsys(ch3, use_ms=True, grad={"H": -1}, ref={"H": 31.7})

    spectrum = simulate_spectrum(spinsys)
    plt.plot(spectrum.ppm["ppm"], spectrum.ppm["real"])
    plt.xlim(8, -2)
    plt.xlabel("$^1$H (ppm)")
    plt.show()
    ```

    `grad` and `ref` convert the DFT shielding $\sigma$ to a chemical shift $\delta$; the reference value depends on your functional and code. See [From DFT to NMR spectrum](user_guide/02_dft_to_simpson.ipynb).

=== "Custom simulation"

    Build a SIMPSON input with a built-in pulse sequence and run it. **Requires SIMPSON on your `PATH`.**

    ```python
    from simpyson import SimpCalc

    spinsys = """
    channels 13C
    nuclei 13C 13C
    shift 1 50p 0 0 0 0 0
    shift 2 20p 0 0 0 0 0
    """

    calc = SimpCalc(
        spinsys,
        pulse_sequence="pulse_90",
        proton_frequency=400e6,
        spin_rate=10000,
        start_operator="Inz",
        detect_operator="Inp",
        method="direct",
        crystal_file="rep100",
        gamma_angles=10,
        np=4096,
        sw=20000,
        verbose=0,
    )
    result = calc.run()
    ```

    `result` is a [`Simpy`][simpyson.simpy.Simpy] object, like the one `read_simp` returns. Use `calc.print()` to see the generated SIMPSON input. See [Custom simulations with SimpCalc](user_guide/03_simulation_calculator.ipynb).

!!! note "Where is SIMPSON needed?"
    Reading and converting files never needs SIMPSON. Anything that runs a simulation (`simulate_spectrum()`, `SimpCalc.run()`) does. See [Installation](installation/install.md#when-do-you-need-simpson).
