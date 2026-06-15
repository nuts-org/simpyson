"""End-to-end tests against a real SIMPSON installation.

Skipped automatically when the ``simpson`` executable is not on PATH, so
they are CI-safe but exercise the full generate -> run -> read pipeline on
developer machines. The expected values were established empirically (see
the PR #24 review): SIMPSON places +delta at +delta*|nu_L| regardless of
the sign of gamma, and the carrier offset follows the sign of gamma.
"""
from __future__ import annotations

import shutil

import numpy as np
import pytest
from simpyson import simulate_spectrum
from simpyson.utils import get_larmor_freq

requires_simpson = pytest.mark.skipif(
    shutil.which('simpson') is None,
    reason='SIMPSON executable not found in PATH',
)


def _peak_ppm(result) -> float:
    ppm = result.ppm
    assert ppm is not None
    return float(ppm['ppm'][np.argmax(ppm['real'])])


@requires_simpson
def test_t1_smoke_13c_single_peak():
    """A single 13C 50p site must read back at 50 ppm."""
    spinsys = """
channels 13C
nuclei 13C
shift 1 50p 0 0 0 0 0
"""
    result = simulate_spectrum(spinsys, proton_frequency=400e6, spin_rate=10000)
    assert _peak_ppm(result) == pytest.approx(50.0, abs=0.5)


@requires_simpson
@pytest.mark.parametrize('nucleus', ['13C', '29Si'])
def test_t2_axis_convention_both_gamma_signs(nucleus):
    """+50p/-100p must read back at +50/-100 ppm for both gamma signs."""
    spinsys = f"""
channels {nucleus}
nuclei {nucleus} {nucleus}
shift 1 50p 0 0 0 0 0
shift 2 -100p 0 0 0 0 0
"""
    result = simulate_spectrum(spinsys, proton_frequency=400e6, spin_rate=10000)
    ppm_axis = result.ppm['ppm']
    real = result.ppm['real']

    # Find the two largest peaks
    order = np.argsort(real)[::-1]
    found = []
    for idx in order:
        p = ppm_axis[idx]
        if all(abs(p - f) > 5.0 for f in found):
            found.append(float(p))
        if len(found) == 2:
            break

    assert sorted(found) == pytest.approx([-100.0, 50.0], abs=0.5)


@requires_simpson
def test_t3_17o_ct_width_no_wraparound():
    """17O CT MAS: sw must cover the second-order pattern; the peak must not
    sit at the spectrum edge (which would indicate wrap-around)."""
    cq = 3.0e6
    spinsys = f"""
channels 17O
nuclei 17O
shift 1 100p 0 0 0 0 0
quadrupole 1 2 {cq} 0.5 0 0 0
"""
    result = simulate_spectrum(spinsys, proton_frequency=800e6, spin_rate=30000)
    spe = result.spe
    nu_l_hz = abs(get_larmor_freq('800.0MHz', '17O')) * 1e6
    assert spe['sw'] >= cq**2 / nu_l_hz

    # Peak must be inside the central 80% of the window
    peak_idx = int(np.argmax(spe['real']))
    n = int(spe['np'])
    assert 0.1 * n < peak_idx < 0.9 * n


@requires_simpson
def test_t6_33s_second_order_shift_field_scaling():
    """The CT displacement below delta_iso must scale as 1/nu_L^2: the
    displacement ratio between 400 and 800 MHz is (800/400)^2 = 4."""
    delta_iso = 335.7
    spinsys = f"""
channels 33S
nuclei 33S
shift 1 {delta_iso}p 0 0 0 0 0
quadrupole 1 2 0.959e6 1.0 0 0 0
"""
    displacements = {}
    for freq in (400e6, 800e6):
        result = simulate_spectrum(spinsys, proton_frequency=freq, spin_rate=20000)
        displacements[freq] = _peak_ppm(result) - delta_iso

    # Second-order shift is negative (peak below delta_iso) at both fields
    assert displacements[400e6] < 0
    assert displacements[800e6] < 0
    ratio = displacements[400e6] / displacements[800e6]
    assert ratio == pytest.approx(4.0, rel=0.25)
