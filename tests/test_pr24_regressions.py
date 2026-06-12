"""CI-safe regression tests for the PR #24 review findings.

Covers the negative-gamma axis/offset conventions (verified empirically
against SIMPSON in a separate session), spectral-width estimation, parsing
helpers, and template generation. No SIMPSON installation required: runs
are exercised through a mocked subprocess.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pytest
from simpyson.calculator import (
    SimpCalc,
    _extract_nucleus,
    _proton_freq_to_b0,
    simulate_spectrum,
)
from simpyson.converter import hz2ppm, ppm2hz
from simpyson.simpy import Simpy
from simpyson.templates import Pulse90
from simpyson.utils import add_spectra, get_larmor_freq

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_fake_spe(path: Path, npoints: int = 32, sw: float = 20000.0) -> None:
    lines = ['SIMP', f'NP={npoints}', f'SW={sw}', 'TYPE=SPE', 'DATA']
    lines += [f'{float(i)} 0.0' for i in range(npoints)]
    lines.append('END')
    path.write_text('\n'.join(lines))


@pytest.fixture
def fake_simpson(monkeypatch):
    """Mock the SIMPSON executable and capture the generated input file."""
    captured = {}

    def fake_run(cmd, **kwargs):
        infile = Path(cmd[1])
        captured['input'] = infile.read_text()
        _write_fake_spe(infile.with_suffix('.spe'))
        return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

    monkeypatch.setattr('simpyson.calculator.shutil.which', lambda _: '/usr/bin/simpson')
    monkeypatch.setattr('simpyson.calculator.subprocess.run', fake_run)
    return captured


def _par_value(input_text: str, name: str) -> float:
    m = re.search(rf'^\s*{name}\s+(\S+)', input_text, re.MULTILINE)
    assert m, f"par entry {name!r} not found in generated input"
    return float(m.group(1))


def _par_variable(input_text: str, name: str) -> float:
    m = re.search(rf'^\s*variable\s+{name}\s+(\S+)', input_text, re.MULTILINE)
    assert m, f"par variable {name!r} not found in generated input"
    return float(m.group(1))


# ---------------------------------------------------------------------------
# Axis sign convention (T2): SIMPSON puts +delta at +delta*|nu_L|
# ---------------------------------------------------------------------------


def test_larmor_freq_signed_for_negative_gamma():
    """get_larmor_freq stays signed -- it is a physical quantity."""
    assert get_larmor_freq('9.4T', '1H') > 0
    assert get_larmor_freq('9.4T', '29Si') < 0


@pytest.mark.parametrize('nucleus', ['1H', '13C', '29Si', '15N'])
def test_hz2ppm_uses_magnitude(nucleus):
    """Positive Hz must map to positive ppm for every gamma sign."""
    assert hz2ppm(1000.0, '400MHz', nucleus) > 0
    assert ppm2hz(50.0, '400MHz', nucleus) > 0


@pytest.mark.parametrize('nucleus', ['13C', '29Si'])
def test_hz_ppm_roundtrip(nucleus):
    hz = np.array([-5000.0, 0.0, 5000.0])
    back = ppm2hz(hz2ppm(hz, '400MHz', nucleus), '400MHz', nucleus)
    np.testing.assert_allclose(back, hz)


# ---------------------------------------------------------------------------
# simulate_spectrum offset/ref signs (T2 fix 2)
# ---------------------------------------------------------------------------

SPINSYS_13C = """
channels 13C
nuclei 13C 13C
shift 1 10p 0 0 0 0 0
shift 2 50p 0 0 0 0 0
"""

SPINSYS_29SI = """
channels 29Si
nuclei 29Si 29Si
shift 1 50p 0 0 0 0 0
shift 2 -100p 0 0 0 0 0
"""


def test_offset_ref_signs_positive_gamma(fake_simpson):
    simulate_spectrum(SPINSYS_13C, proton_frequency=400e6, spin_rate=10000)
    text = fake_simpson['input']
    center_hz = ppm2hz(30.0, '400.0MHz', '13C')  # (10+50)/2
    assert _par_variable(text, 'offset') == pytest.approx(center_hz)
    assert _par_variable(text, 'ref') == pytest.approx(-center_hz)


def test_offset_ref_signs_negative_gamma(fake_simpson):
    """Carrier offset follows sign(gamma); ref is always -center_hz."""
    simulate_spectrum(SPINSYS_29SI, proton_frequency=400e6, spin_rate=10000)
    text = fake_simpson['input']
    center_hz = ppm2hz(-25.0, '400.0MHz', '29Si')  # (50-100)/2, negative
    assert center_hz < 0
    assert _par_variable(text, 'offset') == pytest.approx(-center_hz)  # sign flipped
    assert _par_variable(text, 'ref') == pytest.approx(-center_hz)


# ---------------------------------------------------------------------------
# Spectral width estimation with |nu_L| (finding 2)
# ---------------------------------------------------------------------------


def test_sw_positive_for_negative_gamma_shifts(fake_simpson):
    simulate_spectrum(SPINSYS_29SI, proton_frequency=400e6, spin_rate=10000)
    text = fake_simpson['input']
    sw = _par_value(text, 'sw')
    width_hz = abs(ppm2hz(50.0, '400.0MHz', '29Si') - ppm2hz(-100.0, '400.0MHz', '29Si'))
    assert sw > 0
    assert sw >= 2 * width_hz
    assert sw % 10000 == 0  # rounded up to the spinning rate


def test_sw_floor_static_15n(fake_simpson):
    """Single peak, static: sw must hit the 10-ppm floor at |nu_L|."""
    spinsys = """
channels 15N
nuclei 15N
shift 1 20p 0 0 0 0 0
"""
    simulate_spectrum(spinsys, proton_frequency=800e6, spin_rate=0)
    sw = _par_value(fake_simpson['input'], 'sw')
    nu_l_hz = abs(get_larmor_freq('800.0MHz', '15N')) * 1e6
    assert sw == pytest.approx(nu_l_hz * 10e-6)
    assert sw > 0


def test_sw_quadrupolar_ct_mas_17o(fake_simpson):
    """17O CT MAS: width dominated by Cq^2/|nu_L| second-order broadening."""
    cq = 3.0e6
    spinsys = f"""
channels 17O
nuclei 17O
shift 1 100p 0 0 0 0 0
quadrupole 1 2 {cq} 0.5 0 0 0
"""
    simulate_spectrum(spinsys, proton_frequency=800e6, spin_rate=30000)
    text = fake_simpson['input']
    # quadrupolar nucleus -> central-transition detection by default
    assert 'Inc' in text
    sw = _par_value(text, 'sw')
    nu_l_hz = abs(get_larmor_freq('800.0MHz', '17O')) * 1e6
    assert sw > 0
    assert sw >= cq**2 / nu_l_hz
    assert sw % 30000 == 0


def test_sw_quadrupolar_static_17o(fake_simpson):
    """Quad-only static spinsys: sw == Cq^2/|nu_L| (CT) exactly."""
    cq = 3.0e6
    spinsys = f"""
channels 17O
nuclei 17O
quadrupole 1 2 {cq} 0.5 0 0 0
"""
    simulate_spectrum(spinsys, proton_frequency=800e6, spin_rate=0)
    sw = _par_value(fake_simpson['input'], 'sw')
    nu_l_hz = abs(get_larmor_freq('800.0MHz', '17O')) * 1e6
    assert sw == pytest.approx(max(cq**2 / nu_l_hz, nu_l_hz * 10e-6))
    assert sw > 0


# ---------------------------------------------------------------------------
# _proton_freq_to_b0 (finding 3)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('value', 'expected_mhz'),
    [
        (8e8, 800.0),
        ('8e8', 800.0),
        ('8E8', 800.0),
        ('4.0e8Hz', 400.0),
        ('800MHz', 800.0),
        ('400', 400.0),
    ],
)
def test_proton_freq_to_b0_scientific_notation(value, expected_mhz):
    result = _proton_freq_to_b0(value)
    assert result is not None
    assert result.endswith('MHz')
    assert float(result[:-3]) == pytest.approx(expected_mhz)


# ---------------------------------------------------------------------------
# Newline-eating regexes (finding 5)
# ---------------------------------------------------------------------------


def test_extract_nucleus_multiline():
    assert _extract_nucleus(SPINSYS_29SI) == '29Si'
    assert _extract_nucleus("channels 1H 13C\nnuclei 1H 13C\n") == '1H'


def test_nucleus_from_detect_operator_does_not_leak_lines(fake_simpson):
    """With the old [\\w\\s]+ regex, an out-of-range site index picked tokens
    from the *following* spinsys line (e.g. 'shift') as the nucleus."""
    spinsys = """
channels 1H
nuclei 1H
shift 1 5p 0 0 0 0 0
"""
    result = simulate_spectrum(
        spinsys, proton_frequency=400e6, spin_rate=10000,
        detect_operator='I2p',  # site 2 does not exist
    )
    assert result.nucleus == '1H'  # falls back to the channels line


# ---------------------------------------------------------------------------
# Pulse90 multi-channel padding (T4)
# ---------------------------------------------------------------------------


def test_pulse90_single_channel_unpadded():
    code = Pulse90().generate_code()
    assert 'pulse $par(pH) $par(plH) $par(phH)\n' in code


def test_pulse90_two_channels_padded():
    code = Pulse90(num_channels=2).generate_code()
    assert 'pulse $par(pH) $par(plH) $par(phH) 0 0\n' in code


def test_pulse90_padding_via_simpcalc():
    calc = SimpCalc(
        spinsys="channels 1H 13C\nnuclei 1H 13C\nshift 1 5p 0 0 0 0 0\nshift 2 50p 0 0 0 0 0",
        pulse_sequence='pulse_90',
        proton_frequency=400e6, spin_rate=10000, sw=20000, np=1024,
        start_operator='I1x', detect_operator='I2p', method='direct',
        crystal_file='rep100', gamma_angles=10, verbose=0,
    )
    assert 'pulse $par(pH) $par(plH) $par(phH) 0 0' in str(calc)


# ---------------------------------------------------------------------------
# add_spectra spectral width invariant (finding 4)
# ---------------------------------------------------------------------------


def test_add_spectra_common_sw_is_full_width():
    npoints = 64
    sw = 8000.0
    hz_a = sw * (np.arange(npoints) / npoints - 0.5)
    hz_b = hz_a + 2000.0  # shifted axis forces interpolation path
    rng = np.random.default_rng(0)
    a = Simpy().from_spe(rng.normal(size=npoints), np.zeros(npoints), npoints, sw, hz_a)
    b = Simpy().from_spe(rng.normal(size=npoints), np.zeros(npoints), npoints, sw, hz_b)

    combined = add_spectra([a, b])
    spe = combined.spe
    step = spe['hz'][1] - spe['hz'][0]
    assert spe['sw'] == pytest.approx(step * spe['np'])
