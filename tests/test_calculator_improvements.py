from __future__ import annotations

import pytest
from simpyson.calculator import SimpCalc, _extract_nucleus, _proton_freq_to_b0
from simpyson.converter import ppm2hz
from simpyson.utils import get_larmor_freq


def test_validation_spinsys():
    """Test that spinsys cannot be None."""
    with pytest.raises(ValueError):
        SimpCalc(spinsys=None)


def test_validation_pulse_sequence():
    """Test invalid pulse sequence types."""
    with pytest.raises(ValueError):
        SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence=123)


def test_missing_parameters():
    """Test missing required parameters."""
    calc = SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence="pulse_90")
    with pytest.raises(ValueError, match="Missing required parameters"):
        calc.generate_par()


def test_dry_run():
    """Test dry_run option."""
    calc = SimpCalc(
        spinsys="spinsys { channels 1H }",
        pulse_sequence="pulse_90",
        proton_frequency=400e6,
        spin_rate=10000,
        start_operator="I1z",
        detect_operator="I1p",
        np=1024,
        sw=20000,
        method="direct",
        crystal_file="rep100",
        gamma_angles=10,
        verbose=0,
    )
    # Should not raise FileNotFoundError even if simpson is missing
    cmd = calc.run(dry_run=True)
    assert isinstance(cmd, str)
    assert "simpson" in cmd


# ---------------------------------------------------------------------------
# _proton_freq_to_b0 — unit strings and scientific notation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('value', 'expected_mhz'),
    [
        (400e6,           400.0),
        ('400MHz',        400.0),
        ('400 mhz',       400.0),
        ('400000kHz',     400.0),
        ('0.4GHz',        400.0),
        ('400000000Hz',   400.0),
        ('400',           400.0),
        ('8e8',           800.0),
        ('8E8',           800.0),
        ('4.0e8Hz',       400.0),
    ],
)
def test_proton_freq_to_b0(value, expected_mhz):
    result = _proton_freq_to_b0(value)
    assert result is not None
    assert result.endswith('MHz')
    assert float(result[:-3]) == pytest.approx(expected_mhz)


def test_proton_freq_to_b0_invalid():
    assert _proton_freq_to_b0('not a frequency') is None
    assert _proton_freq_to_b0(None) is None


# ---------------------------------------------------------------------------
# _extract_nucleus — multiline spinsys must not eat adjacent lines
# ---------------------------------------------------------------------------


def test_extract_nucleus_single_channel():
    assert _extract_nucleus("channels 29Si\nnuclei 29Si\nshift 1 50p 0 0 0 0 0") == '29Si'


def test_extract_nucleus_multichannel_returns_first():
    assert _extract_nucleus("channels 1H 13C\nnuclei 1H 13C\n") == '1H'


def test_nucleus_from_detect_operator_does_not_eat_spinsys_lines(fake_simpson):
    """With the old [\\w\\s]+ regex an out-of-range site index could pick up
    the *following* spinsys line (e.g. 'shift') as the nucleus name."""
    from simpyson.calculator import simulate_spectrum
    spinsys = "channels 1H\nnuclei 1H\nshift 1 5p 0 0 0 0 0"
    result = simulate_spectrum(spinsys, proton_frequency=400e6, spin_rate=10000,
                               detect_operator='I2p')
    assert result.nucleus == '1H'
