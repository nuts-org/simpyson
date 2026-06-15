"""Tests for simpyson.converter — hz2ppm, ppm2hz, Larmor frequency."""
from __future__ import annotations

import numpy as np
import pytest

from simpyson.converter import hz2ppm, ppm2hz
from simpyson.utils import get_larmor_freq


def test_larmor_freq_signed_for_negative_gamma():
    """get_larmor_freq stays signed — it is a physical quantity."""
    assert get_larmor_freq('9.4T', '1H') > 0
    assert get_larmor_freq('9.4T', '29Si') < 0


@pytest.mark.parametrize('nucleus', ['1H', '13C', '29Si', '15N'])
def test_hz2ppm_positive_hz_gives_positive_ppm(nucleus):
    """Positive Hz must map to positive ppm regardless of the sign of gamma.

    SIMPSON places a +delta shift at +delta*|nu_L| Hz on its output axis
    (verified empirically with 29Si at 400 MHz — see PR #24 review notes).
    """
    assert hz2ppm(1000.0, '400MHz', nucleus) > 0
    assert ppm2hz(50.0, '400MHz', nucleus) > 0


@pytest.mark.parametrize('nucleus', ['13C', '29Si'])
def test_hz_ppm_roundtrip(nucleus):
    hz = np.array([-5000.0, 0.0, 5000.0])
    back = ppm2hz(hz2ppm(hz, '400MHz', nucleus), '400MHz', nucleus)
    np.testing.assert_allclose(back, hz)
