"""Tests for simpyson.utils — isotope lookups, add_spectra, simple_spinsys."""
from __future__ import annotations

import numpy as np
import pytest

from simpyson.simpy import Simpy
from simpyson.utils import (
    _load_isotope_data,
    add_spectra,
    get_gamma,
    get_larmor_freq,
    get_spin,
)


# ---------------------------------------------------------------------------
# Isotope data helpers
# ---------------------------------------------------------------------------

class TestIsotopeData:
    def test_load_isotope_data_hydrogen(self):
        row = _load_isotope_data("1H")
        assert "Gamma" in row
        assert "Spin" in row

    def test_load_isotope_data_carbon(self):
        row = _load_isotope_data("13C")
        assert row["Gamma"] != 0

    def test_load_isotope_data_invalid(self):
        with pytest.raises(ValueError, match="not found"):
            _load_isotope_data("999Xx")

    def test_get_gamma_hydrogen(self):
        gamma = get_gamma("1H")
        assert isinstance(gamma, float)
        assert gamma > 0

    def test_get_gamma_carbon(self):
        gamma_c = get_gamma("13C")
        gamma_h = get_gamma("1H")
        # Carbon gamma should be smaller than hydrogen
        assert gamma_c < gamma_h

    def test_get_spin_hydrogen(self):
        assert get_spin("1H") == 0.5

    def test_get_spin_sodium(self):
        assert get_spin("23Na") == 1.5

    def test_get_spin_invalid(self):
        with pytest.raises(ValueError, match="not found"):
            get_spin("999Xx")


# ---------------------------------------------------------------------------
# Larmor frequency
# ---------------------------------------------------------------------------

class TestLarmorFreq:
    def test_larmor_freq_tesla(self):
        freq = get_larmor_freq("9.4T", "1H")
        # At 9.4T, proton Larmor should be ~400 MHz
        assert 390 < freq < 410

    def test_larmor_freq_mhz(self):
        freq = get_larmor_freq("400MHz", "1H")
        assert freq == pytest.approx(400.0, rel=0.01)

    def test_larmor_freq_carbon_at_400mhz(self):
        freq_h = get_larmor_freq("400MHz", "1H")
        freq_c = get_larmor_freq("400MHz", "13C")
        # Carbon freq should be ~1/4 of proton
        assert 0.20 < freq_c / freq_h < 0.30

    def test_larmor_freq_invalid_unit(self):
        with pytest.raises(ValueError, match="T or MHz"):
            get_larmor_freq("400kHz", "1H")

    def test_larmor_freq_invalid_nucleus(self):
        with pytest.raises(ValueError, match="not found"):
            get_larmor_freq("400MHz", "999Xx")


# ---------------------------------------------------------------------------
# add_spectra
# ---------------------------------------------------------------------------

def _make_spe(real_values, sw=1000.0, hz=None):
    """Helper to create a Simpy with spectrum data."""
    s = Simpy()
    n = len(real_values)
    s.from_spe(real_values, np.zeros(n), n, sw, hz=hz)
    return s


class TestAddSpectra:
    def test_empty_list(self):
        assert add_spectra([]) is None

    def test_single_spectrum(self):
        s = _make_spe([1, 2, 3])
        result = add_spectra([s])
        assert result is not None
        np.testing.assert_array_equal(result.spe['real'], [1, 2, 3])

    def test_two_spectra_sum(self):
        s1 = _make_spe([1, 2, 3])
        s2 = _make_spe([4, 5, 6])
        result = add_spectra([s1, s2])
        np.testing.assert_array_equal(result.spe['real'], [5, 7, 9])

    def test_add_spectra_does_not_mutate_original(self):
        s1 = _make_spe([1, 2, 3])
        original = s1.spe['real'].copy()
        s2 = _make_spe([10, 20, 30])
        _ = add_spectra([s1, s2])
        np.testing.assert_array_equal(s1.spe['real'], original)

    def test_add_spectra_with_b0_override(self):
        s1 = _make_spe([1, 1])
        s2 = _make_spe([2, 2])
        result = add_spectra([s1, s2], b0="400MHz", nucleus="1H")
        assert result.b0 == "400MHz"
        assert result.nucleus == "1H"

    def test_add_spectra_no_spe_raises(self):
        s = Simpy()  # no data at all
        with pytest.raises(ValueError, match="no frequency-domain data"):
            add_spectra([s])

    def test_add_spectra_different_hz_axes(self):
        """Spectra with different Hz ranges are interpolated onto a common grid."""
        # s1: peak at 100 Hz, covering 0–200 Hz
        hz1 = np.linspace(0, 200, 201)
        real1 = np.zeros(201)
        real1[100] = 1.0  # peak at 100 Hz
        s1 = _make_spe(real1, sw=200.0, hz=hz1)

        # s2: peak at 300 Hz, covering 200–400 Hz (no overlap)
        hz2 = np.linspace(200, 400, 201)
        real2 = np.zeros(201)
        real2[100] = 1.0  # peak at 300 Hz
        s2 = _make_spe(real2, sw=200.0, hz=hz2)

        result = add_spectra([s1, s2])
        # Common grid covers 0–400 Hz
        assert result.spe['hz'][0] == pytest.approx(0.0)
        assert result.spe['hz'][-1] == pytest.approx(400.0)
        # Both peaks should be present; find max positions
        real_combined = result.spe['real']
        hz_combined = result.spe['hz']
        peak_indices = np.where(real_combined > 0.5)[0]
        peak_hz = hz_combined[peak_indices]
        assert any(abs(p - 100.0) < 2.0 for p in peak_hz), "Peak near 100 Hz missing"
        assert any(abs(p - 300.0) < 2.0 for p in peak_hz), "Peak near 300 Hz missing"

    def test_add_spectra_different_hz_does_not_mutate(self):
        """Interpolation path should not mutate originals."""
        hz1 = np.linspace(0, 100, 101)
        real1 = np.ones(101)
        s1 = _make_spe(real1, sw=100.0, hz=hz1)
        original1 = s1.spe['real'].copy()

        hz2 = np.linspace(50, 150, 101)
        real2 = np.ones(101) * 2
        s2 = _make_spe(real2, sw=100.0, hz=hz2)

        _ = add_spectra([s1, s2])
        np.testing.assert_array_equal(s1.spe['real'], original1)


def test_add_spectra_common_sw_is_increment_times_np():
    """sw must equal increment * np, not end-minus-start (one-bin difference)."""
    npoints = 64
    sw = 8000.0
    hz_a = sw * (np.arange(npoints) / npoints - 0.5)
    hz_b = hz_a + 2000.0  # shifted axis forces interpolation path
    rng = np.random.default_rng(0)
    a = _make_spe(rng.normal(size=npoints), sw=sw, hz=hz_a)
    b = _make_spe(rng.normal(size=npoints), sw=sw, hz=hz_b)

    combined = add_spectra([a, b])
    spe = combined.spe
    step = spe['hz'][1] - spe['hz'][0]
    assert spe['sw'] == pytest.approx(step * spe['np'])
