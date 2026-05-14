"""Tests for simpyson.simpy — the Simpy data container."""
from __future__ import annotations

import numpy as np
import pytest

from simpyson.simpy import Simpy


# ---------------------------------------------------------------------------
# Construction and basic properties
# ---------------------------------------------------------------------------

def test_empty_simpy():
    s = Simpy()
    assert s.b0 is None
    assert s.nucleus is None
    assert s.fid is None
    assert s.spe is None
    assert s.ppm is None
    assert s.xreim is None


def test_b0_and_nucleus():
    s = Simpy(b0="400MHz", nucleus="13C")
    assert s.b0 == "400MHz"
    assert s.nucleus == "13C"


# ---------------------------------------------------------------------------
# from_spe / from_fid / from_xreim round-trips
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_spectrum():
    """A Simpy with synthetic spectrum data."""
    s = Simpy()
    npoints = 64
    sw = 10000.0
    hz = np.linspace(-sw / 2, sw / 2, npoints)
    real = np.sin(hz / 100)
    imag = np.cos(hz / 100)
    s.from_spe(real, imag, npoints, sw, hz)
    return s


@pytest.fixture
def simple_fid():
    """A Simpy with synthetic FID data."""
    s = Simpy()
    npoints = 64
    sw = 10000.0
    real = np.random.default_rng(42).standard_normal(npoints)
    imag = np.random.default_rng(43).standard_normal(npoints)
    s.from_fid(real, imag, npoints, sw)
    return s


def test_from_spe_stores_data(simple_spectrum):
    spe = simple_spectrum.spe
    assert spe is not None
    assert len(spe['real']) == 64
    assert spe['sw'] == 10000.0
    assert 'hz' in spe


def test_from_fid_stores_data(simple_fid):
    fid = simple_fid.fid
    assert fid is not None
    assert len(fid['real']) == 64
    assert fid['sw'] == 10000.0
    assert 'time' in fid


def test_from_fid_time_axis_is_ms(simple_fid):
    """Time axis should be in milliseconds (seconds * 1e3)."""
    fid = simple_fid.fid
    npoints = fid['np']
    sw = fid['sw']
    dt = 1.0 / sw
    expected_max_ms = (npoints - 1) * dt * 1e3
    assert fid['time'][-1] == pytest.approx(expected_max_ms, rel=0.01)


def test_from_fid_default_time_matches_explicit():
    """Passing time=None should give the same result as computing manually."""
    npoints = 128
    sw = 20000.0
    real = np.ones(npoints)
    imag = np.zeros(npoints)

    s1 = Simpy()
    s1.from_fid(real, imag, npoints, sw)  # time computed internally

    dt = 1.0 / sw
    time_manual = np.arange(npoints) * dt * 1e3
    s2 = Simpy()
    s2.from_fid(real, imag, npoints, sw, time=time_manual)

    np.testing.assert_allclose(s1.fid['time'], s2.fid['time'])


def test_from_xreim():
    s = Simpy()
    time = np.linspace(0, 1, 32)
    real = np.sin(time)
    imag = np.cos(time)
    s.from_xreim(time, real, imag)
    assert s.xreim is not None
    assert len(s.xreim['time']) == 32


# ---------------------------------------------------------------------------
# Lazy conversions (FID <-> spectrum)
# ---------------------------------------------------------------------------

def test_fid_to_spe_conversion(simple_fid):
    """Accessing .spe on FID-only data should trigger FFT."""
    spe = simple_fid.spe
    assert spe is not None
    assert 'hz' in spe
    assert len(spe['real']) == 64


def test_spe_to_fid_conversion(simple_spectrum):
    """Accessing .fid on spectrum-only data should trigger iFFT."""
    fid = simple_spectrum.fid
    assert fid is not None
    assert 'time' in fid
    assert len(fid['real']) == 64


def test_fid_spe_roundtrip(simple_fid):
    """FID -> SPE -> FID should approximately recover the original."""
    original_real = simple_fid.fid['real'].copy()
    spe = simple_fid.spe  # triggers FFT
    # Access FID through a fresh conversion
    s2 = Simpy()
    s2.from_spe(spe['real'], spe['imag'], spe['np'], spe['sw'], spe['hz'])
    recovered = s2.fid
    np.testing.assert_allclose(recovered['real'], original_real, atol=1e-10)


# ---------------------------------------------------------------------------
# PPM conversion
# ---------------------------------------------------------------------------

def test_ppm_requires_b0_and_nucleus(simple_spectrum):
    """ppm should be None when b0/nucleus are not set."""
    assert simple_spectrum.ppm is None


def test_ppm_computed_when_b0_nucleus_set():
    s = Simpy(b0="400MHz", nucleus="1H")
    npoints = 32
    sw = 5000.0
    hz = np.linspace(-sw / 2, sw / 2, npoints)
    s.from_spe(np.ones(npoints), np.zeros(npoints), npoints, sw, hz)
    ppm = s.ppm
    assert ppm is not None
    assert 'ppm' in ppm
    assert len(ppm['ppm']) == npoints


def test_ppm_invalidated_on_b0_change():
    s = Simpy(b0="400MHz", nucleus="1H")
    s.from_spe(np.ones(8), np.zeros(8), 8, 1000.0)
    _ = s.ppm  # trigger ppm calculation
    assert 'ppm' in s.spe
    s.b0 = "800MHz"
    assert 'ppm' not in s.spe  # should be invalidated


# ---------------------------------------------------------------------------
# Copy
# ---------------------------------------------------------------------------

def test_copy_is_independent(simple_spectrum):
    clone = simple_spectrum.copy()
    assert clone.spe is not None
    # Mutating clone should not affect original
    clone.spe['real'][0] = 999999
    assert simple_spectrum.spe['real'][0] != 999999


# ---------------------------------------------------------------------------
# Write (SIMPSON format)
# ---------------------------------------------------------------------------

def test_write_spe(simple_spectrum, tmp_path):
    outfile = str(tmp_path / "test.spe")
    simple_spectrum.write(outfile, format='spe')
    with open(outfile) as f:
        content = f.read()
    assert content.startswith('SIMP')
    assert 'NP=64' in content
    assert 'SW=10000' in content
    assert 'TYPE=SPE' in content
    assert 'DATA' in content
    assert 'END' in content


def test_write_csv(simple_spectrum, tmp_path):
    outfile = str(tmp_path / "test.csv")
    simple_spectrum.write(outfile, format='csv')
    data = np.loadtxt(outfile, delimiter=',', skiprows=1)
    assert data.shape == (64, 2)


def test_write_unsupported_format(simple_spectrum, tmp_path):
    with pytest.raises(ValueError, match="Unsupported save format"):
        simple_spectrum.write(str(tmp_path / "test.xyz"), format='xyz')


def test_write_no_data(tmp_path):
    s = Simpy()
    with pytest.raises(ValueError, match="No data"):
        s.write(str(tmp_path / "test.csv"), format='csv')


# ---------------------------------------------------------------------------
# Method chaining
# ---------------------------------------------------------------------------

def test_method_chaining():
    s = Simpy()
    result = s.from_spe([1, 2], [0, 0], 2, 100.0)
    assert result is s
    result2 = s.write.__wrapped__(s, "dummy") if hasattr(s.write, '__wrapped__') else None
    # Just verify from_fid also chains
    s2 = Simpy()
    assert s2.from_fid([1, 2], [0, 0], 2, 100.0) is s2
