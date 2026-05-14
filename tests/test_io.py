"""Tests for simpyson.io — reading and writing SIMPSON files."""
from __future__ import annotations

import os

import numpy as np
import pytest

from simpyson.io import read_fid, read_simp, read_spe, read_xreim
from simpyson.simpy import Simpy

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples', 'read')


# ---------------------------------------------------------------------------
# read_simp (format dispatch)
# ---------------------------------------------------------------------------

class TestReadSimp:
    def test_read_spe_file(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.spe')
        data = read_simp(path)
        assert data.spe is not None
        assert data.spe['np'] == 4096
        assert data.spe['sw'] == 10000

    def test_read_fid_file(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.fid')
        data = read_simp(path)
        assert data.fid is not None
        assert data.fid['np'] == 4096
        assert data.fid['sw'] == 10000

    def test_read_spe_with_b0_nucleus(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.spe')
        data = read_simp(path, b0="400MHz", nucleus="1H")
        assert data.b0 == "400MHz"
        assert data.nucleus == "1H"
        assert data.ppm is not None

    def test_unknown_extension_raises(self, tmp_path):
        dummy = tmp_path / "test.xyz"
        dummy.write_text("dummy")
        with pytest.raises(ValueError, match="Cannot determine file format"):
            read_simp(str(dummy))

    def test_missing_file_raises(self):
        with pytest.raises((OSError, FileNotFoundError)):
            read_simp("nonexistent_file.spe")


# ---------------------------------------------------------------------------
# read_spe
# ---------------------------------------------------------------------------

class TestReadSpe:
    def test_basic_read(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.spe')
        data = Simpy()
        read_spe(path, data)
        assert data.spe is not None
        assert len(data.spe['real']) == 4096
        assert data.spe['sw'] == 10000

    def test_hz_axis_centered(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.spe')
        data = Simpy()
        read_spe(path, data)
        hz = data.spe['hz']
        # Hz axis should be roughly centered around zero (no REF in this file)
        assert hz[0] < 0
        assert hz[-1] > 0

    def test_missing_header_raises(self, tmp_path):
        # File with DATA but missing NP/SW
        bad_file = tmp_path / "bad.spe"
        bad_file.write_text("SIMP\nDATA\n1.0 2.0\nEND")
        data = Simpy()
        with pytest.raises(ValueError, match="Missing required header"):
            read_spe(str(bad_file), data)


# ---------------------------------------------------------------------------
# read_fid
# ---------------------------------------------------------------------------

class TestReadFid:
    def test_basic_read(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.fid')
        data = Simpy()
        read_fid(path, data)
        assert data.fid is not None
        assert len(data.fid['real']) == 4096

    def test_time_axis_starts_at_zero(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.fid')
        data = Simpy()
        read_fid(path, data)
        assert data.fid['time'][0] == pytest.approx(0.0)

    def test_time_axis_is_ms(self):
        path = os.path.join(EXAMPLES_DIR, 'ethanol.fid')
        data = Simpy()
        read_fid(path, data)
        fid = data.fid
        expected_max_ms = (fid['np'] - 1) / fid['sw'] * 1e3
        assert fid['time'][-1] == pytest.approx(expected_max_ms, rel=0.01)


# ---------------------------------------------------------------------------
# Write round-trip
# ---------------------------------------------------------------------------

class TestWriteRoundTrip:
    def test_spe_roundtrip(self, tmp_path):
        """Write a .spe and read it back — data should match."""
        original = Simpy()
        npoints = 32
        sw = 5000.0
        real = np.random.default_rng(1).standard_normal(npoints)
        imag = np.random.default_rng(2).standard_normal(npoints)
        original.from_spe(real, imag, npoints, sw)

        outfile = str(tmp_path / "roundtrip.spe")
        original.write(outfile, format='spe')

        loaded = read_simp(outfile)
        np.testing.assert_allclose(loaded.spe['real'], real, atol=1e-6)
        np.testing.assert_allclose(loaded.spe['imag'], imag, atol=1e-6)

    def test_fid_roundtrip(self, tmp_path):
        """Write a .fid and read it back — data should match."""
        original = Simpy()
        npoints = 32
        sw = 5000.0
        real = np.random.default_rng(3).standard_normal(npoints)
        imag = np.random.default_rng(4).standard_normal(npoints)
        original.from_fid(real, imag, npoints, sw)

        outfile = str(tmp_path / "roundtrip.fid")
        original.write(outfile, format='fid')

        loaded = read_simp(outfile)
        np.testing.assert_allclose(loaded.fid['real'], real, atol=1e-6)
        np.testing.assert_allclose(loaded.fid['imag'], imag, atol=1e-6)
