"""Regression tests for SimpCalc.run()'s output-reading path and file round-trips.

These tests mock the SIMPSON executable so the seam between ``calculator``
and ``io`` (where a ``read_simp()`` keyword mismatch previously crashed every
successful run) is exercised in CI without SIMPSON installed.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest
from simpyson.calculator import SimpCalc, _proton_freq_to_b0
from simpyson.io import read_simp
from simpyson.simpy import Simpy

SPINSYS_1H = """
channels 1H
nuclei 1H
shift 1 5p 0 0 0 0 0
"""


def _make_calc(**overrides):
    params = dict(
        spinsys=SPINSYS_1H,
        pulse_sequence='no_pulse',
        proton_frequency=400e6,
        spin_rate=10000,
        start_operator='Inx',
        detect_operator='Inp',
        np=8,
        sw=20000,
        method='direct',
        crystal_file='rep100',
        gamma_angles=10,
        verbose=0,
    )
    params.update(overrides)
    return SimpCalc(**params)


def _write_fake_spe(path: Path, npoints: int = 8, sw: float = 20000.0) -> None:
    lines = ['SIMP', f'NP={npoints}', f'SW={sw}', 'TYPE=SPE', 'DATA']
    lines += [f'{float(i)} {0.0}' for i in range(npoints)]
    lines.append('END')
    path.write_text('\n'.join(lines))


def test_run_reads_output(tmp_path, monkeypatch):
    """run(read_output=True) must return a Simpy object (regression for
    the read_simp(format=...) keyword crash)."""
    infile = tmp_path / 'sim.in'

    def fake_run(cmd, **kwargs):
        # Simulate SIMPSON writing its output next to the input file
        _write_fake_spe(Path(cmd[1]).with_suffix('.spe'))
        return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

    monkeypatch.setattr('simpyson.calculator.shutil.which', lambda _: '/usr/bin/simpson')
    monkeypatch.setattr('simpyson.calculator.subprocess.run', fake_run)

    calc = _make_calc()
    result = calc.run(filepath=str(infile), read_output=True, delete_files=True)

    assert isinstance(result, Simpy)
    assert result.spe is not None
    assert int(result.spe['np']) == 8
    # b0/nucleus auto-derived, so ppm conversion must be available
    assert result.ppm is not None
    # Cleanup must have removed input and output files
    assert not infile.exists()
    assert not infile.with_suffix('.spe').exists()


def test_run_cleanup_keeps_preexisting_files(tmp_path, monkeypatch):
    """delete_files=True must not delete a pre-existing file in cwd that
    happens to share the output file's name."""
    infile = tmp_path / 'sim.in'
    bystander = Path.cwd() / 'sim.spe'
    bystander.write_text('precious user data')

    def fake_run(cmd, **kwargs):
        _write_fake_spe(Path(cmd[1]).with_suffix('.spe'))
        return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

    monkeypatch.setattr('simpyson.calculator.shutil.which', lambda _: '/usr/bin/simpson')
    monkeypatch.setattr('simpyson.calculator.subprocess.run', fake_run)

    try:
        calc = _make_calc()
        calc.run(filepath=str(infile), read_output=True, delete_files=True)
        assert bystander.exists()
        assert bystander.read_text() == 'precious user data'
    finally:
        bystander.unlink(missing_ok=True)


def test_xreim_write_read_roundtrip(tmp_path):
    """xreim files must round-trip through write() and read_simp()."""
    time = np.array([0.0, 1.0, 2.0])
    real = np.array([1.0, 2.0, 3.0])
    imag = np.array([0.1, 0.2, 0.3])

    outfile = tmp_path / 'test.xreim'
    Simpy().from_xreim(time, real, imag).write(str(outfile), format='xreim')

    loaded = read_simp(str(outfile))
    assert loaded.xreim is not None
    np.testing.assert_allclose(loaded.xreim['time'], time)
    np.testing.assert_allclose(loaded.xreim['real'], real)
    np.testing.assert_allclose(loaded.xreim['imag'], imag)


def test_spe_write_preserves_shifted_axis(tmp_path):
    """A non-symmetric Hz axis (e.g. after add_spectra interpolation) must
    survive a write/read round-trip via the REF header."""
    npoints = 16
    sw = 8000.0
    hz = sw * (np.arange(npoints) / npoints - 0.5) + 1234.5  # shifted axis
    real = np.random.default_rng(0).normal(size=npoints)
    imag = np.zeros(npoints)

    outfile = tmp_path / 'shifted.spe'
    Simpy().from_spe(real, imag, npoints, sw, hz).write(str(outfile), format='spe')

    header = outfile.read_text().splitlines()
    assert any(line.startswith('NP=16') for line in header)  # integer NP

    loaded = read_simp(str(outfile))
    np.testing.assert_allclose(loaded.spe['hz'], hz)
    np.testing.assert_allclose(loaded.spe['real'], real)


def test_csdf_write_read_roundtrip(tmp_path):
    """csdf files must round-trip through write() and read_simp()."""
    pytest.importorskip('csdmpy')

    npoints = 16
    sw = 8000.0
    hz = sw * (np.arange(npoints) / npoints - 0.5)
    real = np.random.default_rng(1).normal(size=npoints)
    imag = np.random.default_rng(2).normal(size=npoints)

    outfile = tmp_path / 'test.csdf'
    Simpy().from_spe(real, imag, npoints, sw, hz).write(str(outfile), format='csdf')

    loaded = read_simp(str(outfile))
    assert loaded.spe is not None
    np.testing.assert_allclose(loaded.spe['hz'], hz)
    np.testing.assert_allclose(loaded.spe['real'], real)
    np.testing.assert_allclose(loaded.spe['imag'], imag)
    np.testing.assert_allclose(loaded.spe['sw'], sw)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (400e6, '400.0MHz'),
        ('400MHz', '400.0MHz'),
        ('400 mhz', '400.0MHz'),
        ('400000kHz', '400.0MHz'),
        ('0.4GHz', '400.0MHz'),
        ('400000000Hz', '400.0MHz'),
        ('400', '400.0MHz'),
    ],
)
def test_proton_freq_to_b0_units(value, expected):
    """All Hz-based unit strings must be normalised to MHz."""
    result = _proton_freq_to_b0(value)
    assert result is not None
    assert float(result.replace('MHz', '')) == pytest.approx(float(expected.replace('MHz', '')))
    assert result.endswith('MHz')


def test_proton_freq_to_b0_invalid():
    assert _proton_freq_to_b0('not a frequency') is None
    assert _proton_freq_to_b0(None) is None
