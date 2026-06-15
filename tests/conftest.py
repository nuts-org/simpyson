"""Shared pytest fixtures for the simpyson test suite."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def write_fake_spe(path: Path, npoints: int = 32, sw: float = 20000.0) -> None:
    lines = ['SIMP', f'NP={npoints}', f'SW={sw}', 'TYPE=SPE', 'DATA']
    lines += [f'{float(i)} 0.0' for i in range(npoints)]
    lines.append('END')
    path.write_text('\n'.join(lines))


@pytest.fixture
def fake_simpson(monkeypatch):
    """Mock the SIMPSON executable.

    Intercepts ``subprocess.run`` so no real SIMPSON binary is needed.
    Writes a minimal ``.spe`` file alongside the input and returns a dict
    that the test can inspect: ``captured['input']`` holds the generated
    ``.in`` text.
    """
    captured = {}

    def fake_run(cmd, **kwargs):
        infile = Path(cmd[1])
        captured['input'] = infile.read_text()
        write_fake_spe(infile.with_suffix('.spe'))
        return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

    monkeypatch.setattr('simpyson.calculator.shutil.which', lambda _: '/usr/bin/simpson')
    monkeypatch.setattr('simpyson.calculator.subprocess.run', fake_run)
    return captured
