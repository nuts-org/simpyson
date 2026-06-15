"""Tests for simpyson.converter — hz2ppm, ppm2hz, Larmor frequency, read_vasp."""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pytest

from simpyson.converter import hz2ppm, ppm2hz, read_vasp
from simpyson.utils import get_larmor_freq

OUTCAR = Path(__file__).parent.parent / 'examples' / 'write' / 'AlPO-14.OUTCAR'


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

# Test read_vasp and convertion to magres style
@pytest.mark.skipif(
    not OUTCAR.exists(),
    reason=f'OUTCAR not present at {OUTCAR}',
)
def test_read_vasp_alpo14_structure():
    """read_vasp returns an ASE Atoms object with 48 atoms (8Al + 32O + 8P)."""
    atoms = read_vasp(str(OUTCAR), format='vasp-out')
    syms = atoms.get_chemical_symbols()
    assert len(atoms) == 48
    assert syms.count('Al') == 8
    assert syms.count('O') == 32
    assert syms.count('P') == 8


@pytest.mark.skipif(
    not OUTCAR.exists(),
    reason=f'OUTCAR not present at {OUTCAR}',
)
def test_read_vasp_alpo14_tensor_shapes():
    """ms and efg arrays have shape (n_atoms, 3, 3)."""
    atoms = read_vasp(str(OUTCAR), format='vasp-out')
    ms = atoms.get_array('ms')
    efg = atoms.get_array('efg')
    assert ms.shape == (48, 3, 3)
    assert efg.shape == (48, 3, 3)
