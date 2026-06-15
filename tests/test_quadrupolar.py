from __future__ import annotations

from simpyson.calculator import SimpCalc, simulate_spectrum
from simpyson.utils import get_spin


def test_get_spin():
    """Test get_spin function for various nuclei."""
    assert get_spin('1H') == 0.5
    assert get_spin('13C') == 0.5
    assert get_spin('23Na') == 1.5   # 3/2
    assert get_spin('27Al') == 2.5   # 5/2
    assert get_spin('14N') == 1.0


def test_simulate_spectrum_quadrupolar():
    """Test that simulate_spectrum sets Inc for quadrupolar nuclei."""
    spinsys_23na = """
    channels 23Na
    nuclei 23Na
    shift 1 0 0 0 0 0 0
    """

    original_SimpCalc = SimpCalc
    captured_params = {}

    class MockSimpCalc:
        def __init__(self, spinsys, **kwargs):
            captured_params.update(kwargs)
            self.spinsys = spinsys
            self.parameters = kwargs

        def run(self, **kwargs):
            return "Mock Result"

    import simpyson.calculator
    simpyson.calculator.SimpCalc = MockSimpCalc

    try:
        # Quadrupolar nucleus should get Inc
        simulate_spectrum(spinsys_23na)
        assert captured_params.get('detect_operator') == 'Inc'

        # Explicit override should be respected
        captured_params.clear()
        simulate_spectrum(spinsys_23na, detect_operator='Inp')
        assert captured_params.get('detect_operator') == 'Inp'

        # Spin-1/2 nucleus should keep Inp
        captured_params.clear()
        spinsys_1h = """
        channels 1H
        nuclei 1H
        shift 1 0 0 0 0 0 0
        """
        simulate_spectrum(spinsys_1h)
        assert captured_params.get('detect_operator') == 'Inp'

    finally:
        simpyson.calculator.SimpCalc = original_SimpCalc
