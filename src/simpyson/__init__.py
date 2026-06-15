"""SimPYson: A Pythonic interface for SIMPSON solid-state NMR simulations."""

from __future__ import annotations

from importlib.metadata import version

from simpyson.calculator import SimpCalc, simulate_spectrum
from simpyson.io import read_simp
from simpyson.simpy import Simpy

__version__ = version("simpyson")

__all__ = [
    "SimpCalc",
    "Simpy",
    "read_simp",
    "simulate_spectrum",
]
