"""Physics models and helpers for thermodynamic cycle simulation.

This package provides lightweight interfaces for equation-of-state (EoS)
implementations and per-stroke process models. The default
implementation assumes an ideal gas but the interfaces are designed to
support table-based or real-gas models in the future.
"""

from .eos import EquationOfState, IdealGasEOS, GasProperties
from .state import CylinderState, CylinderDerivatives, StrokeContext
from .combustion import WiebeHeatRelease
from .strokes import (
    StrokeModel,
    IntakeStroke,
    CompressionStroke,
    CombustionStroke,
    ExpansionStroke,
    ExhaustStroke,
)

__all__ = [
    "EquationOfState",
    "IdealGasEOS",
    "GasProperties",
    "CylinderState",
    "CylinderDerivatives",
    "StrokeContext",
    "WiebeHeatRelease",
    "StrokeModel",
    "IntakeStroke",
    "CompressionStroke",
    "CombustionStroke",
    "ExpansionStroke",
    "ExhaustStroke",
]
