"""Shared state containers for physics models."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CylinderState:
    """Minimal state description passed between stroke models.

    Parameters
    ----------
    theta : float
        Crank angle [rad]. Stored here for convenience; solvers may also
        track theta separately but embedding it makes debugging easier.
    volume : float
        Instantaneous cylinder volume [m^3]. Geometry/kinematics drive
        this value and its derivative with respect to crank angle.
    pressure : float
        In-cylinder pressure [Pa]. Should be thermodynamically consistent
        with the temperature, mass, and selected equation of state.
    temperature : float
        Bulk gas temperature [K].
    mass : float
        In-cylinder mass [kg]. May change during intake/exhaust strokes.
    """

    theta: float
    volume: float
    pressure: float
    temperature: float
    mass: float

    @property
    def density(self) -> float:
        if self.volume <= 0:
            raise ValueError("Cylinder volume must be positive to compute density")
        return self.mass / self.volume


@dataclass
class CylinderDerivatives:
    """Derivatives with respect to crank angle [rad]."""

    dm_dtheta: float
    dT_dtheta: float
    dP_dtheta: float


@dataclass
class StrokeContext:
    """Auxiliary inputs shared across stroke models."""

    dvolume_dtheta: float
    wall_heat_loss: float = 0.0  # heat loss rate [J/rad]
    manifold_pressure: float | None = None  # optional boundary pressure [Pa]
    manifold_temperature: float | None = None  # optional boundary temperature [K]
