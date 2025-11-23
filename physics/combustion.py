"""Combustion utilities (Wiebe function)."""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class WiebeHeatRelease:
    """Analytic Wiebe heat-release curve.

    The classic single-parameter Wiebe function is used here to produce
    smooth heat-release profiles. Parameters follow the common notation:

    - `a` and `m` shape the curve (efficiency factor and form exponent).
    - `duration` controls the crank-angle span of the burn [rad].
    - `start_angle` is the crank angle at which the burn commences [rad].
    - `total_heat` is the total energy release per cycle [J].
    """

    a: float
    m: float
    duration: float
    start_angle: float
    total_heat: float

    def mass_fraction_burned(self, theta: float) -> float:
        if theta <= self.start_angle:
            return 0.0
        if theta >= self.start_angle + self.duration:
            return 1.0
        x = (theta - self.start_angle) / self.duration
        return 1.0 - math.exp(-self.a * x ** (self.m + 1.0))

    def dxdtheta(self, theta: float) -> float:
        """Derivative of mass fraction burned with respect to crank angle."""
        if theta <= self.start_angle or theta >= self.start_angle + self.duration:
            return 0.0
        x = (theta - self.start_angle) / self.duration
        coeff = self.a * (self.m + 1.0) * x**self.m / self.duration
        return coeff * math.exp(-self.a * x ** (self.m + 1.0))

    def heat_release_rate(self, theta: float) -> float:
        """Return dQ/dθ [J/rad]."""
        return self.total_heat * self.dxdtheta(theta)
