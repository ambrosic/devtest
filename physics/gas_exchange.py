"""Valve timing and gas-exchange utilities."""
from __future__ import annotations

import math
from dataclasses import dataclass

from .eos import EquationOfState


@dataclass(frozen=True)
class ValveTimingMap:
    """Intake/exhaust valve timing definition for a single cylinder.

    Parameters
    ----------
    ivo : float
        Intake valve opening crank angle [rad].
    ivc : float
        Intake valve closing crank angle [rad].
    evo : float
        Exhaust valve opening crank angle [rad].
    evc : float
        Exhaust valve closing crank angle [rad].
    cycle_period : float
        Crank-angle span for a full cycle [rad]. Defaults to 4π for a
        four-stroke engine.
    """

    ivo: float
    ivc: float
    evo: float
    evc: float
    cycle_period: float = 4.0 * math.pi

    def is_intake_open(self, theta: float) -> bool:
        return self._angle_in_window(theta, self.ivo, self.ivc)

    def is_exhaust_open(self, theta: float) -> bool:
        return self._angle_in_window(theta, self.evo, self.evc)

    def _angle_in_window(self, theta: float, start: float, end: float) -> bool:
        if self.cycle_period <= 0.0:
            raise ValueError("Cycle period must be positive for valve timing evaluation")

        angle = theta % self.cycle_period
        start_angle = start % self.cycle_period
        end_angle = end % self.cycle_period

        if start_angle <= end_angle:
            return start_angle <= angle <= end_angle
        return angle >= start_angle or angle <= end_angle


@dataclass(frozen=True)
class OrificeFlow:
    """Orifice-based mass-flow model using a discharge coefficient."""

    eos: EquationOfState
    discharge_coefficient: float
    effective_area: float
    name: str = "orifice-flow"

    def mass_flow_rate(
        self,
        pressure_upstream: float,
        pressure_downstream: float,
        temperature_upstream: float | None = None,
        density_upstream: float | None = None,
    ) -> float:
        """Return mass flow rate [kg/s] from upstream to downstream."""

        if self.discharge_coefficient <= 0.0 or self.effective_area <= 0.0:
            return 0.0

        delta_p = pressure_upstream - pressure_downstream
        if delta_p <= 0.0:
            return 0.0

        rho = density_upstream
        if rho is None:
            if temperature_upstream is None:
                raise ValueError("Upstream temperature required to compute density")
            cp, cv = self.eos.specific_heats()
            gas_constant = cp - cv
            if gas_constant <= 0.0:
                raise ValueError("Specific gas constant must be positive")
            rho = pressure_upstream / (gas_constant * temperature_upstream)

        return self.discharge_coefficient * self.effective_area * math.sqrt(2.0 * rho * delta_p)

    def mass_flow_per_angle(
        self,
        pressure_upstream: float,
        pressure_downstream: float,
        omega: float,
        temperature_upstream: float | None = None,
        density_upstream: float | None = None,
    ) -> float:
        """Return mass flow per crank angle [kg/rad] from upstream to downstream."""

        if omega <= 0.0:
            raise ValueError("Crankshaft speed must be positive to convert to per-angle flow")
        mdot = self.mass_flow_rate(
            pressure_upstream=pressure_upstream,
            pressure_downstream=pressure_downstream,
            temperature_upstream=temperature_upstream,
            density_upstream=density_upstream,
        )
        return mdot / omega
