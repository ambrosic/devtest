"""Per-stroke process models.

These are light-weight, quasi-thermodynamic models intended to drive a
future crank-angle integrator. Each stroke implements a common interface
so they can be composed into full-cycle simulations.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from .combustion import WiebeHeatRelease
from .eos import EquationOfState
from .state import CylinderDerivatives, CylinderState, StrokeContext


class StrokeModel(Protocol):
    name: str

    def derivatives(
        self, state: CylinderState, context: StrokeContext
    ) -> CylinderDerivatives:
        """Return derivatives with respect to crank angle."""


@dataclass
class IntakeStroke:
    """Simple throttling-based intake model.

    Mass flow is driven by the pressure difference between the intake
    manifold and the cylinder. Enthalpy of the incoming charge is
    approximated with the manifold temperature.
    """

    eos: EquationOfState
    flow_coefficient: float = 1e-6  # kg/(Pa*rad)
    name: str = "intake"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:
        if context.manifold_pressure is None or context.manifold_temperature is None:
            raise ValueError("Intake stroke requires manifold pressure and temperature")

        delta_p = context.manifold_pressure - state.pressure
        dm_dtheta = self.flow_coefficient * delta_p

        # Energy balance: dU = h_in dm - p dV - q_wall
        cp, cv = self.eos.specific_heats()
        h_in = cp * context.manifold_temperature
        dU_dtheta = h_in * dm_dtheta - state.pressure * context.dvolume_dtheta - context.wall_heat_loss
        dT_dtheta = dU_dtheta / (state.mass * cv)
        dP_dtheta = self.eos.pressure(state.density, state.temperature + dT_dtheta * 1e-6)
        dP_dtheta = (dP_dtheta - state.pressure) / 1e-6
        return CylinderDerivatives(dm_dtheta=dm_dtheta, dT_dtheta=dT_dtheta, dP_dtheta=dP_dtheta)


@dataclass
class CompressionStroke:
    """Polytropic compression model."""

    eos: EquationOfState
    polytropic_index: float = 1.32
    name: str = "compression"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:
        dV = context.dvolume_dtheta
        dT_dtheta = -(self.polytropic_index - 1.0) * state.temperature / state.volume * dV
        dP_dtheta = -self.polytropic_index * state.pressure / state.volume * dV
        return CylinderDerivatives(dm_dtheta=0.0, dT_dtheta=dT_dtheta, dP_dtheta=dP_dtheta)


@dataclass
class CombustionStroke:
    """Closed-volume combustion with Wiebe heat release."""

    eos: EquationOfState
    wiebe: WiebeHeatRelease
    name: str = "combustion"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:
        cp, cv = self.eos.specific_heats()
        dQ_dtheta = self.wiebe.heat_release_rate(state.theta)
        dT_dtheta = (dQ_dtheta - state.pressure * context.dvolume_dtheta - context.wall_heat_loss) / (state.mass * cv)
        dP_dtheta = self.eos.pressure(state.density, state.temperature + dT_dtheta * 1e-6)
        dP_dtheta = (dP_dtheta - state.pressure) / 1e-6
        return CylinderDerivatives(dm_dtheta=0.0, dT_dtheta=dT_dtheta, dP_dtheta=dP_dtheta)


@dataclass
class ExpansionStroke:
    """Isentropic expansion model."""

    eos: EquationOfState
    name: str = "expansion"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:
        gamma = self.eos.gamma()
        dV = context.dvolume_dtheta
        dT_dtheta = -(gamma - 1.0) * state.temperature / state.volume * dV
        dP_dtheta = -gamma * state.pressure / state.volume * dV
        return CylinderDerivatives(dm_dtheta=0.0, dT_dtheta=dT_dtheta, dP_dtheta=dP_dtheta)


@dataclass
class ExhaustStroke:
    """Simple blowdown/exhaust model."""

    eos: EquationOfState
    flow_coefficient: float = 1e-6
    name: str = "exhaust"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:
        if context.manifold_pressure is None or context.manifold_temperature is None:
            raise ValueError("Exhaust stroke requires manifold pressure and temperature")

        delta_p = state.pressure - context.manifold_pressure
        dm_dtheta = -self.flow_coefficient * delta_p

        cp, cv = self.eos.specific_heats()
        h_out = cp * state.temperature
        dU_dtheta = -h_out * dm_dtheta - state.pressure * context.dvolume_dtheta - context.wall_heat_loss
        dT_dtheta = dU_dtheta / (state.mass * cv)
        dP_dtheta = self.eos.pressure(state.density, state.temperature + dT_dtheta * 1e-6)
        dP_dtheta = (dP_dtheta - state.pressure) / 1e-6
        return CylinderDerivatives(dm_dtheta=dm_dtheta, dT_dtheta=dT_dtheta, dP_dtheta=dP_dtheta)
