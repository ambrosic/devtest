"""Loss and boundary models for thermodynamic cycle simulation.

The submodels defined here are intentionally composable so the main cycle
integrator can mix and match them depending on configuration or desired
fidelity. Each model expresses its outputs in per–crank-angle units to
align with the integrator's independent variable (radians), leaving the
integration loop to decide how to combine them with stroke-level
derivatives.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol

from .state import CylinderState


# Heat transfer -----------------------------------------------------------------


@dataclass(frozen=True)
class HeatTransferConditions:
    """Geometric and thermal inputs for heat-transfer correlations.

    Parameters
    ----------
    surface_area:
        Instantaneous gas-side surface area [m^2].
    bore:
        Cylinder bore [m]; used by the Woschni correlation for scaling.
    gas_velocity:
        Characteristic gas velocity [m/s] capturing swirl/tumble and
        piston-induced motion.
    wall_temperature:
        Effective wall temperature [K] (piston, liner, and head lumped
        into a single value).
    omega:
        Crankshaft angular speed [rad/s]; used to convert power [W] to
        per-crank-angle energy [J/rad].
    """

    surface_area: float
    bore: float
    gas_velocity: float
    wall_temperature: float
    omega: float


class HeatTransferModel(Protocol):
    """Protocol for heat-transfer correlations."""

    name: str

    def heat_loss_rate(self, state: CylinderState, conditions: HeatTransferConditions) -> float:
        """Return wall heat-loss rate [J/rad] for the given state."""


@dataclass
class WoschniHeatTransfer:
    """Classic Woschni correlation for convective heat transfer.

    The formulation follows the widely used Woschni-style scaling:

    .. math::
        h = c_1 p^{0.8} T^{-0.53} w^{0.8} B^{-0.2}

    where ``p`` is pressure [Pa], ``T`` temperature [K], ``w`` a
    characteristic gas velocity [m/s], and ``B`` the bore [m]. The
    resulting heat-transfer coefficient ``h`` has units W/(m²·K). Heat
    loss per crank angle is derived by multiplying ``h`` by the exposed
    surface area and temperature difference, then dividing by the
    crankshaft speed to convert W to J/rad.
    """

    coefficient: float = 3.26
    exponent_pressure: float = 0.8
    exponent_temperature: float = -0.53
    exponent_velocity: float = 0.8
    bore_exponent: float = -0.2
    name: str = "woschni"

    def heat_loss_rate(self, state: CylinderState, conditions: HeatTransferConditions) -> float:
        if conditions.omega <= 0.0:
            raise ValueError("Crankshaft speed must be positive for heat-transfer calculations")
        if conditions.gas_velocity <= 0.0:
            raise ValueError("Gas velocity must be positive for Woschni correlation")
        if conditions.surface_area <= 0.0:
            raise ValueError("Surface area must be positive for heat-transfer calculations")
        if conditions.bore <= 0.0:
            raise ValueError("Bore must be positive for heat-transfer calculations")

        h_coeff = self.coefficient
        h_coeff *= state.pressure ** self.exponent_pressure
        h_coeff *= state.temperature ** self.exponent_temperature
        h_coeff *= conditions.gas_velocity ** self.exponent_velocity
        h_coeff *= conditions.bore ** self.bore_exponent

        heat_flux = h_coeff * (state.temperature - conditions.wall_temperature)
        heat_power = heat_flux * conditions.surface_area
        return heat_power / conditions.omega


# Friction -----------------------------------------------------------------------


class FrictionModel(Protocol):
    """Protocol for computing friction mean effective pressure (FMEP)."""

    name: str

    def fmep(self, mean_piston_speed: float, peak_pressure: float | None = None) -> float:
        """Return friction mean effective pressure [Pa]."""


@dataclass
class ChenFlynnFriction:
    """Simple Chen–Flynn style friction model.

    The Chen–Flynn formulation expresses FMEP as a polynomial of mean
    piston speed with an optional peak-pressure term:

    .. math::
        FMEP = a + b C_m + c C_m^2 + d p_{max}

    where ``C_m`` is mean piston speed [m/s] and ``p_max`` is peak
    in-cylinder pressure [Pa]. Coefficients are configurable so the model
    can be tuned to measured data.
    """

    a: float = 25_000.0
    b: float = 2_500.0
    c: float = 150.0
    d: float = 0.0
    name: str = "chen-flynn"

    def fmep(self, mean_piston_speed: float, peak_pressure: float | None = None) -> float:
        fmep_value = self.a + self.b * mean_piston_speed + self.c * mean_piston_speed**2
        if peak_pressure is not None:
            fmep_value += self.d * peak_pressure
        return fmep_value


# Blow-by ------------------------------------------------------------------------


@dataclass(frozen=True)
class BlowByConditions:
    """Inputs describing blow-by leakage to the crankcase.

    Parameters
    ----------
    crankcase_pressure:
        Crankcase pressure [Pa]; determines the leakage driving force.
    discharge_coefficient:
        Dimensionless discharge coefficient for the leakage path.
    leakage_area:
        Effective leakage area [m^2].
    omega:
        Crankshaft speed [rad/s] for converting kg/s to kg/rad.
    """

    crankcase_pressure: float
    discharge_coefficient: float
    leakage_area: float
    omega: float


class BlowByModel(Protocol):
    """Protocol for mass leakage models."""

    name: str

    def mass_leakage(self, state: CylinderState, conditions: BlowByConditions) -> float:
        """Return mass leakage rate [kg/rad] (negative for loss)."""


@dataclass
class OrificeBlowByLeakage:
    """Orifice-style blow-by model driven by pressure differential."""

    name: str = "orifice-blow-by"

    def mass_leakage(self, state: CylinderState, conditions: BlowByConditions) -> float:
        if conditions.omega <= 0.0:
            raise ValueError("Crankshaft speed must be positive for blow-by calculations")
        if conditions.leakage_area <= 0.0:
            return 0.0

        delta_p = state.pressure - conditions.crankcase_pressure
        if delta_p <= 0.0:
            return 0.0

        rho = state.density
        mdot = conditions.discharge_coefficient * conditions.leakage_area * math.sqrt(2.0 * rho * delta_p)
        return -mdot / conditions.omega


# Pumping losses -----------------------------------------------------------------


@dataclass(frozen=True)
class PumpingConditions:
    """Boundary pressures and kinematics for pumping work calculation."""

    intake_pressure: float
    exhaust_pressure: float
    dvolume_dtheta: float


class PumpingLossModel(Protocol):
    """Protocol for pumping work estimates."""

    name: str

    def pumping_work(self, conditions: PumpingConditions) -> float:
        """Return pumping work per crank angle [J/rad]."""


@dataclass
class MeanPressurePumpingLoss:
    """Pumping loss using the average of intake and exhaust pressures."""

    name: str = "mean-pressure-pumping"

    def pumping_work(self, conditions: PumpingConditions) -> float:
        boundary_pressure = 0.5 * (conditions.intake_pressure + conditions.exhaust_pressure)
        return boundary_pressure * conditions.dvolume_dtheta


# Aggregation --------------------------------------------------------------------


@dataclass(frozen=True)
class LossTerms:
    """Bundle of per-angle loss contributions."""

    heat_loss: float = 0.0
    mass_leakage: float = 0.0
    pumping_work: float = 0.0
    friction_mep: float = 0.0


@dataclass(frozen=True)
class LossEvaluationContext:
    """Optional inputs required to evaluate the configured loss models."""

    heat_transfer: HeatTransferConditions | None = None
    blow_by: BlowByConditions | None = None
    pumping: PumpingConditions | None = None
    mean_piston_speed: float = 0.0
    peak_pressure: float | None = None


@dataclass
class CycleLossModels:
    """Composable container of loss submodels for the cycle integrator."""

    heat_transfer: HeatTransferModel
    friction: FrictionModel
    blow_by: BlowByModel | None = None
    pumping: PumpingLossModel | None = None

    def evaluate(self, state: CylinderState, context: LossEvaluationContext) -> LossTerms:
        heat_loss = 0.0
        mass_leakage = 0.0
        pumping_work = 0.0

        if context.heat_transfer is not None:
            heat_loss = self.heat_transfer.heat_loss_rate(state, context.heat_transfer)

        if self.blow_by is not None and context.blow_by is not None:
            mass_leakage = self.blow_by.mass_leakage(state, context.blow_by)

        if self.pumping is not None and context.pumping is not None:
            pumping_work = self.pumping.pumping_work(context.pumping)

        friction_mep = self.friction.fmep(context.mean_piston_speed, context.peak_pressure)

        return LossTerms(
            heat_loss=heat_loss,
            mass_leakage=mass_leakage,
            pumping_work=pumping_work,
            friction_mep=friction_mep,
        )

