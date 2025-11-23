"""Unit tests for loss submodels."""
from __future__ import annotations

import math

import pytest

from physics import (
    BlowByConditions,
    ChenFlynnFriction,
    CycleLossModels,
    HeatTransferConditions,
    LossEvaluationContext,
    MeanPressurePumpingLoss,
    OrificeBlowByLeakage,
    PumpingConditions,
    WoschniHeatTransfer,
    CylinderState,
)


def test_woschni_heat_transfer_matches_manual_calculation():
    state = CylinderState(theta=0.0, volume=5.0e-4, pressure=3.0e6, temperature=900.0, mass=0.0002)
    conditions = HeatTransferConditions(
        surface_area=0.03,
        bore=0.086,
        gas_velocity=25.0,
        wall_temperature=450.0,
        omega=200.0,
    )
    model = WoschniHeatTransfer()

    rate = model.heat_loss_rate(state, conditions)

    h_coeff = model.coefficient
    h_coeff *= state.pressure ** model.exponent_pressure
    h_coeff *= state.temperature ** model.exponent_temperature
    h_coeff *= conditions.gas_velocity ** model.exponent_velocity
    h_coeff *= conditions.bore ** model.bore_exponent
    expected = h_coeff * (state.temperature - conditions.wall_temperature) * conditions.surface_area / conditions.omega

    assert rate == pytest.approx(expected, rel=1e-6)


def test_chen_flynn_fmep_scaling():
    friction = ChenFlynnFriction(a=20_000.0, b=3_000.0, c=200.0, d=1e-6)
    fmep = friction.fmep(mean_piston_speed=5.0, peak_pressure=6.0e6)
    expected = 20_000.0 + 3_000.0 * 5.0 + 200.0 * 25.0 + 1e-6 * 6.0e6
    assert fmep == pytest.approx(expected)


def test_orifice_blow_by_zero_when_not_positive_delta_p():
    state = CylinderState(theta=0.0, volume=4.0e-4, pressure=1.0e5, temperature=300.0, mass=0.0001)
    conditions = BlowByConditions(crankcase_pressure=1.1e5, discharge_coefficient=0.7, leakage_area=2e-6, omega=200.0)
    model = OrificeBlowByLeakage()
    assert model.mass_leakage(state, conditions) == 0.0


def test_orifice_blow_by_mass_flow_direction():
    state = CylinderState(theta=0.0, volume=4.0e-4, pressure=2.0e5, temperature=320.0, mass=0.0001)
    conditions = BlowByConditions(crankcase_pressure=1.0e5, discharge_coefficient=0.6, leakage_area=1.5e-6, omega=100.0)
    model = OrificeBlowByLeakage()
    mass_rate = model.mass_leakage(state, conditions)
    assert mass_rate < 0.0

    rho = state.density
    delta_p = state.pressure - conditions.crankcase_pressure
    expected = -conditions.discharge_coefficient * conditions.leakage_area * math.sqrt(2.0 * rho * delta_p) / conditions.omega
    assert mass_rate == pytest.approx(expected)


def test_mean_pressure_pumping_loss_averages_boundaries():
    model = MeanPressurePumpingLoss()
    conditions = PumpingConditions(intake_pressure=1.0e5, exhaust_pressure=1.3e5, dvolume_dtheta=-1.5e-6)
    work = model.pumping_work(conditions)
    expected = 0.5 * (1.0e5 + 1.3e5) * (-1.5e-6)
    assert work == pytest.approx(expected)


def test_cycle_loss_models_bundle_outputs():
    state = CylinderState(theta=0.0, volume=5.0e-4, pressure=3.0e6, temperature=900.0, mass=0.0002)
    losses = CycleLossModels(
        heat_transfer=WoschniHeatTransfer(),
        friction=ChenFlynnFriction(),
        blow_by=OrificeBlowByLeakage(),
        pumping=MeanPressurePumpingLoss(),
    )

    heat_ctx = HeatTransferConditions(
        surface_area=0.03,
        bore=0.086,
        gas_velocity=20.0,
        wall_temperature=400.0,
        omega=200.0,
    )
    blow_by_ctx = BlowByConditions(
        crankcase_pressure=1.05e5,
        discharge_coefficient=0.65,
        leakage_area=1.5e-6,
        omega=200.0,
    )
    pumping_ctx = PumpingConditions(intake_pressure=1.0e5, exhaust_pressure=1.3e5, dvolume_dtheta=-1.0e-6)

    context = LossEvaluationContext(
        heat_transfer=heat_ctx,
        blow_by=blow_by_ctx,
        pumping=pumping_ctx,
        mean_piston_speed=5.0,
        peak_pressure=3.2e6,
    )

    terms = losses.evaluate(state, context)

    assert terms.heat_loss > 0.0  # positive magnitude for heat rejected to the walls
    assert terms.mass_leakage <= 0.0
    assert terms.pumping_work != 0.0
    assert terms.friction_mep > 0.0
