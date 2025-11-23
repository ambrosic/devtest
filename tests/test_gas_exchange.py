import math

import pytest

from physics import CylinderState, GasProperties, IdealGasEOS, IntakeStroke, OrificeFlow, StrokeContext, ValveTimingMap


def test_valve_timing_map_handles_wraparound():
    timings = ValveTimingMap(ivo=3.5 * math.pi, ivc=0.5 * math.pi, evo=math.pi, evc=2.5 * math.pi)
    assert timings.is_intake_open(3.75 * math.pi)
    assert timings.is_intake_open(0.25 * math.pi)
    assert not timings.is_intake_open(math.pi)
    assert timings.is_exhaust_open(1.25 * math.pi)
    assert not timings.is_exhaust_open(3.0 * math.pi)


def test_orifice_flow_mass_flow_per_angle_matches_manual_calc():
    properties = GasProperties(gas_constant=287.0, gamma=1.4)
    eos = IdealGasEOS(properties)
    orifice = OrificeFlow(eos=eos, discharge_coefficient=0.8, effective_area=2.0e-4)

    p_up = 2.0e5
    p_down = 1.0e5
    temp_up = 300.0
    omega = 200.0

    rho = p_up / ((properties.cp - properties.cv) * temp_up)
    expected = 0.8 * 2.0e-4 * math.sqrt(2.0 * rho * (p_up - p_down)) / omega

    assert orifice.mass_flow_per_angle(p_up, p_down, omega, temperature_upstream=temp_up) == pytest.approx(expected)


def test_intake_stroke_respects_valve_timing_and_orifice_flow():
    properties = GasProperties(gas_constant=287.0, gamma=1.35)
    eos = IdealGasEOS(properties)
    orifice = OrificeFlow(eos=eos, discharge_coefficient=0.7, effective_area=1.0e-4)
    timings = ValveTimingMap(ivo=0.0, ivc=math.pi, evo=2.0 * math.pi, evc=3.0 * math.pi)

    stroke = IntakeStroke(eos=eos, orifice_flow=orifice)
    state = CylinderState(theta=0.25 * math.pi, volume=5.0e-4, pressure=9.0e4, temperature=310.0, mass=0.0002)
    context = StrokeContext(
        dvolume_dtheta=0.0,
        manifold_pressure=1.2e5,
        manifold_temperature=320.0,
        wall_heat_loss=0.0,
        valve_timings=timings,
        omega=150.0,
    )

    derivatives = stroke.derivatives(state, context)
    cp, cv = eos.specific_heats()
    gas_constant = cp - cv
    rho_manifold = context.manifold_pressure / (gas_constant * context.manifold_temperature)
    expected_dm = 0.7 * 1.0e-4 * math.sqrt(2.0 * rho_manifold * (context.manifold_pressure - state.pressure)) / context.omega

    assert derivatives.dm_dtheta == pytest.approx(expected_dm)

    closed_state = CylinderState(theta=1.25 * math.pi, volume=5.0e-4, pressure=9.0e4, temperature=310.0, mass=0.0002)
    closed_context = StrokeContext(
        dvolume_dtheta=0.0,
        manifold_pressure=1.2e5,
        manifold_temperature=320.0,
        wall_heat_loss=0.0,
        valve_timings=timings,
        omega=150.0,
    )

    closed_derivatives = stroke.derivatives(closed_state, closed_context)
    assert closed_derivatives.dm_dtheta == 0.0
