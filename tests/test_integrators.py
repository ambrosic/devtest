import math
from dataclasses import dataclass

import pytest

from integrators import (
    AdaptiveSettings,
    CrankAngleIntegrator,
    EulerStepper,
    EventBoundary,
    IntegrationResult,
    RK4Stepper,
)
from physics.state import CylinderDerivatives, CylinderState, StrokeContext
from physics.strokes import StrokeModel
from physics.gas_exchange import ValveTimingMap


@dataclass
class DummyStroke(StrokeModel):
    name: str = "dummy"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:  # type: ignore[override]
        return CylinderDerivatives(dm_dtheta=0.0, dT_dtheta=context.wall_heat_loss, dP_dtheta=0.0)


@dataclass
class LinearTemperatureStroke(StrokeModel):
    name: str = "linear-temp"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:  # type: ignore[override]
        return CylinderDerivatives(dm_dtheta=0.0, dT_dtheta=state.temperature, dP_dtheta=0.0)


@dataclass
class GradientStroke(StrokeModel):
    pressure_gradient: float
    temperature_gradient: float
    name: str = "gradient"

    def derivatives(self, state: CylinderState, context: StrokeContext) -> CylinderDerivatives:  # type: ignore[override]
        return CylinderDerivatives(
            dm_dtheta=0.0, dT_dtheta=self.temperature_gradient, dP_dtheta=self.pressure_gradient
        )


def make_state(theta: float = 0.0, temperature: float = 300.0, pressure: float = 1e5) -> CylinderState:
    return CylinderState(theta=theta, volume=0.001, pressure=pressure, temperature=temperature, mass=0.01)


def make_context(**kwargs) -> StrokeContext:
    defaults = dict(dvolume_dtheta=0.0, wall_heat_loss=0.0, manifold_pressure=None, manifold_temperature=None)
    defaults.update(kwargs)
    return StrokeContext(**defaults)


def test_respects_valve_event_boundary():
    timings = ValveTimingMap(ivo=0.3, ivc=1.0, evo=2.0, evc=3.0)
    context = make_context(valve_timings=timings)
    stroke = DummyStroke()

    integrator = CrankAngleIntegrator(EulerStepper(), AdaptiveSettings(max_step=0.5, min_step=1e-4))
    result = integrator.integrate(stroke, make_state(), theta_end=0.8, base_context=context)

    first_step = result.states[1].theta - result.states[0].theta
    assert math.isclose(first_step, 0.3, rel_tol=0, abs_tol=1e-9)
    assert [event.name for event in result.events] == ["intake-open"]


def test_adaptive_step_reduces_on_gradients():
    stroke = GradientStroke(pressure_gradient=5e6, temperature_gradient=0.0)
    context = make_context()
    settings = AdaptiveSettings(max_step=0.2, min_step=1e-4, pressure_tolerance=1e3, safety_factor=0.5)
    integrator = CrankAngleIntegrator(EulerStepper(), settings=settings)

    result = integrator.integrate(stroke, make_state(), theta_end=0.2, base_context=context)
    step_size = result.states[1].theta - result.states[0].theta
    assert step_size < settings.max_step
    assert step_size >= settings.min_step


def test_rk4_tracks_exponential_temperature():
    stroke = LinearTemperatureStroke()
    context = make_context()
    settings = AdaptiveSettings(max_step=0.1, min_step=0.05, pressure_tolerance=1e9, temperature_tolerance=1e9)
    integrator = CrankAngleIntegrator(RK4Stepper(), settings=settings)

    result = integrator.integrate(stroke, make_state(temperature=1.0), theta_end=0.5, base_context=context)
    final_temp = result.states[-1].temperature
    expected = math.exp(0.5)
    assert pytest.approx(expected, rel=1e-3) == final_temp
