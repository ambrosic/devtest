"""Crank-angle integrators with plug-in steppers."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, Protocol, Sequence

from physics.combustion import WiebeHeatRelease
from physics.gas_exchange import ValveTimingMap
from physics.state import CylinderDerivatives, CylinderState, StrokeContext
from physics.strokes import StrokeModel


class Stepper(Protocol):
    """Protocol for plug-in integration steppers."""

    name: str

    def step(
        self,
        rhs: Callable[[CylinderState, StrokeContext], CylinderDerivatives],
        state: CylinderState,
        step: float,
        context_provider: Callable[[float, CylinderState], StrokeContext],
    ) -> CylinderState:
        ...


@dataclass(frozen=True)
class AdaptiveSettings:
    """Controls for adaptive step sizing."""

    max_step: float = 0.05
    min_step: float = 1e-5
    pressure_tolerance: float = 5e3
    temperature_tolerance: float = 2.5
    safety_factor: float = 0.5
    event_tolerance: float = 1e-9


@dataclass(frozen=True)
class EventBoundary:
    """Represents a crank-angle event that must be hit exactly."""

    angle: float
    name: str


@dataclass
class IntegrationResult:
    """Container for integration traces."""

    states: list[CylinderState]
    events: list[EventBoundary] = field(default_factory=list)


class EulerStepper:
    """Simple forward-Euler stepper."""

    name = "euler"

    def step(
        self,
        rhs: Callable[[CylinderState, StrokeContext], CylinderDerivatives],
        state: CylinderState,
        step: float,
        context_provider: Callable[[float, CylinderState], StrokeContext],
    ) -> CylinderState:
        context = context_provider(state.theta, state)
        derivs = rhs(state, context)
        return _advance_state(state, derivs, step, context)


class RK4Stepper:
    """Classical 4th-order Runge–Kutta stepper."""

    name = "rk4"

    def step(
        self,
        rhs: Callable[[CylinderState, StrokeContext], CylinderDerivatives],
        state: CylinderState,
        step: float,
        context_provider: Callable[[float, CylinderState], StrokeContext],
    ) -> CylinderState:
        half_step = step * 0.5

        ctx1 = context_provider(state.theta, state)
        k1 = rhs(state, ctx1)

        s2 = _advance_state(state, k1, half_step, ctx1)
        ctx2 = context_provider(s2.theta, s2)
        k2 = rhs(s2, ctx2)

        s3 = _advance_state(state, k2, half_step, ctx2)
        ctx3 = context_provider(s3.theta, s3)
        k3 = rhs(s3, ctx3)

        s4 = _advance_state(state, k3, step, ctx3)
        ctx4 = context_provider(s4.theta, s4)
        k4 = rhs(s4, ctx4)

        dm = (k1.dm_dtheta + 2 * k2.dm_dtheta + 2 * k3.dm_dtheta + k4.dm_dtheta) / 6.0
        dT = (k1.dT_dtheta + 2 * k2.dT_dtheta + 2 * k3.dT_dtheta + k4.dT_dtheta) / 6.0
        dP = (k1.dP_dtheta + 2 * k2.dP_dtheta + 2 * k3.dP_dtheta + k4.dP_dtheta) / 6.0

        return _advance_state(
            state,
            CylinderDerivatives(dm_dtheta=dm, dT_dtheta=dT, dP_dtheta=dP),
            step,
            ctx1,
        )


class CrankAngleIntegrator:
    """Adaptive crank-angle integrator with pluggable steppers."""

    def __init__(
        self,
        stepper: Stepper,
        settings: AdaptiveSettings | None = None,
        context_provider: Callable[[float, CylinderState], StrokeContext] | None = None,
    ) -> None:
        self.stepper = stepper
        self.settings = settings or AdaptiveSettings()
        self.context_provider = context_provider

    def integrate(
        self,
        stroke: StrokeModel,
        initial_state: CylinderState,
        theta_end: float,
        base_context: StrokeContext,
        extra_events: Sequence[EventBoundary] | None = None,
    ) -> IntegrationResult:
        """Integrate the stroke until ``theta_end``."""

        ctx_provider = self.context_provider or _constant_context_provider(base_context)
        rhs = lambda s, c: stroke.derivatives(s, c)
        boundaries = self._collect_boundaries(stroke, base_context, extra_events, theta_end)

        theta = initial_state.theta
        state = initial_state
        states = [initial_state]
        hits: list[EventBoundary] = []

        while theta < theta_end - self.settings.event_tolerance:
            context = ctx_provider(theta, state)
            derivs = rhs(state, context)
            step = self._propose_step(derivs)

            remaining = theta_end - theta
            step = min(step, remaining)

            upcoming = self._next_boundary(theta, boundaries)
            crossed_event = None
            if upcoming is not None:
                distance = upcoming.angle - theta
                if distance < self.settings.event_tolerance:
                    crossed_event = upcoming
                else:
                    step = min(step, distance)
                    if abs(step - distance) <= self.settings.event_tolerance:
                        crossed_event = upcoming

            state = self.stepper.step(rhs, state, step, ctx_provider)
            theta = state.theta
            states.append(state)

            if crossed_event is not None and (not hits or hits[-1] != crossed_event):
                hits.append(crossed_event)

        return IntegrationResult(states=states, events=hits)

    def _propose_step(self, derivs: CylinderDerivatives) -> float:
        settings = self.settings
        grad_pressure = abs(derivs.dP_dtheta) / max(settings.pressure_tolerance, 1e-12)
        grad_temp = abs(derivs.dT_dtheta) / max(settings.temperature_tolerance, 1e-12)
        gradient = max(grad_pressure, grad_temp, 1e-12)
        proposed = settings.safety_factor / gradient
        return _clamp(proposed, settings.min_step, settings.max_step)

    def _collect_boundaries(
        self,
        stroke: StrokeModel,
        base_context: StrokeContext,
        extra_events: Sequence[EventBoundary] | None,
        theta_end: float,
    ) -> list[EventBoundary]:
        boundaries: list[EventBoundary] = []

        if base_context.valve_timings is not None:
            timings: ValveTimingMap = base_context.valve_timings
            boundaries.extend(
                [
                    EventBoundary(angle=timings.ivo, name="intake-open"),
                    EventBoundary(angle=timings.ivc, name="intake-close"),
                    EventBoundary(angle=timings.evo, name="exhaust-open"),
                    EventBoundary(angle=timings.evc, name="exhaust-close"),
                ]
            )

        wiebe: WiebeHeatRelease | None = getattr(stroke, "wiebe", None)
        if isinstance(wiebe, WiebeHeatRelease):
            boundaries.append(EventBoundary(angle=wiebe.start_angle, name="ignition-start"))
            boundaries.append(
                EventBoundary(angle=wiebe.start_angle + wiebe.duration, name="ignition-end")
            )

        if extra_events:
            boundaries.extend(extra_events)

        upper = theta_end + self.settings.event_tolerance
        filtered = [b for b in boundaries if b.angle <= upper]
        filtered.sort(key=lambda b: b.angle)
        return filtered

    def _next_boundary(
        self, theta: float, boundaries: Sequence[EventBoundary]
    ) -> EventBoundary | None:
        for boundary in boundaries:
            if boundary.angle > theta + self.settings.event_tolerance:
                return boundary
        return None


def _advance_state(
    state: CylinderState,
    derivs: CylinderDerivatives,
    step: float,
    context: StrokeContext,
) -> CylinderState:
    return CylinderState(
        theta=state.theta + step,
        volume=state.volume + context.dvolume_dtheta * step,
        pressure=state.pressure + derivs.dP_dtheta * step,
        temperature=state.temperature + derivs.dT_dtheta * step,
        mass=state.mass + derivs.dm_dtheta * step,
    )


def _constant_context_provider(context: StrokeContext) -> Callable[[float, CylinderState], StrokeContext]:
    return lambda _theta, _state: context


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))
