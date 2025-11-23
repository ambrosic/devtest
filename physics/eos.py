"""Equation-of-state utilities.

The interfaces in this module are intentionally narrow so they can be
backed by ideal-gas relations for now and upgraded to real-gas models in
the future without affecting callers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


class EquationOfState(Protocol):
    """Minimal interface for thermodynamic properties.

    Implementations are expected to be pure (no internal caches or global
    state) so they can be swapped easily for alternative formulations.
    """

    name: str

    def pressure(self, density: float, temperature: float) -> float:
        """Return pressure [Pa] given density [kg/m^3] and temperature [K]."""

    def temperature(self, density: float, pressure: float) -> float:
        """Return temperature [K] given density [kg/m^3] and pressure [Pa]."""

    def internal_energy(self, temperature: float) -> float:
        """Specific internal energy [J/kg] as a function of temperature."""

    def temperature_from_energy(self, energy: float) -> float:
        """Inverse of :func:`internal_energy` for convenience."""

    def specific_heats(self) -> tuple[float, float]:
        """Return (cp, cv) [J/(kg*K)]."""

    def gamma(self) -> float:
        """Heat capacity ratio cp/cv (dimensionless)."""


@dataclass(frozen=True)
class GasProperties:
    """Constant property set for an ideal gas model."""

    gas_constant: float  # specific gas constant R [J/(kg*K)]
    gamma: float  # ratio of specific heats cp/cv

    @property
    def cp(self) -> float:
        return self.gamma * self.gas_constant / (self.gamma - 1.0)

    @property
    def cv(self) -> float:
        return self.gas_constant / (self.gamma - 1.0)


class IdealGasEOS:
    """Ideal-gas implementation of :class:`EquationOfState`."""

    def __init__(self, properties: GasProperties, name: str | None = None) -> None:
        self.properties = properties
        self.name = name or "ideal-gas"

    def pressure(self, density: float, temperature: float) -> float:
        return density * self.properties.gas_constant * temperature

    def temperature(self, density: float, pressure: float) -> float:
        if density <= 0:
            raise ValueError("Density must be positive to compute temperature")
        return pressure / (density * self.properties.gas_constant)

    def internal_energy(self, temperature: float) -> float:
        return self.properties.cv * temperature

    def temperature_from_energy(self, energy: float) -> float:
        return energy / self.properties.cv

    def specific_heats(self) -> tuple[float, float]:
        return (self.properties.cp, self.properties.cv)

    def gamma(self) -> float:
        return self.properties.gamma
