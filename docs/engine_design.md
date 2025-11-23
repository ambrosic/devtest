# Thermodynamic Cycle Support and Simulation Architecture

## Supported Cycles and Required Outputs
- **Cycles:** Otto and Diesel cycles are supported initially, with the architecture structured to add future cycles (e.g., Miller/Atkinson) by extending the combustion/heat-release models.
- **Outputs:**
  - Pressure–Volume (P–V) diagram per cylinder and cycle.
  - Pressure–Crank Angle (P–θ) trace with resolved compression, combustion, expansion, and exhaust phases.
  - Aggregate power and torque at the crankshaft, including mean effective pressure (IMEP/FMEP/BMEP) derivations.
  - Cycle thermal efficiency (indicated and brake), fuel flow, and specific fuel consumption.
  - Diagnostic metrics: peak pressure/temperature, knock margin (where modeled), and mass balance closure.

## Stepping Strategy and Cylinder State
- **Stepping decision:** Use adaptive crank-angle stepping with embedded error control to capture rapid combustion heat-release while keeping cost low during smoother intake/exhaust periods. Clamp step sizes between a small minimum (e.g., 0.05°) and a larger cruise step (e.g., 1–2°) to avoid missed events and stiff instabilities. Offer an optional fixed-step mode for deterministic regression tests.
- **Cylinder state dataclass (conceptual):**
  ```python
  @dataclass
  class CylinderState:
      theta: float          # crank angle, radians
      volume: float         # instantaneous cylinder volume, m^3
      pressure: float       # cylinder pressure, Pa
      temperature: float    # bulk-gas temperature, K
      mass: float           # in-cylinder mass, kg
      gamma: float          # ratio of specific heats, updated with composition/temperature
      heat_release: float   # cumulative released heat, J
  ```
  State evolution couples mass/energy balances with heat-release and boundary heat-transfer models; gamma and temperature update via working-fluid property tables/ideal-gas relations.

## Architecture Overview
```mermaid
graph TD
  Config[Config Parsing]\n(YAML/TOML) --> Core
  Core[Physics Core]\n(cycle models, thermodynamics) --> Integrators
  Integrators[Integrators]\n(adaptive RK, events) --> State[State Manager]\n(cylinder state, caches)
  State --> Reports[Reporting]\n(power, torque, efficiency)
  State --> Viz[Visualization]\n(P-V, P-θ plots)
  Core --> Reports
  Config --> Reports
  Config --> Viz
  Reports --> Outputs[Artifacts]\n(CSV, JSON, PDF)
  Viz --> Outputs
```
- **Physics core:** Implements cycle-specific heat-release, gas exchange, friction, and property models; exposes RHS functions for integrators.
- **Integrators:** Adaptive crank-angle integrators with event hooks for valve timing, combustion start/stop, and volume changes; provide deterministic fixed-step mode for testing.
- **Config parsing:** Validates engine geometry, operating conditions, fuel properties, and solver tolerances; yields strongly typed configs for the core.
- **Reporting:** Aggregates traces into cycle metrics (IMEP/BMEP, power/torque, efficiencies) and exports tabular outputs.
- **Visualization:** Generates P–V and P–θ plots plus summary dashboards for comparisons across operating points or cycles.
