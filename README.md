
# Orbital Toolbox

Orbital Toolbox is an educational Python project for basic two-body orbital mechanics calculations, validation, and visualization.

## Purpose

I built this project as part of my academic rebuild in aerospace engineering. The goal is to strengthen my understanding of orbital mechanics, math, physics, and programming by implementing core astrodynamics concepts myself and validating them feature by feature. 

## Current Capabilities

- COE to RV conversion
- RV to COE conversion
- two-body propagation
- Euler and RK4 integration
- trajectory plotting
- orbital diagnostics such as energy and angular momentum
- impulsive burns
- inclination-change maneuver
- Hohmann transfer calculations
- validation checks

## Assumptions

- two-body gravity
- point-mass central body
- impulsive burns
- fixed timestep propagation
- km, km/s, seconds
- radians internally, degrees for display/output

## Validation

The project includes `validation.py`, which runs sanity checks for the core toolbox functions. It currently checks COE/RV conversion, propagation behavior, burn behavior, Hohmann transfer behavior, node calculations, and inclination-change behavior.

Run from the project root:

```powershell
python validation.py
```

## Examples

Example scripts are located in the `examples/` folder.

Run from the project root:

```powershell
python examples/example_coe_2_rv.py
```
