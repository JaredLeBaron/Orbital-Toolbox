
# Orbital Toolbox Documentation

## Overview

This toolbox provides tools for orbital mechanics including:

- COE ↔ RV conversion
- Orbit propagation
- Impulsive maneuvers
- Pure inclination-change burn calculation
- Transfer analysis, currently Hohmann transfer
- Validation and visualization

---

## Assumptions

- Two-body gravity
- Units:
  - Distance: km
  - Velocity: km/s
  - Time: seconds
  - Angles: radians internally, degrees for output/display

---

## File Structure

read_me.md  
Short project overview, quick-start instructions, current capabilities, and roadmap.

constants.py  
Stores physical constants such as gravitational parameters and body radii.

orbital_elements.py  
COE ↔ RV conversions and perifocal transformations, and ascending/descending node state calculations.

rotations.py  
Rotation matrices and frame transformations.

integrators.py  
Numerical integrators, currently Euler and RK4.

dynamics.py  
Acceleration models, currently two-body gravity.

propagate.py  
Orbit propagation and history tracking.

maneuvers.py  
Impulsive burns, direction vectors, and inclination-change calculations.

transfers.py  
Transfer calculations, currently Hohmann transfer.

diagnostics.py  
Error calculations and formatted summaries.

plotting.py  
Visualization utilities for single-orbit results and comparison results.

validation.py  
Verification tests for correctness. This file can be run directly as the toolbox health check.

docs/  
Technical documentation for the toolbox.

docs/Orbital_Toolbox_Docs.md  
Detailed documentation covering assumptions, equations, workflows, validation, limitations, and planned extensions.

examples/  
Polished demonstration scripts for toolbox features. Currently reserved for future examples.

scratch/  
Temporary experiments, throwaway tests, and sandbox scripts. Code in this folder is not considered part of the stable toolbox.

tests/  
Future automated test files.

---

## Core Workflows

### COE → RV

COE → perifocal state → ECI rotation → r_vec, v_vec

### Propagation

r_vec, v_vec → integrator → time history → final state

### Maneuver

direction → delta-v vector → new state → propagate

For a pure inclination-change calculation:

r_vec, v_vec, i_initial, i_final → rotated velocity vector → delta-v vector and magnitude

### Transfer

r1, r2 → transfer design → delta-v values → maneuver → propagation → validation

---

## Data Flow

The toolbox follows a consistent data flow:

COE → RV → Propagation → Diagnostics → Visualization

- Orbital elements are converted into Cartesian state vectors.
- State vectors are propagated over time using numerical integrators.
- Orbital quantities such as semi-major axis, eccentricity, angular momentum, and energy are recorded at each timestep.
- Results are analyzed, validated, and plotted.

---

## Design Principles

- Modular structure: physics, integration, maneuvers, transfers, diagnostics, and visualization are separated.
- Low-level functions remain flexible and reusable.
- Higher-level validation functions combine lower-level tools into complete checks.
- Interchangeable integrators allow Euler and RK4 comparison.
- Clear unit consistency: km, km/s, seconds.
- Functions are designed for reuse and extension.

---

## Key Functions

propagate_orbit  
Core propagation wrapper. Advances a state vector using one selected integrator and records history.

two_body_propagation_euler_rk4  
Comparison wrapper that propagates the same state using Euler and RK4.

rk4_orbit_step  
Fourth-order Runge-Kutta integrator used for accurate orbit propagation.

euler_orbit_step  
Forward Euler integrator used mainly for comparison.

two_body_acceleration  
Computes gravitational acceleration using the two-body model.

orbital_elements_from_rv  
Computes orbital elements from position and velocity vectors. Also returns ascending and descending node true anomalies,
position vectors, and velocity vectors when the nodes are defined.

coe_to_rv  
Converts classical orbital elements into Cartesian position and velocity vectors.

apply_impulsive_delta_v  
Applies an instantaneous velocity change to a state vector.

prograde_unit, radial_unit, normal_unit  
Generate common burn-direction unit vectors.

inclination_change  
Computes a pure plane-change burn by rotating the velocity vector about the current radius vector. Returns the delta-v vector, new velocity vector, delta-v magnitude, and scalar plane-change delta-v check.

hohmann_transfer_requirements  
Computes delta-v, transfer time, and related quantities for a Hohmann transfer.

run_all_validations  
Runs the full validation suite and reports pass/fail status.

---

## Current Capabilities

- COE ↔ RV conversion
- Ascending and descending node true anomaly, position vector, and velocity vector calculation
- Two-body orbit propagation
- Euler and RK4 integrators
- Impulsive maneuver modeling
- Direction-vector construction for prograde, retrograde, radial, inward radial, normal, and anti-normal directions
- Pure inclination-change burn calculation using velocity-vector rotation
- Hohmann transfer calculation and validation
- Validation of numerical accuracy
- Visualization of trajectories, histories, and error plots
- Central body sphere plotting using stored body radii

---

## Core Equations

These equations summarize the main physics and numerical methods used throughout the toolbox.

### Two-Body Gravitational Acceleration

a_vec = -mu * r_vec / |r_vec|^3

### Specific Angular Momentum

h_vec = r_vec × v_vec

h = |h_vec|

h_hat = h_vec / h

### Node Vector

n_vec = k_hat × h_vec

n = |n_vec|

The node vector points toward the ascending node.

n_hat = n_vec / n

The descending node points in the opposite direction:

r_hat_des = -n_hat

### Eccentricity Vector

e_vec = (1/mu) * ((v^2 - mu/r) * r_vec - (r_vec · v_vec) * v_vec)

e = |e_vec|

### Specific Orbital Energy

epsilon = v^2 / 2 - mu / r

### Semi-Major Axis

a = -mu / (2 * epsilon)

### Orbital Period

T = 2 * pi * sqrt(a^3 / mu)

### Perifocal Position

p = a * (1 - e^2)

r = p / (1 + e * cos(nu))

r_pf = [r cos(nu), r sin(nu), 0]

### Perifocal Velocity

v_pf = sqrt(mu / p) * [-sin(nu), e + cos(nu), 0]

### Ascending and Descending Nodes

The argument of latitude is:

u = omega + nu

At the ascending node:

u = 0

nu_asc = -omega

At the descending node:

u = pi

nu_des = pi - omega

The orbital radius at each node is found using:

r = p / (1 + e cos(nu))

The ascending node position vector is:

r_vec_asc = r_asc * n_hat

The descending node position vector is:

r_vec_des = -r_des * n_hat

The velocity vector at a node can be written using radial and transverse components:

v_vec = v_r * r_hat + v_theta * theta_hat

where:

r_hat = r_vec / |r_vec|

theta_hat = h_hat × r_hat

v_r = (mu / h) * e * sin(nu)

v_theta = h / r

### Euler Integration

v_{n+1} = v_n + a_n * dt

r_{n+1} = r_n + v_{n+1} * dt

### Runge-Kutta 4

y_{n+1} = y_n + (1/6)(k1 + 2k2 + 2k3 + k4)

### Pure Inclination Change

The inclination-change helper models an instantaneous plane-change burn at the current position.
It rotates the velocity vector about the radius direction using Rodrigues' rotation formula.

r_hat = r_vec / |r_vec|

delta_i = i_final - i_initial

v_new =
v cos(delta_i)
+ (r_hat x v) sin(delta_i)
+ r_hat (r_hat dot v) (1 - cos(delta_i))

delta_v_vec = v_new - v

delta_v_mag = |delta_v_vec|

For a pure plane change that preserves speed, the scalar check is:

delta_v = 2 v sin(|delta_i| / 2)

---

## Hohmann Transfer

The toolbox includes a Hohmann transfer module for computing and validating two-impulse transfers between circular, coplanar orbits.

### Assumptions

- Initial and final orbits are circular.
- Orbits are coplanar.
- Burns are impulsive.
- Burns are tangential.
- Inputs `r1` and `r2` are scalar orbital radii measured from the central body center.
- Inputs are radii, not altitudes.

---

### Method

The transfer consists of two burns:

1. Burn 1 at `r1`  
   Tangential burn to enter the transfer ellipse.

2. Burn 2 at `r2`  
   Tangential burn to circularize into the final orbit.

For a raising transfer, both burns are prograde.  
For a lowering transfer, both burns are retrograde.

---

### Hohmann Equations

Circular velocity:

v_c = sqrt(mu / r)

Transfer orbit semi-major axis:

a_t = (r1 + r2) / 2

Transfer velocity from the vis-viva equation:

v_t = sqrt(mu * (2/r - 1/a_t))

Delta-v:

dv1 = v_t(r1) - v_c(r1)

dv2 = v_c(r2) - v_t(r2)

Total delta-v:

dv_total = |dv1| + |dv2|

Transfer time:

t_transfer = pi * sqrt(a_t^3 / mu)

---

### Function Output

`hohmann_transfer_requirements()` returns a dictionary containing:

- `r1`
- `r2`
- `mu`
- `a_transfer`
- `v_circular_1`
- `v_circular_2`
- `v_transfer_at_r1`
- `v_transfer_at_r2`
- `dv1`
- `dv2`
- `total_dv`
- `transfer_time`
- `orbit_1_period`
- `orbit_2_period`
- `transfer_type`
- `notes`


## Orbital Elements Output

`orbital_elements_from_rv()` returns a dictionary containing:

- `r_vec`, `v_vec`
- `r`, `v`
- `h_vec`, `h`
- `n_vec`, `n`
- `e_vec`, `e`
- `orbital_type`
- `eps`, `energy`
- `a`
- `rp`, `ra`
- `T`, `period`
- `i_rad`, `i_deg`
- `RAAN_rad`, `RAAN_deg`
- `omega_rad`, `omega_deg`
- `argp_rad`, `argp_deg`
- `nu_rad`, `nu_deg`
- `nu_asc_rad`, `nu_asc_deg`
- `nu_des_rad`, `nu_des_deg`
- `r_vec_asc`
- `r_vec_des`
- `v_vec_asc`
- `v_vec_des`

When the nodes are undefined, the node-related outputs are set to `None`.

---

### Validation Approach

The transfer is validated using full propagation:

1. Start in circular orbit at `r1`.
2. Apply first burn `dv1`.
3. Propagate for the transfer time.
4. Apply second burn `dv2`.
5. Propagate the circularized final orbit.

The following are checked:

- Final transfer radius ≈ `r2`
- Final semi-major axis ≈ `r2`
- Final eccentricity ≈ 0

Example result:

PASS: Hohmann transfer validation successful.

---

### Limitations

- Only valid for circular, coplanar orbits.
- Does not handle:
  - elliptical initial or final orbits
  - combined transfer and inclination-change optimization
  - phasing or timing constraints
  - finite burn durations
  - perturbations

- Not a general transfer solver.
- Lambert solving is not yet implemented.
- Ascending and descending node states are undefined for equatorial orbits because the node vector has zero magnitude. 
For circular inclined orbits, the node direction is defined, but true anomaly relative to periapsis is convention-based because periapsis is undefined.

---

## Propagation Output Structure

`propagate_orbit()` returns a standardized dictionary containing:

- `r_initial`
- `v_initial`
- `r_final`
- `v_final`
- `x`, `y`, `z`
- `vx`, `vy`, `vz`
- `time`
- `energy`
- `a`
- `e`
- `h`
- `initial_elements`
- `final_elements`
- `dt`
- `t_total`
- `mu`
- `step_function`

This makes downstream plotting, diagnostics, and validation more predictable.

---

## Plotting System

The plotting module supports two main result types:

### Single-Orbit Plotting

These functions expect output directly from `propagate_orbit()`:

- `plot_2d_orbit`
- `plot_3d_orbit`
- `plot_history`
- `plot_standard_orbit_results`

### Comparison Plotting

These functions expect output from `two_body_propagation_euler_rk4()`:

- `plot_2d_comparison`
- `plot_3d_comparison`
- `plot_energy_comparison`
- `plot_semi_major_axis_comparison`
- `plot_eccentricity_comparison`
- `plot_angular_momentum_comparison`
- `plot_energy_error_comparison`
- `plot_angular_momentum_error_comparison`
- `plot_comparison_results`
- `plot_error_results`

The plotting module also includes:

- `plot_sphere`
- `set_axes_equal_3d`

These are used to show the central body and preserve accurate 3D visual scaling.

---

## Validation System

The toolbox includes a validation runner:

run_all_validations()

Current validation checks include:

- COE → RV conversion
- COE → RV → COE round trip
- RK4 one-orbit conservation
- impulsive burn behavior
- Hohmann transfer propagation and circularization

`validation.py` can be run directly from the project root.

On Windows PowerShell:

```powershell
py validation.py
```

## Limitations (Global)

- Two-body only
- No perturbations
- Fixed timestep integration
- Impulsive burn assumption
- No finite burn modeling
- No atmospheric drag
- No third-body gravity
- No Lambert solver yet
- No target phasing yet
- No time/date system yet

---

## Planned Extensions

- Inclination-change validation and combined plane-change transfer tools
- Bi-elliptic transfers
- J2 perturbation modeling
- Atmospheric drag
- Third-body gravity
- Adaptive timestep integrators
- Lambert solver
- Phasing and rendezvous tools
- Mission-planning wrapper functions
- Polished example scripts
- Automated test files
