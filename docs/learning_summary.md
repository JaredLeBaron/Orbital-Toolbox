# Learning Summary

This project is part of an academic rebuild focused on taking ownership of
orbital mechanics fundamentals through code, validation, and documentation.

The goal is not only to produce working Python scripts. The goal is to explain
the assumptions, math, physics, numerical methods, validation checks, and
limitations behind each feature.

## Current Focus

The most recent mastery pass focused on the pure inclination-change helper in
`maneuvers.py`.

Key takeaways:

- An impulsive burn changes velocity instantly while position is held fixed.
- A pure plane change rotates the velocity vector rather than changing its
  magnitude.
- The orbital plane is defined by `r_vec` and `v_vec`; changing `v_vec` changes
  `h_vec = r_vec x v_vec`, which changes the plane orientation.
- The burn vector is `delta_v_vec = v_new - v_old`.
- The scalar plane-change check for the simple equal-speed case is
  `delta_v = 2 * v * sin(abs(delta_i) / 2)`.
- The helper assumes the burn state is physically appropriate, typically at a
  node. That assumption should be understood by the user.

## Validation Mindset

The validation suite is used as a learning and correctness tool. A passing
validation does not mean the toolbox is a production astrodynamics library; it
means the implemented assumptions are being checked against known relationships
and controlled cases.

Current validation coverage includes:

- COE to RV conversion checks
- COE to RV to COE round-trip checks
- ascending and descending node state checks
- RK4 one-orbit conservation checks
- impulsive burn behavior checks
- Hohmann transfer propagation checks
- inclination-change speed, delta-v, and final-inclination checks

## Ownership Goal

My goal is to work through each feature until I can explain it in:

- plain English
- equations
- inputs and outputs
- assumptions and limitations
- validation checks
- edge cases and failure modes

This is an ongoing process. Some parts of the toolbox are more mature than others, and the project is being reviewed
feature by feature as part of my academic rebuild.
