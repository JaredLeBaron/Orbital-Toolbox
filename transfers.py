
# Jared LeBaron
# Transfers

"""
Transfer utilities.

This module provides functions for computing idealized orbital transfer
requirements.

Assumptions:
- Position/radius units: km
- Velocity units: km/s
- Time units: seconds
- Two-body gravity
"""

import numpy as np


def hohmann_transfer_requirements(r1, r2, mu):
    """
    Compute ideal two-impulse Hohmann transfer requirements between two
    circular, coplanar orbits.

    Assumptions:
    - r1 and r2 are scalar orbital radii measured from the central body's center.
    - Initial and final orbits are circular and coplanar.
    - Burns are impulsive and tangential.
    - Units must be consistent, typically km, km/s, and seconds.

    Parameters
    ----------
    r1 : float
        Initial circular orbit radius [km].
    r2 : float
        Final circular orbit radius [km].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].

    Returns
    -------
    requirements : dict
        Dictionary containing circular velocities, transfer velocities,
        delta-v values, transfer time, and orbit periods.
    """

    if np.ndim(r1) != 0 or np.ndim(r2) != 0:
        raise ValueError("r1 and r2 must be scalar orbital radii, not vectors.")

    if np.ndim(mu) != 0:
        raise ValueError("mu must be a scalar gravitational parameter.")

    r1 = float(r1)
    r2 = float(r2)
    mu = float(mu)

    if r1 <= 0:
        raise ValueError("r1 must be greater than 0.")

    if r2 <= 0:
        raise ValueError("r2 must be greater than 0.")

    if mu <= 0:
        raise ValueError("mu must be greater than 0.")

    orbit_1_period = 2 * np.pi * np.sqrt(r1**3 / mu)
    orbit_2_period = 2 * np.pi * np.sqrt(r2**3 / mu)

    v_circular_1 = np.sqrt(mu / r1)
    v_circular_2 = np.sqrt(mu / r2)

    if r1 == r2:
        return {
            "r1": r1,
            "r2": r2,
            "mu": mu,
            "a_transfer": r1,
            "v_circular_1": v_circular_1,
            "v_circular_2": v_circular_2,
            "v_transfer_at_r1": v_circular_1,
            "v_transfer_at_r2": v_circular_2,
            "dv1": 0.0,
            "dv2": 0.0,
            "total_dv": 0.0,
            "transfer_time": 0.0,
            "orbit_1_period": orbit_1_period,
            "orbit_2_period": orbit_2_period,
            "transfer_type": "none",
            "notes": "Initial and final radii are equal; no Hohmann transfer required.",
        }

    a_transfer = (r1 + r2) / 2

    v_transfer_at_r1 = np.sqrt(mu * (2 / r1 - 1 / a_transfer))
    v_transfer_at_r2 = np.sqrt(mu * (2 / r2 - 1 / a_transfer))

    dv1 = v_transfer_at_r1 - v_circular_1
    dv2 = v_circular_2 - v_transfer_at_r2

    total_dv = abs(dv1) + abs(dv2)

    transfer_time = np.pi * np.sqrt(a_transfer**3 / mu)

    if r2 > r1:
        transfer_type = "raise"
    else:
        transfer_type = "lower"

    return {
        "r1": r1,
        "r2": r2,
        "mu": mu,
        "a_transfer": a_transfer,
        "v_circular_1": v_circular_1,
        "v_circular_2": v_circular_2,
        "v_transfer_at_r1": v_transfer_at_r1,
        "v_transfer_at_r2": v_transfer_at_r2,
        "dv1": dv1,
        "dv2": dv2,
        "total_dv": total_dv,
        "transfer_time": transfer_time,
        "orbit_1_period": orbit_1_period,
        "orbit_2_period": orbit_2_period,
        "transfer_type": transfer_type,
        "notes": "Ideal Hohmann transfer for circular, coplanar orbits.",
    }