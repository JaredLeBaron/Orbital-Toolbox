
# Jared LeBaron
# Dynamics

"""
Physical models for orbital motion.

This module defines acceleration models used in orbit propagation.
Currently includes the two-body gravitational model, which assumes
a point-mass central body.

Assumptions:
- Two-body central gravity
- Distance in km
- Time in seconds
- Acceleration in km/s^2
"""

import numpy as np


def two_body_acceleration(r_vec, mu):
    """
    Compute gravitational acceleration under the two-body assumption.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector relative to the central body [km].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].

    Returns
    -------
    a_vec : ndarray (3,)
        Acceleration vector [km/s^2].

    Notes
    -----
    - Uses Newton's law of gravitation:
          a_vec = -mu * r_vec / |r_vec|^3
    - Assumes the central body is fixed and spherically symmetric.
    """

    r_vec = np.array(r_vec, dtype=float)

    if r_vec.shape != (3,):
        raise ValueError("r_vec must be a 3-element position vector.")

    if np.ndim(mu) != 0:
        raise ValueError("mu must be a scalar gravitational parameter.")

    mu = float(mu)

    if mu <= 0:
        raise ValueError("mu must be greater than 0.")

    r_norm = np.linalg.norm(r_vec)

    if r_norm == 0:
        raise ValueError("Position magnitude is zero. Cannot compute gravitational acceleration.")

    return -mu * r_vec / r_norm**3