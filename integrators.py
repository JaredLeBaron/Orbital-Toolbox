
# Jared LeBaron
# Integrators


"""
Numerical integration methods for orbit propagation.

This module provides one-step integration functions that advance a
Cartesian state (position and velocity) forward in time using a
specified timestep.

These functions are designed to be passed into higher-level propagation
wrappers (e.g., propagate_orbit).

Assumptions:
- State vectors are in km and km/s
- Time is in seconds
- Acceleration is provided by a dynamics model (e.g., two-body gravity)
"""

# Imports
from dynamics import two_body_acceleration


def euler_orbit_step(r_vec, v_vec, dt, mu):

    """
    Advance the state by one timestep using the forward Euler method.

    Parameters
    ----------
    r_vec : ndarray (3,)
        Current position vector [km].
    v_vec : ndarray (3,)
        Current velocity vector [km/s].
    dt : float
        Timestep [s].
    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    r_new : ndarray (3,)
        Updated position vector [km].
    v_new : ndarray (3,)
        Updated velocity vector [km/s].

    Notes
    -----
    - First-order method (low accuracy).
    - Energy is not conserved well over long propagation times.
    - Primarily useful for demonstration or comparison.
    """

    a_vec = two_body_acceleration(r_vec, mu)
    v_new = v_vec + a_vec * dt
    r_new = r_vec + v_new * dt
    return r_new, v_new


def rk4_orbit_step(r_vec, v_vec, dt, mu):

    """
    Advance the state by one timestep using fourth-order Runge-Kutta (RK4).

    Parameters
    ----------
    r_vec : ndarray (3,)
        Current position vector [km].
    v_vec : ndarray (3,)
        Current velocity vector [km/s].
    dt : float
        Timestep [s].
    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    r_new : ndarray (3,)
        Updated position vector [km].
    v_new : ndarray (3,)
        Updated velocity vector [km/s].

    Notes
    -----
    - Fourth-order method with significantly improved accuracy over Euler.
    - Provides good energy conservation for moderate timesteps.
    - Standard choice for orbit propagation in this toolbox.
    """

    # RK4 stages: evaluate slope (r', v') at intermediate points

    # k1
    a1 = two_body_acceleration(r_vec, mu)
    k1_r = dt * v_vec
    k1_v = dt * a1

    # k2
    r_k2 = r_vec + 0.5 * k1_r
    v_k2 = v_vec + 0.5 * k1_v
    a2 = two_body_acceleration(r_k2, mu)
    k2_r = dt * v_k2
    k2_v = dt * a2

    # k3
    r_k3 = r_vec + 0.5 * k2_r
    v_k3 = v_vec + 0.5 * k2_v
    a3 = two_body_acceleration(r_k3, mu)
    k3_r = dt * v_k3
    k3_v = dt * a3

    # k4
    r_k4 = r_vec + k3_r
    v_k4 = v_vec + k3_v
    a4 = two_body_acceleration(r_k4, mu)
    k4_r = dt * v_k4
    k4_v = dt * a4

    # Final weighted update
    r_new = r_vec + (1.0 / 6.0) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)
    v_new = v_vec + (1.0 / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

    return r_new, v_new