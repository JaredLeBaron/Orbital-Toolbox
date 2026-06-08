
# Jared LeBaron
# Propagator

"""
Propagation functions for two-body orbital motion.

This module runs orbital simulations and returns organized propagation
results including state histories, diagnostic quantities, and orbital
elements.

Assumptions:
- Position units: km
- Velocity units: km/s
- Time units: seconds
- Two-body gravity
"""

# Imports
import numpy as np

from orbital_elements import orbital_elements_from_rv
from integrators import rk4_orbit_step, euler_orbit_step


def initialize_history():
    """
    Initialize the history dictionary used during propagation.

    Returns
    -------
    history : dict
        Empty lists for position, velocity, time, energy, semi-major axis,
        eccentricity, and angular momentum history.
    """

    return {
        "x": [],
        "y": [],
        "z": [],
        "vx": [],
        "vy": [],
        "vz": [],
        "time": [],
        "energy": [],
        "a": [],
        "e": [],
        "h": [],
    }


def record_history(history, r_vec, v_vec, time_value, mu):
    """
    Record the current state and orbital diagnostics into the history dictionary.

    Parameters
    ----------
    history : dict
        Propagation history dictionary.
    r_vec : ndarray (3,)
        Current position vector [km].
    v_vec : ndarray (3,)
        Current velocity vector [km/s].
    time_value : float
        Current simulation time [s].
    mu : float
        Gravitational parameter [km^3/s^2].
    """

    elements = orbital_elements_from_rv(r_vec, v_vec, mu)

    history["x"].append(r_vec[0])
    history["y"].append(r_vec[1])
    history["z"].append(r_vec[2])
    history["vx"].append(v_vec[0])
    history["vy"].append(v_vec[1])
    history["vz"].append(v_vec[2])
    history["time"].append(time_value)
    history["energy"].append(elements["eps"])
    history["a"].append(elements["a"])
    history["e"].append(elements["e"])
    history["h"].append(elements["h"])


def propagate_orbit(r_vec, v_vec, dt, t_total, mu, step_function):
    """
    Propagate an orbital state forward in time using a specified numerical integrator.

    Parameters
    ----------
    r_vec : array_like (3,)
        Initial position vector [km].
    v_vec : array_like (3,)
        Initial velocity vector [km/s].
    dt : float
        Timestep size [s].
    t_total : float
        Total propagation time [s].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].
    step_function : callable
        Numerical integrator function, such as euler_orbit_step or rk4_orbit_step.

    Returns
    -------
    results : dict
        Dictionary containing initial state, final state, state histories,
        diagnostic histories, and initial/final orbital elements.

    Notes
    -----
    - This function returns one propagation result for one selected integrator.
    - Euler/RK4 comparison should be handled by a separate wrapper function.
    """

    r_initial = np.array(r_vec, dtype=float)
    v_initial = np.array(v_vec, dtype=float)

    if r_initial.shape != (3,):
        raise ValueError("r_vec must be a 3-element position vector.")

    if v_initial.shape != (3,):
        raise ValueError("v_vec must be a 3-element velocity vector.")

    if dt <= 0:
        raise ValueError("dt must be greater than 0.")

    if t_total <= 0:
        raise ValueError("t_total must be greater than 0.")

    if mu <= 0:
        raise ValueError("mu must be greater than 0.")

    r_current = r_initial.copy()
    v_current = v_initial.copy()

    initial_elements = orbital_elements_from_rv(r_initial, v_initial, mu)
    history = initialize_history()

    time_current = 0.0
    record_history(history, r_current, v_current, time_current, mu)

    while time_current < t_total - 1e-12:
        dt_step = min(dt, t_total - time_current)
        r_current, v_current = step_function(r_current, v_current, dt_step, mu)
        time_current += dt_step
        record_history(history, r_current, v_current, time_current, mu)

    final_elements = orbital_elements_from_rv(r_current, v_current, mu)

    return {
        "r_initial": r_initial,
        "v_initial": v_initial,
        "r_final": r_current,
        "v_final": v_current,
        "x": history["x"],
        "y": history["y"],
        "z": history["z"],
        "vx": history["vx"],
        "vy": history["vy"],
        "vz": history["vz"],
        "time": history["time"],
        "energy": history["energy"],
        "a": history["a"],
        "e": history["e"],
        "h": history["h"],
        "initial_elements": initial_elements,
        "final_elements": final_elements,
        "dt": dt,
        "t_total": t_total,
        "mu": mu,
        "step_function": step_function.__name__,
    }


def two_body_propagation_euler_rk4(r_vec, v_vec, dt, t_total, mu):
    """
    Propagate the same initial orbit using both Euler and RK4 integrators.

    Parameters
    ----------
    r_vec : array_like (3,)
        Initial position vector [km].
    v_vec : array_like (3,)
        Initial velocity vector [km/s].
    dt : float
        Timestep size [s].
    t_total : float
        Total propagation time [s].
    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    results : dict
        Dictionary containing full Euler and RK4 propagation results.

    Notes
    -----
    - This function is mainly for comparison and validation.
    - For normal propagation, use propagate_orbit directly.
    """

    initial_elements = orbital_elements_from_rv(
        np.array(r_vec, dtype=float),
        np.array(v_vec, dtype=float),
        mu,
    )

    euler_results = propagate_orbit(r_vec, v_vec, dt, t_total, mu, euler_orbit_step)
    rk4_results = propagate_orbit(r_vec, v_vec, dt, t_total, mu, rk4_orbit_step)

    return {
        "dt": dt,
        "t_total": t_total,
        "mu": mu,
        "initial_elements": initial_elements,
        "euler": euler_results,
        "rk4": rk4_results,
    }
