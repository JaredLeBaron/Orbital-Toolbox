
# Jared LeBaron
# Maneuvers

"""
Impulsive maneuver utilities.

This module provides tools for applying instantaneous velocity changes
(delta-v) to Cartesian state vectors and for constructing common burn
directions such as prograde, radial, and normal.

All direction functions return unit vectors in the ECI frame.

Assumptions:
- Position in km
- Velocity and delta-v in km/s
- Burns are impulsive, meaning position does not change during the burn
"""

import numpy as np


def _as_vector3(vec, name):
    """
    Convert input to a 3-element NumPy vector.

    Parameters
    ----------
    vec : array_like
        Input vector.
    name : str
        Name of the vector for error messages.

    Returns
    -------
    ndarray (3,)
        Vector as a NumPy array of floats.
    """

    vec = np.array(vec, dtype=float)

    if vec.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    return vec


def apply_impulsive_delta_v(r_vec, v_vec, delta_v_vec):
    """
    Apply an instantaneous velocity change to the current state.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector before the burn [km].
    v_vec : array_like (3,)
        Velocity vector before the burn [km/s].
    delta_v_vec : array_like (3,)
        Impulsive velocity change vector [km/s].

    Returns
    -------
    r_new : ndarray (3,)
        Position vector after the burn [km].
    v_new : ndarray (3,)
        Velocity vector after the burn [km/s].

    Notes
    -----
    - Position is unchanged because the burn is assumed instantaneous.
    """

    r_vec = _as_vector3(r_vec, "r_vec")
    v_vec = _as_vector3(v_vec, "v_vec")
    delta_v_vec = _as_vector3(delta_v_vec, "delta_v_vec")

    r_new = r_vec.copy()
    v_new = v_vec + delta_v_vec

    return r_new, v_new


def prograde_unit(v_vec):
    """
    Return the unit vector in the direction of velocity.

    Parameters
    ----------
    v_vec : array_like (3,)
        Velocity vector [km/s].

    Returns
    -------
    ndarray (3,)
        Unit vector in the prograde direction.
    """

    v_vec = _as_vector3(v_vec, "v_vec")
    v_mag = np.linalg.norm(v_vec)

    if v_mag == 0:
        raise ValueError("Velocity magnitude is zero. Cannot define prograde direction.")

    return v_vec / v_mag


def retrograde_unit(v_vec):
    """
    Return the unit vector opposite the direction of velocity.

    Parameters
    ----------
    v_vec : array_like (3,)
        Velocity vector [km/s].

    Returns
    -------
    ndarray (3,)
        Unit vector in the retrograde direction.
    """

    return -prograde_unit(v_vec)


def radial_unit(r_vec):
    """
    Return the unit vector in the radial direction from the central body.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector [km].

    Returns
    -------
    ndarray (3,)
        Unit vector in the outward radial direction.
    """

    r_vec = _as_vector3(r_vec, "r_vec")
    r_mag = np.linalg.norm(r_vec)

    if r_mag == 0:
        raise ValueError("Position magnitude is zero. Cannot define radial direction.")

    return r_vec / r_mag


def inward_radial_unit(r_vec):
    """
    Return the unit vector opposite the radial direction.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector [km].

    Returns
    -------
    ndarray (3,)
        Unit vector in the inward radial direction.
    """

    return -radial_unit(r_vec)


def normal_unit(r_vec, v_vec):
    """
    Return the unit vector normal to the orbital plane.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector [km].
    v_vec : array_like (3,)
        Velocity vector [km/s].

    Returns
    -------
    ndarray (3,)
        Unit vector perpendicular to the orbital plane.
    """

    r_vec = _as_vector3(r_vec, "r_vec")
    v_vec = _as_vector3(v_vec, "v_vec")

    h_vec = np.cross(r_vec, v_vec)
    h_mag = np.linalg.norm(h_vec)

    if h_mag == 0:
        raise ValueError("Angular momentum magnitude is zero. Cannot define normal direction.")

    return h_vec / h_mag


def anti_normal_unit(r_vec, v_vec):
    """
    Return the unit vector opposite the orbital angular momentum direction.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector [km].
    v_vec : array_like (3,)
        Velocity vector [km/s].

    Returns
    -------
    ndarray (3,)
        Unit vector opposite the normal direction.
    """

    return -normal_unit(r_vec, v_vec)


def delta_v_from_direction(delta_v_mag, direction_unit):
    """
    Construct a delta-v vector from a magnitude and direction.

    Parameters
    ----------
    delta_v_mag : float
        Signed velocity change magnitude [km/s].
    direction_unit : array_like (3,)
        Direction vector. Expected to be a unit vector.

    Returns
    -------
    ndarray (3,)
        Delta-v vector [km/s].

    Notes
    -----
    - If delta_v_mag is negative, the burn is applied opposite the given direction.
    """

    direction_unit = _as_vector3(direction_unit, "direction_unit")

    direction_mag = np.linalg.norm(direction_unit)

    if direction_mag == 0:
        raise ValueError("direction_unit magnitude is zero. Cannot construct delta-v.")

    # Normalize defensively so small user/input errors do not scale the burn incorrectly.
    direction_unit = direction_unit / direction_mag

    return delta_v_mag * direction_unit

def inclination_change(r_vec, v_vec, i_initial, i_final):
    """
    Performs a pure inclination-change burn by rotating the velocity vector
    about the current position vector.

    Inputs:
        r_vec     : position vector at burn location
        v_vec     : velocity vector before burn
        i_initial : initial inclination angle in radians
        i_final   : final inclination angle in radians

    Returns:
        delta_v_vec : required burn vector
        v_vec_new   : velocity vector after inclination change
        delta_v_mag : magnitude of required delta-v
        delta_v_scalar_check : delta v scalar check
    """

    r_vec = _as_vector3(r_vec, "r_vec")
    v_vec = _as_vector3(v_vec, "v_vec")

    r_mag = np.linalg.norm(r_vec)
    v_mag = np.linalg.norm(v_vec)

    if r_mag == 0:
        raise ValueError("r_vec magnitude cannot be zero.")

    r_hat = r_vec / r_mag

    i_delta = i_final - i_initial

    # Rodrigues' rotation formula:
    # rotates v_vec about r_hat by i_delta
    v_vec_new = (
        v_vec * np.cos(i_delta)
        + np.cross(r_hat, v_vec) * np.sin(i_delta)
        + r_hat * np.dot(r_hat, v_vec) * (1 - np.cos(i_delta))
    )

    delta_v_vec = v_vec_new - v_vec
    delta_v_mag = np.linalg.norm(delta_v_vec)

    # Optional scalar check:
    # For a pure plane change with same speed:
    # delta_v_mag should equal 2 * v_mag * sin(abs(i_delta)/2)
    delta_v_scalar_check = 2 * v_mag * np.sin(abs(i_delta) / 2)

    return delta_v_vec, v_vec_new, delta_v_mag, delta_v_scalar_check

