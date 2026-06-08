
# Jared LeBaron
# Orbital Elements

"""
Orbital element conversion utilities.

This module converts between Cartesian state vectors and classical orbital
elements. It also builds position/velocity vectors in the perifocal frame
and rotates them into the ECI frame.

Assumptions:
- Two-body orbital mechanics
- Distance in km
- Velocity in km/s
- Time in seconds
- Angles in radians internally, degrees included for output convenience
"""

import numpy as np

from rotations import perifocal_to_eci_matrix


def _as_vector3(vec, name):
    """
    Convert input to a 3-element NumPy vector.

    Parameters
    ----------
    vec : array_like
        Input vector.
    name : str
        Name used in error messages.

    Returns
    -------
    ndarray (3,)
        Vector as a NumPy array of floats.
    """

    vec = np.array(vec, dtype=float)

    if vec.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    return vec


def _as_scalar(value, name):
    """
    Convert input to a scalar float.

    Parameters
    ----------
    value : float
        Input value.
    name : str
        Name used in error messages.

    Returns
    -------
    float
        Scalar value as a float.
    """

    if np.ndim(value) != 0:
        raise ValueError(f"{name} must be a scalar.")

    return float(value)


def orbital_elements_from_rv(r_vec, v_vec, mu):
    """
    Compute classical orbital elements and related quantities from a Cartesian state.

    Given position and velocity vectors, this function computes the classical
    orbital elements along with auxiliary vectors and derived quantities.

    Parameters
    ----------
    r_vec : array_like (3,)
        Position vector [km].
    v_vec : array_like (3,)
        Velocity vector [km/s].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].

    Returns
    -------
    results : dict
        Dictionary containing position/velocity magnitudes, angular momentum,
        node vector, eccentricity vector, orbital type, specific energy,
        semi-major axis, periapsis/apoapsis distances, period, angular
        elements, ascending/descending node position and velocity vectors.

    Notes
    -----
    - Assumes a two-body central gravitational field.
    - Special handling is included for circular and equatorial edge cases.
    - For circular or equatorial orbits, some angular elements are undefined
      mathematically and are set to 0.0 by convention.
    """

    r_vec = _as_vector3(r_vec, "r_vec")
    v_vec = _as_vector3(v_vec, "v_vec")
    mu = _as_scalar(mu, "mu")

    if mu <= 0:
        raise ValueError("mu must be greater than 0.")

    tol = 1e-12

    i_hat = np.array([1.0, 0.0, 0.0])
    k_hat = np.array([0.0, 0.0, 1.0])

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    if r <= tol:
        raise ValueError("Position magnitude is zero. Cannot compute orbital elements.")


    # Angular Momentum
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)

    if h <= tol:
        raise ValueError("Angular momentum magnitude is zero. Orbit plane is undefined.")

    # Inclination
    i = np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0))

    # Node Vector
    n_vec = np.cross(k_hat, h_vec)
    n = np.linalg.norm(n_vec)

    # Right ascending of ascending node
    if n > tol:
        raan = np.arccos(np.clip(np.dot(n_vec, i_hat) / n, -1.0, 1.0))

        if n_vec[1] < 0:
            raan = 2 * np.pi - raan
    else:
        raan = 0.0

    # Eccentricity
    e_vec = (1 / mu) * ((v**2 - mu / r) * r_vec - np.dot(r_vec, v_vec) * v_vec)
    e = np.linalg.norm(e_vec)

    if e < 1e-10:
        orbital_type = "Circular"
    elif e < 1 - tol:
        orbital_type = "Elliptical"
    elif abs(e - 1) < tol:
        orbital_type = "Parabolic"
    else:
        orbital_type = "Hyperbolic"

    # Argument of periapsis
    if n > tol and e > tol:
        argp = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1.0, 1.0))

        if e_vec[2] < 0:
            argp = 2 * np.pi - argp
    else:
        argp = 0.0

    if e > tol:
        nu = np.arccos(np.clip(np.dot(e_vec, r_vec) / (e * r), -1.0, 1.0))

        if np.dot(r_vec, v_vec) < 0:
            nu = 2 * np.pi - nu
    else:
        nu = 0.0

    # Specific Energy
    eps = v**2 / 2 - mu / r

    # Semi-major axis
    if abs(eps) > tol:
        a = -mu / (2 * eps)
    else:
        a = np.inf

    if np.isfinite(a):
        rp = a * (1 - e)
    else:
        rp = np.nan

    # Period
    if orbital_type == "Elliptical":
        ra = a * (1 + e)
        T = 2 * np.pi * np.sqrt(a**3 / mu)
    else:
        ra = np.inf
        T = np.inf

    # Ascending and Descending Nodes
    #
    # The node vector points toward the ascending node.
    # The descending node is in the opposite direction.
    #
    # Notes:
    # - Nodes are undefined for equatorial orbits because n_vec = 0.
    # - For circular inclined orbits, argp is convention-based because periapsis is undefined.

    r_vec_asc = None
    r_vec_des = None
    v_vec_asc = None
    v_vec_des = None

    nu_asc = None
    nu_des = None

    if n > tol and np.isfinite(a):
        h_hat = h_vec / h

        # True anomaly at ascending and descending nodes
        nu_asc = (-argp) % (2 * np.pi)
        nu_des = (np.pi - argp) % (2 * np.pi)

        p = a * (1 - e ** 2)

        # Ascending node position
        r_asc = p / (1 + e * np.cos(nu_asc))
        r_vec_asc = r_asc * n_vec / n

        r_hat_asc = r_vec_asc / np.linalg.norm(r_vec_asc)
        theta_hat_asc = np.cross(h_hat, r_hat_asc)

        v_rad_asc = (mu / h) * e * np.sin(nu_asc)
        v_theta_asc = h / r_asc

        v_vec_asc = v_rad_asc * r_hat_asc + v_theta_asc * theta_hat_asc

        # Descending node position
        r_des = p / (1 + e * np.cos(nu_des))
        r_vec_des = -r_des * n_vec / n

        r_hat_des = r_vec_des / np.linalg.norm(r_vec_des)
        theta_hat_des = np.cross(h_hat, r_hat_des)

        v_rad_des = (mu / h) * e * np.sin(nu_des)
        v_theta_des = h / r_des

        v_vec_des = v_rad_des * r_hat_des + v_theta_des * theta_hat_des



    return {
        "r_vec": r_vec,
        "v_vec": v_vec,
        "r": r,
        "v": v,

        "h_vec": h_vec,
        "h": h,

        "n_vec": n_vec,
        "n": n,

        "e_vec": e_vec,
        "e": e,

        "orbital_type": orbital_type,

        "eps": eps,
        "energy": eps,

        "a": a,
        "rp": rp,
        "ra": ra,

        "T": T,
        "period": T,

        "i_rad": i,
        "i_deg": np.degrees(i),

        "RAAN_rad": raan,
        "RAAN_deg": np.degrees(raan),

        "omega_rad": argp,
        "omega_deg": np.degrees(argp),

        "argp_rad": argp,
        "argp_deg": np.degrees(argp),

        "nu_rad": nu,
        "nu_deg": np.degrees(nu),

        "nu_asc_rad": nu_asc,
        "nu_asc_deg": None if nu_asc is None else np.degrees(nu_asc),

        "nu_des_rad": nu_des,
        "nu_des_deg": None if nu_des is None else np.degrees(nu_des),

        "r_vec_asc": r_vec_asc,
        "r_vec_des": r_vec_des,

        "v_vec_asc": v_vec_asc,
        "v_vec_des": v_vec_des
    }


def coe_to_perifocal_rv(a, e, nu, mu):
    """
    Construct position and velocity vectors in the perifocal frame.

    Parameters
    ----------
    a : float
        Semi-major axis [km].
    e : float
        Eccentricity [-].
    nu : float
        True anomaly [rad].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].

    Returns
    -------
    r_pf : ndarray (3,)
        Position vector in the perifocal frame [km].
    v_pf : ndarray (3,)
        Velocity vector in the perifocal frame [km/s].

    Notes
    -----
    - The returned z-components are zero because the perifocal frame lies in
      the orbital plane.
    - Parabolic orbits are not supported by this semi-major-axis formulation.
    """

    a = _as_scalar(a, "a")
    e = _as_scalar(e, "e")
    nu = _as_scalar(nu, "nu")
    mu = _as_scalar(mu, "mu")

    if mu <= 0:
        raise ValueError("mu must be greater than 0.")

    if e < 0:
        raise ValueError("e must be greater than or equal to 0.")

    if abs(e - 1.0) < 1e-12:
        raise ValueError("Parabolic orbits are not supported by this formulation.")

    p = a * (1 - e**2)

    if p <= 0:
        raise ValueError("Semi-latus rectum p must be greater than 0. Check a and e.")

    denominator = 1 + e * np.cos(nu)

    if abs(denominator) < 1e-12:
        raise ValueError("Radius denominator is near zero. State is singular.")

    r = p / denominator

    r_pf = np.array([
        r * np.cos(nu),
        r * np.sin(nu),
        0.0,
    ])

    v_pf = np.sqrt(mu / p) * np.array([
        -np.sin(nu),
        e + np.cos(nu),
        0.0,
    ])

    return r_pf, v_pf


def perifocal_rv_to_eci(r_pf, v_pf, i, raan, argp):
    """
    Rotate position and velocity vectors from the perifocal frame to ECI.

    Parameters
    ----------
    r_pf : array_like (3,)
        Position vector in the perifocal frame [km].
    v_pf : array_like (3,)
        Velocity vector in the perifocal frame [km/s].
    i : float
        Inclination [rad].
    raan : float
        Right ascension of the ascending node [rad].
    argp : float
        Argument of periapsis [rad].

    Returns
    -------
    r_eci : ndarray (3,)
        Position vector in the ECI frame [km].
    v_eci : ndarray (3,)
        Velocity vector in the ECI frame [km/s].

    Notes
    -----
    - Uses the 3-1-3 rotation sequence:
      Rz(raan) @ Rx(i) @ Rz(argp).
    """

    r_pf = _as_vector3(r_pf, "r_pf")
    v_pf = _as_vector3(v_pf, "v_pf")

    i = _as_scalar(i, "i")
    raan = _as_scalar(raan, "raan")
    argp = _as_scalar(argp, "argp")

    q_rotation = perifocal_to_eci_matrix(raan, i, argp)

    r_eci = q_rotation @ r_pf
    v_eci = q_rotation @ v_pf

    return r_eci, v_eci


def coe_to_rv(a, e, i, raan, argp, nu, mu):
    """
    Convert classical orbital elements into ECI position and velocity vectors.

    This is a high-level wrapper that first constructs the state vector in the
    perifocal frame, then rotates that state into the inertial ECI frame.

    Parameters
    ----------
    a : float
        Semi-major axis [km].
    e : float
        Eccentricity [-].
    i : float
        Inclination [rad].
    raan : float
        Right ascension of the ascending node [rad].
    argp : float
        Argument of periapsis [rad].
    nu : float
        True anomaly [rad].
    mu : float
        Gravitational parameter of the central body [km^3/s^2].

    Returns
    -------
    r_eci : ndarray (3,)
        Position vector in the ECI frame [km].
    v_eci : ndarray (3,)
        Velocity vector in the ECI frame [km/s].

    Notes
    -----
    - Angles must be provided in radians.
    - Parabolic orbits are not supported by the current semi-major-axis
      formulation.
    """

    r_pf, v_pf = coe_to_perifocal_rv(a, e, nu, mu)
    r_eci, v_eci = perifocal_rv_to_eci(r_pf, v_pf, i, raan, argp)

    return r_eci, v_eci