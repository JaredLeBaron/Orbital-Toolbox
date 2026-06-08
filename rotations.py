
# Jared LeBaron
# Rotations

"""
Rotation matrix utilities for coordinate transformations.

This module provides basic rotation matrices about principal axes
and constructs compound rotations used in orbital mechanics, such as
transforming from the perifocal frame to the ECI frame.

Assumptions:
- Angles are in radians
- Right-handed coordinate system
"""

import numpy as np


def _as_scalar_angle(theta, name):
    """
    Convert an input angle to a scalar float.

    Parameters
    ----------
    theta : float
        Input angle [rad].
    name : str
        Name used in error messages.

    Returns
    -------
    float
        Angle as a scalar float.
    """

    if np.ndim(theta) != 0:
        raise ValueError(f"{name} must be a scalar angle in radians.")

    return float(theta)


def rotation_matrix_x(theta):
    """
    Construct a rotation matrix for rotation about the x-axis.

    Parameters
    ----------
    theta : float
        Rotation angle [rad].

    Returns
    -------
    ndarray (3, 3)
        Rotation matrix.
    """

    theta = _as_scalar_angle(theta, "theta")

    return np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(theta), -np.sin(theta)],
        [0.0, np.sin(theta),  np.cos(theta)],
    ])


def rotation_matrix_y(theta):
    """
    Construct a rotation matrix for rotation about the y-axis.

    Parameters
    ----------
    theta : float
        Rotation angle [rad].

    Returns
    -------
    ndarray (3, 3)
        Rotation matrix.
    """

    theta = _as_scalar_angle(theta, "theta")

    return np.array([
        [ np.cos(theta), 0.0, np.sin(theta)],
        [0.0, 1.0, 0.0],
        [-np.sin(theta), 0.0, np.cos(theta)],
    ])


def rotation_matrix_z(theta):
    """
    Construct a rotation matrix for rotation about the z-axis.

    Parameters
    ----------
    theta : float
        Rotation angle [rad].

    Returns
    -------
    ndarray (3, 3)
        Rotation matrix.
    """

    theta = _as_scalar_angle(theta, "theta")

    return np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta),  np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])


def perifocal_to_eci_matrix(raan, i, argp):
    """
    Construct the rotation matrix from the perifocal frame to the ECI frame.

    This applies the 3-1-3 Euler rotation sequence:

        Rz(raan) @ Rx(i) @ Rz(argp)

    Parameters
    ----------
    raan : float
        Right ascension of ascending node [rad].
    i : float
        Inclination [rad].
    argp : float
        Argument of periapsis [rad].

    Returns
    -------
    ndarray (3, 3)
        Rotation matrix from perifocal frame to ECI frame.

    Notes
    -----
    - The transformation maps vectors from the orbital plane to the inertial frame.
    - Angles must be provided in radians.
    """

    raan = _as_scalar_angle(raan, "raan")
    i = _as_scalar_angle(i, "i")
    argp = _as_scalar_angle(argp, "argp")

    return rotation_matrix_z(raan) @ rotation_matrix_x(i) @ rotation_matrix_z(argp)
