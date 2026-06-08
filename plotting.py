
# Jared LeBaron
# Plotting

"""
Visualization utilities for orbital simulations.

This module provides plotting functions for visualizing orbital trajectories,
diagnostic histories, and comparison results between integrators.

Assumptions:
- Position in km
- Time in seconds
- Matplotlib used for visualization

Design:
- Functions ending in "_comparison" expect results from
  two_body_propagation_euler_rk4(), which contains "euler" and "rk4".
- Single-orbit plotting functions expect direct results from propagate_orbit().
"""

import matplotlib.pyplot as plt
import numpy as np

from diagnostics import angular_momentum_error, energy_error


# ============================================================
# Single propagation result plotters
# Expected input: output from propagate_orbit()
# ============================================================

def plot_2d_orbit(results, label="Orbit"):
    """
    Plot a single propagated orbit in the x-y plane.

    Parameters
    ----------
    results : dict
        Output dictionary from propagate_orbit.
    label : str
        Plot label.
    """

    plt.figure()
    plt.plot(results["x"], results["y"], label=label)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("2D Orbit Trajectory")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_3d_orbit(results, label="Orbit", show_body=False, body_radius=6378):
    """
    Plot a single propagated orbit in 3D.

    Parameters
    ----------
    results : dict
        Output dictionary from propagate_orbit.
    label : str
        Plot label.
    show_body : bool
        If True, plot a central body sphere.
    body_radius : float
        Radius of the central body [km].
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if show_body:
        plot_sphere(ax, radius=body_radius)

    ax.plot(results["x"], results["y"], results["z"], label=label)

    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_zlabel("z (km)")
    ax.set_title("3D Orbit Trajectory")
    ax.set_box_aspect([1, 1, 1])
    set_axes_equal_3d(ax)
    ax.legend()
    plt.show()


def plot_history(results, key, ylabel=None, title=None):
    """
    Plot a single recorded quantity versus time.

    Parameters
    ----------
    results : dict
        Output dictionary from propagate_orbit.
    key : str
        History key to plot, such as 'energy', 'a', 'e', or 'h'.
    ylabel : str, optional
        Y-axis label. If None, uses the key name.
    title : str, optional
        Plot title. If None, generates a default title.
    """

    if key not in results:
        raise KeyError(f"'{key}' not found in propagation results.")

    if ylabel is None:
        ylabel = key

    if title is None:
        title = f"{key} vs Time"

    plt.figure()
    plt.plot(results["time"], results[key])
    plt.xlabel("time (s)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.show()


def plot_standard_orbit_results(results, label="Orbit", show_body=False, body_radius=6378):
    """
    Generate standard plots for a single propagation result.

    Parameters
    ----------
    results : dict
        Output dictionary from propagate_orbit.
    label : str
        Plot label.
    show_body : bool
        If True, plot central body in the 3D trajectory plot.
    body_radius : float
        Radius of the central body [km].
    """

    plot_2d_orbit(results, label=label)
    plot_3d_orbit(results, label=label, show_body=show_body, body_radius=body_radius)
    plot_history(results, "energy", ylabel="Specific Orbital Energy", title="Specific Orbital Energy vs Time")
    plot_history(results, "a", ylabel="Semi-major Axis (km)", title="Semi-major Axis vs Time")
    plot_history(results, "e", ylabel="Eccentricity", title="Eccentricity vs Time")
    plot_history(results, "h", ylabel="Angular Momentum", title="Angular Momentum vs Time")


# ============================================================
# Euler vs RK4 comparison plotters
# Expected input: output from two_body_propagation_euler_rk4()
# ============================================================

def plot_2d_comparison(results):
    """
    Plot Euler and RK4 trajectories projected onto the x-y plane.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    plt.figure()
    plt.plot(euler["x"], euler["y"], label="Euler")
    plt.plot(rk4["x"], rk4["y"], label="RK4")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Orbit Trajectory Comparison (Euler vs RK4)")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_3d_comparison(results, show_body=False, body_radius=6378):
    """
    Plot Euler and RK4 trajectories in 3D.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.
    show_body : bool
        If True, plot central body.
    body_radius : float
        Radius of central body [km].
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if show_body:
        plot_sphere(ax, radius=body_radius)

    ax.plot(euler["x"], euler["y"], euler["z"], label="Euler")
    ax.plot(rk4["x"], rk4["y"], rk4["z"], label="RK4")

    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_zlabel("z (km)")
    ax.set_title("3D Orbit Trajectory Comparison (Euler vs RK4)")
    ax.set_box_aspect([1, 1, 1])
    set_axes_equal_3d(ax)
    ax.legend()
    plt.show()


def plot_energy_comparison(results):
    """
    Plot specific orbital energy over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    plt.figure()
    plt.plot(euler["time"], euler["energy"], label="Euler")
    plt.plot(rk4["time"], rk4["energy"], label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Specific Orbital Energy")
    plt.title("Specific Orbital Energy vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_semi_major_axis_comparison(results):
    """
    Plot semi-major axis over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    plt.figure()
    plt.plot(euler["time"], euler["a"], label="Euler")
    plt.plot(rk4["time"], rk4["a"], label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Semi-major axis (km)")
    plt.title("Semi-major Axis vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_eccentricity_comparison(results):
    """
    Plot eccentricity over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    plt.figure()
    plt.plot(euler["time"], euler["e"], label="Euler")
    plt.plot(rk4["time"], rk4["e"], label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Eccentricity")
    plt.title("Eccentricity vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_angular_momentum_comparison(results):
    """
    Plot angular momentum magnitude over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    plt.figure()
    plt.plot(euler["time"], euler["h"], label="Euler")
    plt.plot(rk4["time"], rk4["h"], label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Angular Momentum")
    plt.title("Angular Momentum vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_angular_momentum_error_comparison(results):
    """
    Plot relative angular momentum error over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    h_err_euler = angular_momentum_error(euler["h"])
    h_err_rk4 = angular_momentum_error(rk4["h"])

    plt.figure()
    plt.plot(euler["time"], h_err_euler, label="Euler")
    plt.plot(rk4["time"], h_err_rk4, label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Relative Angular Momentum Error")
    plt.title("Angular Momentum Error vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()


def plot_energy_error_comparison(results):
    """
    Plot relative energy error over time for Euler and RK4.
    """

    euler = results["euler"]
    rk4 = results["rk4"]

    e_err_euler = energy_error(euler["energy"])
    e_err_rk4 = energy_error(rk4["energy"])

    plt.figure()
    plt.plot(euler["time"], e_err_euler, label="Euler")
    plt.plot(rk4["time"], e_err_rk4, label="RK4")
    plt.xlabel("time (s)")
    plt.ylabel("Relative Energy Error")
    plt.title("Energy Error vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_error_results(results):
    """
    Generate standard error plots for Euler and RK4 propagation.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.
    """

    plot_energy_error_comparison(results)
    plot_angular_momentum_error_comparison(results)



def plot_comparison_results(results, show_body=False, body_radius=6378, show_errors=False):
    """
    Generate standard comparison plots for Euler and RK4 propagation.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.
    show_body : bool
        If True, plot central body in the 3D trajectory plot.
    body_radius : float
        Radius of the central body [km].
    show_errors : bool
        If True, also plot relative energy and angular momentum errors.
    """

    plot_2d_comparison(results)
    plot_3d_comparison(results, show_body=show_body, body_radius=body_radius)
    plot_energy_comparison(results)
    plot_semi_major_axis_comparison(results)
    plot_eccentricity_comparison(results)
    plot_angular_momentum_comparison(results)

    if show_errors:
        plot_energy_error_comparison(results)
        plot_angular_momentum_error_comparison(results)


# ============================================================
# 3D plotting helpers
# ============================================================

def plot_sphere(ax, radius=6378, center=(0, 0, 0), resolution=15, alpha=0.3, label="Earth"):
    """
    Plot a low-resolution sphere on a 3D axis.

    Parameters
    ----------
    ax : matplotlib 3D axis
        Axis to plot on.
    radius : float
        Sphere radius [km].
    center : tuple
        Sphere center coordinates [km].
    resolution : int
        Number of sample points used in each angular direction.
    alpha : float
        Surface transparency.
    label : str
        Label for legend.
    """

    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)

    x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = radius * np.outer(np.ones_like(u), np.cos(v)) + center[2]

    ax.plot_surface(x, y, z, alpha=alpha, label=label)


def set_axes_equal_3d(ax):
    """
    Set equal scaling for a 3D matplotlib axis.

    This prevents spheres from looking like ellipsoids and keeps orbital
    geometry visually accurate.
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    x_middle = sum(x_limits) / 2
    y_middle = sum(y_limits) / 2
    z_middle = sum(z_limits) / 2

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_vector(ax, origin, vector, scale=1, dimensions="2d", label=None):
    if dimensions == "3d":
        ax.quiver(
            origin[0], origin[1], origin[2],
            scale * vector[0], scale * vector[1], scale * vector[2],
            label=label
        )

    elif dimensions == "2d":
        ax.quiver(
            origin[0], origin[1],
            scale * vector[0], scale * vector[1],
            angles="xy",
            scale_units="xy",
            scale=1,
            label=label
        )

    else:
        raise ValueError("dimensions must be '2d' or '3d'")

    # TODO:
    # Improve vector overlay scaling.
    # Current issue: fixed visual scale does not generalize across orbit sizes.
    # Possible solution: normalized display vectors based on axis/orbit size.


# ============================================================
# Backward-compatible aliases
# These keep older code from breaking immediately.
# Eventually, update old code to use the newer names above.
# ============================================================

plot_2d_trajectory = plot_2d_comparison
plot_3d_trajectory = plot_3d_comparison
plot_energy = plot_energy_comparison
plot_semi_major_axis = plot_semi_major_axis_comparison
plot_eccentricity = plot_eccentricity_comparison
plot_angular_momentum = plot_angular_momentum_comparison
plot_angular_momentum_error = plot_angular_momentum_error_comparison
plot_energy_error = plot_energy_error_comparison
plot_results = plot_comparison_results
plot_error_results = plot_error_results