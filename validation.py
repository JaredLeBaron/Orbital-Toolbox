
# Jared LeBaron
# Validation

"""
Validation and sanity-check routines for the orbital toolbox.

This module runs numerical and physics-based checks for core toolbox
capabilities, including COE-to-RV conversion, round-trip element recovery,
RK4 propagation conservation, and impulsive maneuver behavior.

These routines are intended for development verification, not production use.

"""

import numpy as np

from constants import MU_BODIES, RADIUS_BODIES
from orbital_elements import coe_to_rv, orbital_elements_from_rv
from propagate import propagate_orbit
from integrators import rk4_orbit_step
from maneuvers import prograde_unit, normal_unit, delta_v_from_direction, apply_impulsive_delta_v, inclination_change
from transfers import hohmann_transfer_requirements


def print_scalar_check(name, original, recovered, units=""):
    """
    Print a formatted comparison between a reference scalar and computed scalar.
    """
    error = recovered - original
    print(
        f"{name:<16}: original = {original:.12f} {units} | "
        f"recovered = {recovered:.12f} {units} | "
        f"error = {error:.3e} {units}"
    )


def print_angle_check(name, original_rad, recovered_deg):
    """
    Print a formatted comparison between a reference angle and computed angle.
    """
    original_deg = np.rad2deg(original_rad)
    error = recovered_deg - original_deg
    print(
        f"{name:<16}: original = {original_deg:.6f} deg | "
        f"recovered = {recovered_deg:.6f} deg | "
        f"error = {error:.3e} deg"
    )


def coe_to_rv_validation():
    """
    Validate COE-to-RV conversion using radius and vis-viva speed checks.
    """
    mu = MU_BODIES["Earth"]

    a = 12000
    e = 0.3
    i = np.deg2rad(40)
    raan = np.deg2rad(30)
    argp = np.deg2rad(60)
    nu = np.deg2rad(45)

    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)

    r_mag = np.linalg.norm(r_vec)
    v_mag = np.linalg.norm(v_vec)

    p = a * (1 - e**2)
    r_expected = p / (1 + e * np.cos(nu))
    v_expected = np.sqrt(mu * (2 / r_expected - 1 / a))

    print("\nCOE → RV Validation")
    print("-------------------")
    print_scalar_check("r magnitude", r_expected, r_mag, "km")
    print_scalar_check("v magnitude", v_expected, v_mag, "km/s")

    r_error = abs(r_mag - r_expected)
    v_error = abs(v_mag - v_expected)

    r_pass = r_error < 1e-10
    v_pass = v_error < 1e-10

    validation_passed = r_pass and v_pass

    if validation_passed:
        print("PASS: COE to RV validation successful.")
    else:
        print("FAIL: COE to RV validation outside tolerance.")

    return validation_passed


def round_trip_validation():
    """
    Validate COE-to-RV-to-COE round-trip consistency.
    """
    mu = MU_BODIES["Earth"]

    a = 12000
    e = 0.3
    i = np.deg2rad(40)
    raan = np.deg2rad(30)
    argp = np.deg2rad(60)
    nu = np.deg2rad(45)

    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)
    elements = orbital_elements_from_rv(r_vec, v_vec, mu)

    print("\nCOE → RV → COE Validation")
    print("-------------------------")
    print_scalar_check("a", a, elements["a"], "km")
    print_scalar_check("e", e, elements["e"])
    print_angle_check("i", i, elements["i_deg"])
    print_angle_check("RAAN", raan, elements["RAAN_deg"])
    print_angle_check("argp", argp, elements["omega_deg"])
    print_angle_check("nu", nu, elements["nu_deg"])

    a_pass = abs(elements["a"] - a) < 1e-8
    e_pass = abs(elements["e"] - e) < 1e-12
    i_pass = abs(elements["i_rad"] - i) < 1e-12
    raan_pass = abs(elements["RAAN_rad"] - raan) < 1e-12
    argp_pass = abs(elements["argp_rad"] - argp) < 1e-12
    nu_pass = abs(elements["nu_rad"] - nu) < 1e-12

    validation_passed = (
        a_pass and
        e_pass and
        i_pass and
        raan_pass and
        argp_pass and
        nu_pass
    )

    if validation_passed:
        print("PASS: COE round-trip validation successful.")
    else:
        print("FAIL: COE round-trip validation outside tolerance.")

    return validation_passed


def rk4_propagation_validation():
    """
    Validate RK4 propagation over one orbit by checking conservation of
    semi-major axis, eccentricity, energy, and angular momentum.
    """
    mu = MU_BODIES["Earth"]

    a = 12000
    e = 0.3
    i = np.deg2rad(40)
    raan = np.deg2rad(30)
    argp = np.deg2rad(60)
    nu = np.deg2rad(45)

    dt = 10
    T = 2 * np.pi * np.sqrt(a**3 / mu)
    t_total = T

    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)

    rk4_results = propagate_orbit(
        r_vec,
        v_vec,
        dt,
        t_total,
        mu,
        rk4_orbit_step
    )

    initial = orbital_elements_from_rv(r_vec, v_vec, mu)
    final = rk4_results["final_elements"]

    print("\nRK4 One-Orbit Propagation Validation")
    print("------------------------------------")
    print_scalar_check("a", initial["a"], final["a"], "km")
    print_scalar_check("e", initial["e"], final["e"])
    print_scalar_check("eps", initial["eps"], final["eps"], "km^2/s^2")
    print_scalar_check("h", initial["h"], final["h"], "km^2/s")

    print(f"\nPeriod: {T:.6f} s")
    print(f"Time propagated: {t_total:.6f} s")


    a_pass = abs(final["a"] - initial["a"]) < 1e-4
    e_pass = abs(final["e"] - initial["e"]) < 1e-8
    eps_pass = abs(final["eps"] - initial["eps"]) < 1e-8
    h_pass = abs(final["h"] - initial["h"]) < 1e-4

    print("PASS: semi-major axis conserved" if a_pass else "FAIL: semi-major axis drift too large")
    print("PASS: eccentricity conserved" if e_pass else "FAIL: eccentricity drift too large")
    print("PASS: energy conserved" if eps_pass else "FAIL: energy drift too large")
    print("PASS: angular momentum conserved" if h_pass else "FAIL: angular momentum drift too large")

    validation_passed = a_pass and e_pass and eps_pass and h_pass

    if validation_passed:
        print("PASS: RK4 propagation validation successful.")
    else:
        print("FAIL: RK4 propagation validation outside tolerance.")

    return validation_passed


def burn_validation():
    """
    Validate impulsive maneuver behavior and post-burn propagation consistency.

    Checks:
    - Prograde burn increases semi-major axis, energy, and angular momentum.
    - Retrograde burn decreases semi-major axis, energy, and angular momentum.
    - Normal burn increases inclination.
    - Post-burn propagation conserves the new orbit's a, e, eps, and h.
    """

    mu = MU_BODIES["Earth"]

    a = 7000
    e = 0
    i = 0
    raan = 0
    argp = 0
    nu = 0

    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)
    elements_before = orbital_elements_from_rv(r_vec, v_vec, mu)

    print("\nBurn Validation")
    print("---------------")

    all_burns_passed = True

    for name in ["Prograde", "Retrograde", "Normal"]:

        if name == "Prograde":
            direction = prograde_unit(v_vec)
            dv = 0.1

        elif name == "Retrograde":
            direction = -prograde_unit(v_vec)
            dv = 0.1

        elif name == "Normal":
            direction = normal_unit(r_vec, v_vec)
            dv = 0.1

        delta_v_vec = delta_v_from_direction(dv, direction)

        r_new, v_new = apply_impulsive_delta_v(r_vec, v_vec, delta_v_vec)
        elements_after = orbital_elements_from_rv(r_new, v_new, mu)

        if name == "Prograde":
            behavior_pass = (
                elements_after["a"] > elements_before["a"] and
                elements_after["eps"] > elements_before["eps"] and
                elements_after["h"] > elements_before["h"]
            )

        elif name == "Retrograde":
            behavior_pass = (
                elements_after["a"] < elements_before["a"] and
                elements_after["eps"] < elements_before["eps"] and
                elements_after["h"] < elements_before["h"]
            )

        elif name == "Normal":
            behavior_pass = elements_after["i_rad"] > elements_before["i_rad"]

        print(f"\n{name} Burn")
        print("------------------------")
        print_scalar_check("a", elements_before["a"], elements_after["a"], "km")
        print_scalar_check("e", elements_before["e"], elements_after["e"])
        print_scalar_check("eps", elements_before["eps"], elements_after["eps"], "km^2/s^2")
        print_scalar_check("h", elements_before["h"], elements_after["h"], "km^2/s")
        print_angle_check("i", elements_before["i_rad"], elements_after["i_deg"])
        print_angle_check("RAAN", elements_before["RAAN_rad"], elements_after["RAAN_deg"])
        print_angle_check("argp", elements_before["omega_rad"], elements_after["omega_deg"])
        print_angle_check("nu", elements_before["nu_rad"], elements_after["nu_deg"])

        T_new = 2 * np.pi * np.sqrt(elements_after["a"] ** 3 / mu)

        prop_results = propagate_orbit(
            r_new,
            v_new,
            dt=10,
            t_total=T_new,
            mu=mu,
            step_function=rk4_orbit_step,
        )

        final = prop_results["final_elements"]

        a_conserved = abs(final["a"] - elements_after["a"]) < 1e-4
        e_conserved = abs(final["e"] - elements_after["e"]) < 1e-8
        eps_conserved = abs(final["eps"] - elements_after["eps"]) < 1e-8
        h_conserved = abs(final["h"] - elements_after["h"]) < 1e-4

        propagation_pass = (
            a_conserved and
            e_conserved and
            eps_conserved and
            h_conserved
        )

        burn_passed = behavior_pass and propagation_pass
        all_burns_passed = all_burns_passed and burn_passed

        print(f"\n{name} behavior check: {'PASS' if behavior_pass else 'FAIL'}")
        print(f"{name} post-burn propagation check: {'PASS' if propagation_pass else 'FAIL'}")

        print("\nPost-Burn Propagation Check")
        print("---------------------------")
        print_scalar_check("a", elements_after["a"], final["a"], "km")
        print_scalar_check("e", elements_after["e"], final["e"])
        print_scalar_check("eps", elements_after["eps"], final["eps"], "km^2/s^2")
        print_scalar_check("h", elements_after["h"], final["h"], "km^2/s")

    if all_burns_passed:
        print("\nPASS: burn validation successful.")
    else:
        print("\nFAIL: burn validation outside tolerance.")

    return all_burns_passed



def hohmann_transfer_requirements_validation(plot=True):
    """
    Validate the Hohmann transfer implementation by:

    1. Creating an initial circular orbit.
    2. Applying the first prograde burn.
    3. Propagating for the transfer time.
    4. Applying the circularization burn.
    5. Checking that the final orbit matches the target circular orbit.

    This validation assumes circular, coplanar orbits and uses RK4 propagation.

    Parameters
    ----------
    plot : bool
        If True, plot the initial orbit, transfer orbit, and target final orbit.

    Returns
    -------
    validation_passed : bool
        True if all validation checks pass.
    """

    mu = MU_BODIES["Earth"]

    r1 = 7000
    r2 = 8000
    dt = 1  # Small timestep used for validation accuracy

    reqs = hohmann_transfer_requirements(r1, r2, mu)

    v1 = reqs["v_circular_1"]
    v2 = reqs["v_circular_2"]
    dv1 = reqs["dv1"]
    dv2 = reqs["dv2"]
    total_dv = reqs["total_dv"]
    time_transfer = reqs["transfer_time"]
    orbit_1_period = reqs["orbit_1_period"]
    orbit_2_period = reqs["orbit_2_period"]

    # Initial circular orbit starts at +x and moves +y
    r1_vec = np.array([r1, 0.0, 0.0])
    v1_vec = np.array([0.0, v1, 0.0])

    # After half the Hohmann transfer, spacecraft reaches -x side
    r2_vec = np.array([-r2, 0.0, 0.0])
    v2_vec = np.array([0.0, -v2, 0.0])

    # Reference circular orbits
    initial_results = propagate_orbit(
        r1_vec,
        v1_vec,
        dt,
        orbit_1_period,
        mu,
        rk4_orbit_step,
    )

    final_results = propagate_orbit(
        r2_vec,
        v2_vec,
        dt,
        orbit_2_period,
        mu,
        rk4_orbit_step,
    )

    # Burn 1: enter transfer orbit
    burn1_vec = np.array([0.0, dv1, 0.0])

    r_after_burn1, v_after_burn1 = apply_impulsive_delta_v(
        r1_vec,
        v1_vec,
        burn1_vec,
    )

    # Propagate transfer orbit
    transfer_results = propagate_orbit(
        r_after_burn1,
        v_after_burn1,
        dt,
        time_transfer,
        mu,
        rk4_orbit_step,
    )

    # State at end of transfer
    r_transfer_final = transfer_results["r_final"]
    v_transfer_final = transfer_results["v_final"]

    # Burn 2: circularize at the final orbit
    # At -x, prograde direction is -y
    burn2_vec = np.array([0.0, -dv2, 0.0])

    r_after_burn2, v_after_burn2 = apply_impulsive_delta_v(
        r_transfer_final,
        v_transfer_final,
        burn2_vec,
    )

    circularized_results = propagate_orbit(
        r_after_burn2,
        v_after_burn2,
        dt,
        orbit_2_period,
        mu,
        rk4_orbit_step,
    )

    if plot:
        import matplotlib.pyplot as plt
        from plotting import plot_sphere, set_axes_equal_3d

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        plot_sphere(ax, radius=RADIUS_BODIES["Earth"])

        ax.plot(
            initial_results["x"],
            initial_results["y"],
            initial_results["z"],
            label="Initial Orbit",
        )

        ax.plot(
            transfer_results["x"],
            transfer_results["y"],
            transfer_results["z"],
            label="Transfer Orbit",
        )

        ax.plot(
            final_results["x"],
            final_results["y"],
            final_results["z"],
            label="Target Final Orbit",
        )

        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_zlabel("z (km)")
        ax.set_title("Hohmann Transfer Validation")
        ax.set_box_aspect([1, 1, 1])
        ax.legend()
        set_axes_equal_3d(ax)
        plt.show()

    # Validation quantities
    transfer_final_radius = np.linalg.norm(r_transfer_final)
    final_elements = circularized_results["final_elements"]

    radius_error = transfer_final_radius - r2
    a_error = final_elements["a"] - r2
    final_e = final_elements["e"]

    radius_tol = 1e-5
    a_tol = 1e-5
    e_tol = 1e-5

    radius_pass = abs(radius_error) < radius_tol
    a_pass = abs(a_error) < a_tol
    e_pass = final_e < e_tol

    validation_passed = radius_pass and a_pass and e_pass

    print("----- Hohmann Transfer Validation -----")
    print(f"dv1: {dv1:.8f} km/s")
    print(f"dv2: {dv2:.8f} km/s")
    print(f"total dv: {total_dv:.8f} km/s")
    print(f"transfer time: {time_transfer:.8f} s")
    print()
    print(f"Transfer final radius: {transfer_final_radius:.8f} km")
    print(f"Target radius: {r2:.8f} km")
    print(f"Radius error: {radius_error:.8e} km")
    print(f"Radius check: {'PASS' if radius_pass else 'FAIL'}")
    print()
    print(f"Final semi-major axis: {final_elements['a']:.8f} km")
    print(f"Semi-major axis error: {a_error:.8e} km")
    print(f"Semi-major axis check: {'PASS' if a_pass else 'FAIL'}")
    print()
    print(f"Final eccentricity: {final_e:.8e}")
    print(f"Eccentricity check: {'PASS' if e_pass else 'FAIL'}")
    print()

    if validation_passed:
        print("PASS: Hohmann transfer validation successful.")
    else:
        print("FAIL: Hohmann transfer validation outside tolerance.")

    return validation_passed


def node_state_validation():
    """
    Validate ascending and descending node state calculations.

    This checks that the node states returned by orbital_elements_from_rv()
    match independent COE -> RV calculations at the returned node true anomalies.
    """

    mu = MU_BODIES["Earth"]

    # Use an inclined, non-circular orbit.
    # Avoid circular/equatorial edge cases for this validation.
    a = 12000.0
    e = 0.3
    i = np.radians(40.0)
    raan = np.radians(30.0)
    argp = np.radians(60.0)
    nu = np.radians(45.0)

    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)

    elements = orbital_elements_from_rv(r_vec, v_vec, mu)

    nu_asc = elements["nu_asc_rad"]
    nu_des = elements["nu_des_rad"]

    r_vec_asc = elements["r_vec_asc"]
    v_vec_asc = elements["v_vec_asc"]

    r_vec_des = elements["r_vec_des"]
    v_vec_des = elements["v_vec_des"]

    if (
        nu_asc is None
        or nu_des is None
        or r_vec_asc is None
        or v_vec_asc is None
        or r_vec_des is None
        or v_vec_des is None
    ):
        print("FAIL: Node states were not computed.")
        return False

    r_check_asc, v_check_asc = coe_to_rv(a, e, i, raan, argp, nu_asc, mu)
    r_check_des, v_check_des = coe_to_rv(a, e, i, raan, argp, nu_des, mu)

    r_error_asc = np.linalg.norm(r_vec_asc - r_check_asc)
    v_error_asc = np.linalg.norm(v_vec_asc - v_check_asc)

    r_error_des = np.linalg.norm(r_vec_des - r_check_des)
    v_error_des = np.linalg.norm(v_vec_des - v_check_des)

    tolerance_r = 1e-6      # km
    tolerance_v = 1e-9      # km/s

    asc_pass = r_error_asc < tolerance_r and v_error_asc < tolerance_v
    des_pass = r_error_des < tolerance_r and v_error_des < tolerance_v

    if asc_pass and des_pass:
        return True

    print("FAIL: Node state validation failed.")
    print(f"Ascending r error [km]: {r_error_asc:.12e}")
    print(f"Ascending v error [km/s]: {v_error_asc:.12e}")
    print(f"Descending r error [km]: {r_error_des:.12e}")
    print(f"Descending v error [km/s]: {v_error_des:.12e}")

    return False



def inclination_change_validation():
    """
    Validate pure inclination-change behavior.

    Checks that the rotated velocity preserves speed, produces the expected
    plane-change delta-v, and reaches the target inclination.
    """

    mu = MU_BODIES["Earth"]

    a = 7000.0
    e = 0.0
    i_initial = 0.0
    raan = 0.0
    argp = 0.0
    nu = 0.0

    i_final = np.deg2rad(45.0)

    r_vec, v_vec = coe_to_rv(a, e, i_initial, raan, argp, nu, mu)

    delta_v_vec, v_vec_new, delta_v_mag, delta_v_scalar_check = inclination_change(
        r_vec,
        v_vec,
        i_initial,
        i_final,
    )

    elements_after = orbital_elements_from_rv(r_vec, v_vec_new, mu)

    v_initial = np.linalg.norm(v_vec)
    v_final = np.linalg.norm(v_vec_new)

    delta_i = i_final - i_initial
    delta_v_expected = 2 * v_initial * np.sin(abs(delta_i) / 2)

    speed_pass = np.isclose(v_initial, v_final, atol=1e-10)
    delta_v_pass = np.isclose(delta_v_mag, delta_v_expected, atol=1e-10)
    scalar_check_pass = np.isclose(delta_v_mag, delta_v_scalar_check, atol=1e-10)
    inclination_pass = np.isclose(elements_after["i_rad"], i_final, atol=1e-10)

    validation_passed = (
        speed_pass and
        delta_v_pass and
        scalar_check_pass and
        inclination_pass
    )

    print("\nInclination Change Validation")
    print("-----------------------------")
    print(f"Initial inclination : {np.degrees(i_initial):.6f} deg")
    print(f"Target inclination  : {np.degrees(i_final):.6f} deg")
    print(f"Recovered inclination: {elements_after['i_deg']:.6f} deg")
    print(f"Speed preserved     : {'PASS' if speed_pass else 'FAIL'}")
    print(f"Delta-v equation    : {'PASS' if delta_v_pass else 'FAIL'}")
    print(f"Scalar check        : {'PASS' if scalar_check_pass else 'FAIL'}")
    print(f"Inclination target  : {'PASS' if inclination_pass else 'FAIL'}")

    if validation_passed:
        print("PASS: inclination-change validation successful.")
    else:
        print("FAIL: inclination-change validation outside tolerance.")

    return validation_passed






def run_all_validations():
    results = {
        "coe_to_rv": coe_to_rv_validation(),
        "round_trip": round_trip_validation(),
        "nodes": node_state_validation(),
        "rk4": rk4_propagation_validation(),
        "burns": burn_validation(),
        "hohmann": hohmann_transfer_requirements_validation(plot=False),
        "inclination_change": inclination_change_validation(),
    }

    print("\n==========================")
    print("Validation Summary")
    print("==========================")

    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"{name:<15}: {status}")

    all_passed = all(results.values())

    print("--------------------------")
    print(f"Overall        : {'PASS' if all_passed else 'FAIL'}")

    return results


if __name__ == "__main__":

    run_all_validations()
