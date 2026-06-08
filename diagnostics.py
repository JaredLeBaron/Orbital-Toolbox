
# Jared LeBaron
# Diagnostics

"""
Utility functions for inspecting and summarizing orbital propagation results.

This module computes relative error histories, builds method comparison
summaries, and prints formatted orbital diagnostics. It is used mainly for
debugging, validation, and interpreting numerical propagation behavior.

Assumptions:
- Inputs come from the toolbox propagation and orbital element functions.
- Distances are in km.
- Velocities are in km/s.
- Time is in seconds.
"""


def relative_error(values):
    """
    Compute relative error over time compared to the initial value.

    Parameters
    ----------
    values : sequence of float
        Quantity history.

    Returns
    -------
    list of float
        Relative error history compared to the initial value.
    """

    initial_value = values[0]

    if initial_value == 0:
        raise ValueError("Initial value is zero. Cannot compute relative error.")

    return [(value - initial_value) / initial_value for value in values]


def angular_momentum_error(h_values):
    """
    Compute relative angular momentum error over time.

    Parameters
    ----------
    h_values : sequence of float
        Angular momentum magnitude history [km^2/s].

    Returns
    -------
    list of float
        Relative error history compared to the initial angular momentum.
    """

    return relative_error(h_values)


def energy_error(energy_values):
    """
    Compute relative specific energy error over time.

    Parameters
    ----------
    energy_values : sequence of float
        Specific orbital energy history [km^2/s^2].

    Returns
    -------
    list of float
        Relative error history compared to the initial specific energy.
    """

    return relative_error(energy_values)


def method_summary(result):
    """
    Build a summary dictionary comparing initial and final orbital quantities
    for one propagation result.

    Parameters
    ----------
    result : dict
        Output dictionary from propagate_orbit.

    Returns
    -------
    summary : dict
        Summary of initial, final, and change values for key orbital quantities.
    """

    initial = result["initial_elements"]
    final = result["final_elements"]

    return {
        "r_initial": result["r_initial"],
        "v_initial": result["v_initial"],
        "r_final": result["r_final"],
        "v_final": result["v_final"],

        "initial_a": initial["a"],
        "initial_e": initial["e"],
        "initial_eps": initial["eps"],
        "initial_h": initial["h"],

        "final_a": final["a"],
        "final_e": final["e"],
        "final_eps": final["eps"],
        "final_h": final["h"],

        "delta_a": final["a"] - initial["a"],
        "delta_e": final["e"] - initial["e"],
        "delta_eps": final["eps"] - initial["eps"],
        "delta_h": final["h"] - initial["h"],
    }


def build_comparison_summary(results):
    """
    Build Euler and RK4 method summaries from propagation comparison results.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.

    Returns
    -------
    summary : dict
        Dictionary containing summary data for Euler and RK4 propagation.
    """

    return {
        "euler": method_summary(results["euler"]),
        "rk4": method_summary(results["rk4"]),
    }


def print_method_summary(method_name, summary):
    """
    Print a formatted summary for one propagation method.

    Parameters
    ----------
    method_name : str
        Name of the propagation method.
    summary : dict
        Summary dictionary produced by method_summary.
    """

    print(f"\nUsing {method_name} Method")
    print("--------------------------")
    print("Initial r_vec =", summary["r_initial"])
    print("Initial v_vec =", summary["v_initial"])
    print("Final r_vec   =", summary["r_final"])
    print("Final v_vec   =", summary["v_final"])

    print("\nInitial quantities:")
    print(f"a   = {summary['initial_a']:.12f} km")
    print(f"e   = {summary['initial_e']:.12f}")
    print(f"eps = {summary['initial_eps']:.12f} km^2/s^2")
    print(f"h   = {summary['initial_h']:.12f} km^2/s")

    print("\nFinal quantities:")
    print(f"a   = {summary['final_a']:.12f} km")
    print(f"e   = {summary['final_e']:.12f}")
    print(f"eps = {summary['final_eps']:.12f} km^2/s^2")
    print(f"h   = {summary['final_h']:.12f} km^2/s")

    print("\nChanges:")
    print(f"delta a   = {summary['delta_a']:.12e} km")
    print(f"delta e   = {summary['delta_e']:.12e}")
    print(f"delta eps = {summary['delta_eps']:.12e} km^2/s^2")
    print(f"delta h   = {summary['delta_h']:.12e} km^2/s")


def print_propagation_summary(result, method_name="Propagation"):
    """
    Print a formatted summary for one propagation result.

    Parameters
    ----------
    result : dict
        Output dictionary from propagate_orbit.
    method_name : str
        Name to display for the propagation method.
    """

    summary = method_summary(result)
    print_method_summary(method_name, summary)


def print_comparison_summary(results):
    """
    Print formatted Euler and RK4 propagation summaries.

    Parameters
    ----------
    results : dict
        Output dictionary from two_body_propagation_euler_rk4.
    """

    summary = build_comparison_summary(results)
    print_method_summary("Euler", summary["euler"])
    print_method_summary("Runge-Kutta 4", summary["rk4"])


def print_orbital_elements(elements):
    """
    Print a formatted summary of orbital elements.

    Parameters
    ----------
    elements : dict
        Dictionary returned by orbital_elements_from_rv.
    """

    print("\nOrbital Elements")
    print("----------------")
    print(f"Orbital type : {elements['orbital_type']}")
    print(f"a            : {elements['a']:.6f} km")
    print(f"e            : {elements['e']:.6f}")
    print(f"i            : {elements['i_deg']:.6f} deg")
    print(f"RAAN         : {elements['RAAN_deg']:.6f} deg")
    print(f"omega / argp : {elements['omega_deg']:.6f} deg")
    print(f"nu           : {elements['nu_deg']:.6f} deg")
    print(f"rp           : {elements['rp']:.6f} km")
    print(f"ra           : {elements['ra']:.6f} km")
    print(f"eps          : {elements['eps']:.6f} km^2/s^2")
    print(f"h            : {elements['h']:.6f} km^2/s")


# Backward-compatible aliases
build_summary = build_comparison_summary
print_summary = print_comparison_summary