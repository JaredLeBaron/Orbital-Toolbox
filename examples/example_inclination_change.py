
import sys
from pathlib import Path

import numpy as np

project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from constants import MU_BODIES
from maneuvers import inclination_change
from orbital_elements import coe_to_rv, orbital_elements_from_rv


def main():
    mu = MU_BODIES["Earth"]

    a = 7000.0
    e = 0.0
    i_initial = 0.0
    raan = 0.0
    argp = 0.0
    nu = 0.0

    i_final = np.radians(45.0)

    r_vec, v_vec = coe_to_rv(a, e, i_initial, raan, argp, nu, mu)

    delta_v_vec, v_vec_new, delta_v_mag, scalar_check = inclination_change(
        r_vec,
        v_vec,
        i_initial,
        i_final,
    )

    elements_after = orbital_elements_from_rv(r_vec, v_vec_new, mu)

    print("Inclination Change Example")
    print("--------------------------")
    print(f"Initial inclination [deg]  : {np.degrees(i_initial):.6f}")
    print(f"Target inclination [deg]   : {np.degrees(i_final):.6f}")
    print(f"Recovered inclination [deg]: {elements_after['i_deg']:.6f}")
    print()
    print(f"Initial speed [km/s]       : {np.linalg.norm(v_vec):.8f}")
    print(f"Final speed [km/s]         : {np.linalg.norm(v_vec_new):.8f}")
    print()
    print(f"delta-v vector [km/s]      : {delta_v_vec}")
    print(f"delta-v magnitude [km/s]   : {delta_v_mag:.8f}")
    print(f"scalar check [km/s]        : {scalar_check:.8f}")


if __name__ == "__main__":
    main()