
import sys
from pathlib import Path

import numpy as np


project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from constants import MU_BODIES
from orbital_elements import coe_to_rv, orbital_elements_from_rv


def main():
    mu = MU_BODIES["Earth"]


    a = 12000
    e = 0.3
    i = np.radians(40)
    raan = np.radians(30)
    argp = np.radians(60)
    nu = np.radians(45)


    r_vec, v_vec = coe_to_rv(a, e, i, raan, argp, nu, mu)


    elements = orbital_elements_from_rv(r_vec, v_vec, mu)

    nu_asc = elements["nu_asc_rad"]
    nu_des = elements["nu_des_rad"]

    r_vec_asc = elements["r_vec_asc"]
    v_vec_asc = elements["v_vec_asc"]

    r_vec_des = elements["r_vec_des"]
    v_vec_des = elements["v_vec_des"]


    r_check_asc, v_check_asc = coe_to_rv(a, e, i, raan, argp, nu_asc, mu)
    r_check_des, v_check_des = coe_to_rv(a, e, i, raan, argp, nu_des, mu)


    r_error_asc = np.linalg.norm(r_vec_asc - r_check_asc)
    v_error_asc = np.linalg.norm(v_vec_asc - v_check_asc)

    r_error_des = np.linalg.norm(r_vec_des - r_check_des)
    v_error_des = np.linalg.norm(v_vec_des - v_check_des)

    print("Ascending Node")
    print("----------------")
    print(f"nu_asc [deg]: {np.degrees(nu_asc):.8f}")
    print(f"r error [km]: {r_error_asc:.12e}")
    print(f"v error [km/s]: {v_error_asc:.12e}")
    print()

    print("Descending Node")
    print("----------------")
    print(f"nu_des [deg]: {np.degrees(nu_des):.8f}")
    print(f"r error [km]: {r_error_des:.12e}")
    print(f"v error [km/s]: {v_error_des:.12e}")
    print()

    tolerance_r = 1e-6
    tolerance_v = 1e-9

    asc_pass = r_error_asc < tolerance_r and v_error_asc < tolerance_v
    des_pass = r_error_des < tolerance_r and v_error_des < tolerance_v

    if asc_pass and des_pass:
        print("PASS: Node state calculation matches coe_to_rv.")
    else:
        print("FAIL: Node state calculation does not match coe_to_rv.")


if __name__ == "__main__":
    main()