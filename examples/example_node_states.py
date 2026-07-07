import sys
from pathlib import Path

import numpy as np

project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from constants import MU_BODIES
from orbital_elements import coe_to_rv, orbital_elements_from_rv


def main():
    mu = MU_BODIES["Earth"]

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

    r_check_asc, v_check_asc = coe_to_rv(a, e, i, raan, argp, nu_asc, mu)
    r_check_des, v_check_des = coe_to_rv(a, e, i, raan, argp, nu_des, mu)

    r_error_asc = np.linalg.norm(r_vec_asc - r_check_asc)
    v_error_asc = np.linalg.norm(v_vec_asc - v_check_asc)

    r_error_des = np.linalg.norm(r_vec_des - r_check_des)
    v_error_des = np.linalg.norm(v_vec_des - v_check_des)

    print("Node State Example")
    print("------------------")
    print(f"Ascending node true anomaly [deg]  : {np.degrees(nu_asc):.8f}")
    print(f"Ascending node r error [km]        : {r_error_asc:.12e}")
    print(f"Ascending node v error [km/s]      : {v_error_asc:.12e}")
    print()
    print(f"Descending node true anomaly [deg] : {np.degrees(nu_des):.8f}")
    print(f"Descending node r error [km]       : {r_error_des:.12e}")
    print(f"Descending node v error [km/s]     : {v_error_des:.12e}")


if __name__ == "__main__":
    main()
