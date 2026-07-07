
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

    print("COE to RV Example")
    print("-----------------")
    print(f"r_vec [km]   : {r_vec}")
    print(f"v_vec [km/s] : {v_vec}")
    print()
    print("Recovered elements")
    print(f"a [km]       : {elements['a']:.6f}")
    print(f"e            : {elements['e']:.12f}")
    print(f"i [deg]      : {elements['i_deg']:.6f}")
    print(f"RAAN [deg]   : {elements['RAAN_deg']:.6f}")
    print(f"argp [deg]   : {elements['argp_deg']:.6f}")
    print(f"nu [deg]     : {elements['nu_deg']:.6f}")


if __name__ == "__main__":
    main()