
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from constants import MU_BODIES
from transfers import hohmann_transfer_requirements


def main():
    mu = MU_BODIES["Earth"]

    r1 = 7000.0
    r2 = 8000.0

    reqs = hohmann_transfer_requirements(r1, r2, mu)

    print("Hohmann Transfer Example")
    print("------------------------")
    print(f"Transfer type        : {reqs['transfer_type']}")
    print(f"Initial radius [km]  : {reqs['r1']:.2f}")
    print(f"Final radius [km]    : {reqs['r2']:.2f}")
    print(f"Transfer a [km]      : {reqs['a_transfer']:.2f}")
    print()
    print(f"v circular 1 [km/s]  : {reqs['v_circular_1']:.8f}")
    print(f"v transfer 1 [km/s]  : {reqs['v_transfer_at_r1']:.8f}")
    print(f"dv1 [km/s]           : {reqs['dv1']:.8f}")
    print()
    print(f"v transfer 2 [km/s]  : {reqs['v_transfer_at_r2']:.8f}")
    print(f"v circular 2 [km/s]  : {reqs['v_circular_2']:.8f}")
    print(f"dv2 [km/s]           : {reqs['dv2']:.8f}")
    print()
    print(f"total delta-v [km/s] : {reqs['total_dv']:.8f}")
    print(f"transfer time [s]    : {reqs['transfer_time']:.8f}")


if __name__ == "__main__":
    main()