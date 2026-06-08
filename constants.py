
# Jared LeBaron
# Constants

"""
Physical constants used by the orbital toolbox.

This module stores commonly used physical constants such as gravitational
parameters and planetary radii.

Units:
- Gravitational parameter (mu): km^3/s^2
- Radius: km
"""

# Standard gravitational parameters [km^3/s^2]
MU_BODIES = {
    "Sun":     1.3271244004127942e11,
    "Mercury": 22031.868551,
    "Venus":   324858.592000,
    "Earth":   398600.435507,
    "Moon":    4902.800118,
    "Mars":    42828.375816,
    "Jupiter": 126712764.100000,
    "Saturn":  37940584.841800,
    "Uranus":  5794556.400000,
    "Neptune": 6836527.100580,
    "Pluto":   975.500000,
}


# Mean body radii [km]
RADIUS_BODIES = {
    "Sun":     695700.0,
    "Mercury": 2439.7,
    "Venus":   6051.8,
    "Earth":   6378.137,
    "Moon":    1737.4,
    "Mars":    3389.5,
    "Jupiter": 69911.0,
    "Saturn":  58232.0,
    "Uranus":  25362.0,
    "Neptune": 24622.0,
    "Pluto":   1188.3,
}