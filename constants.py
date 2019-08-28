import math

#: Gravitational constant
G_M3_KGS2 = 6.67428e-11
#: Planck's constant
HPLANCK_M2KG_S = 6.62607004e-34
#: Speed of light
CLIGHT_M_S = 2.998e8

#: Solar mass
MSUN_KG = 1.988416e30
#: Solar radius
RSUN_M = 6.957e8
#: Astronomical Unit
AU_M = 1.49598e11
#: Parsec
PARSEC_M = 3.086e16
#: Earth mass
MEARTH_KG = 5.9722e24
#: Earth radius
REARTH_M = 6.3781e6
#: Jovial mass
MJUPITER_KG = 1.89813e27
# Jovial radius
RJUPITER_M = 7.1492e7

#: Seconds in 1 day
DAY_S = 86400.
#: Seconds in 1 minute
MINUTE_S = 60
#: Seconds in 1 hour
HOUR_S = 1440

#: 1 AU in Earth radii
AU_REARTH = AU_M / REARTH_M
#: 1 solar mass in Earth masses
MSUN_MEARTH = MSUN_KG / MEARTH_KG
#: 1 solar radius in Earth radii
RSUN_REARTH = RSUN_M / REARTH_M
#: Gravitational constant in Earth units (Rearth^3 * Mearth^-1 * day^-2)
G_EARTH = 11468

#: 1 radian in deg
RAD_DEG = math.pi / 180
