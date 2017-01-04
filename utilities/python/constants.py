#
"""
    Useful physical constants with units given in the comments.
"""

#  radiation and atomic physics:
clight = 2.99792e10        # speed of light [cm s**-1]
SBc  = 5.670373e-5         # Stefan-Boltzmann constant [erg s**-1 cm**-2 K**-4]
arad = 4.0 * SBc / clight  # radiation constant [erg cm**-3 K**-4]
kBoltz = 1.3807e-16        # Boltzmann's constant [erg/K]
mpro = 1.6726e-24          # proton mass [g]
mamu = 1.66054e-24         # atomic mass unit [g]
hPlanck = 6.6260755e-27    # Planck's constant [erg s]

#  astronomical:
au    = 1.49598e13         # astronomical unit [cm]
m_sun  = 1.98892e33        # mass of the Sun [g]
r_sun  = 6.9598e10         # radius of the Sun [cm]
l_sun  = 3.839e33          # luminosity of the Sun [erg s**-1]
gcnst = 6.67e-8            # gravitational constant [cm**3 g**-1 s**-2]
pc    = 3.0856776e18       # parsec [cm]
m_earth = 5.972e27         # mass of the Earth [g]
r_earth = 6.371e8          # radius of the Earth [cm]

#  conversion factors:
eV_to_K = 1.1604505e9      # electron volts to Kelvin
Habing = 1.6e-3            # CGS units of flux [erg cm**-3 s**-1] to Habing units
yr_to_s = 3.156e7          # years to seconds
