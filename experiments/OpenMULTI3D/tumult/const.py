"""Physical constants (CGS) - the single source for the whole package.

These are FROZEN legacy values, deliberately not astropy/CODATA-current: the atomic
physics that consumes them (partition functions and Saha in atom.py, Einstein
coefficients in trans.py, level energies in level.py) was validated against the Fortran
side (MULTI/dispatch) with exactly these numbers. Updating them to modern precision is a
physics change, not a refactor - if that is ever wanted, regenerate the golden files in
tests/ deliberately and revalidate against the Fortran output.

Unit-aware code should keep using astropy.constants directly; this module exists for the
raw-CGS numerics that must stay bit-compatible with the legacy chain.
"""
import numpy as np

EE = 1.602189E-12    # electron charge equivalent energy of 1 eV [erg]   (astropy: e/erg)
HH = 6.626176E-27    # Planck constant [erg s]                           (astropy: h.cgs)
CC = 2.99792458E10   # speed of light [cm/s]                             (exact, matches c.cgs)
EM = 9.109534E-28    # electron mass [g]                                 (astropy: m_e.cgs)
UU = 1.6605655E-24   # atomic mass unit [g]                              (astropy: u.cgs)
BK = 1.380662E-16    # Boltzmann constant [erg/K]                        (astropy: k_B.cgs)
PI = 3.14159265359   # historic truncation of pi - keep for bit-compat (np.pi differs!)

EC = 4.80325E-10     # electron charge [esu]
MH = 1.6726e-24      # proton mass [g]
XMH = 1.67339E-24    # mass of the HYDROGEN ATOM [g] (proton + electron; != MH)
SBOL = 5.6704e-5     # Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]
RYDV = 13.605        # Rydberg energy [eV]

CC_KMS = 299792.458  # speed of light [km/s]

# derived combinations (wavelengths in AA where the 1e8/1e24 factors appear)
HCE = HH * CC / EE * 1e8   # hc/e with lambda in AA -> eV
HC2 = 2. * HH * CC * 1e24  # 2hc with lambda in AA (Planck / Einstein-B prefactor)
HCK = HH * CC / BK * 1e8   # hc/k with lambda in AA -> K
EK = EE / BK               # eV -> K

# Saha-equation prefactor. NOTE: written exactly as atom.py historically computed it
# (.../np.sqrt(2*PI*EM)/np.sqrt(BK), NOT .../np.sqrt(2*PI*EM*BK)) - the two differ in the
# last ULP and the partition-function golden (tests/golden/atom.npz) pins this one.
saha_const = 0.5 * (HH / np.sqrt(2. * PI * EM) / np.sqrt(BK)) ** 3
