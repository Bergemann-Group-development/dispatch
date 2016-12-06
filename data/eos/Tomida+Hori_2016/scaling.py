# =================================================================
# Code units
# =================================================================
import constants as const

l = const.R_Earth                                 # length
m = const.M_Earth                                 # mass
d = m/l**3                                        # density
u = const.kms                                     # velocity
time = l/u                                        # time
p = d*u**2                                        # pressure
e = m*u**2                                        # energy
f = p*u                                           # energy flux = P u
grav = const.grav*(m/l**3)*time**2                # gravity
temp = const.mu * (const.amu/m)/(const.K_B/e)     # temperature for given mu
sigma = const.sigma/m*time**3*temp**4             # stefan - boltzmann constant
#M_p  = const.cgs_M_p/m
#amu  = const.amu/m
#K_B  = const.K_B /e * temp
#h_p  = const.h_p /(e * time)
#c    = const.c / l * time
#AU   = const.cgs_AU / l
#pc   = const.cgs_pc / l
#sec_per_year = const.sec_per_year / time
#Myr  = const.cgs_Myr / time
#M_eart = const.M_Earth / m
#M_Sun  = const.M_Sun / m
#kms    = const.cgs_kms / l * time
#R_Earth = const.cgs_R_Earth / l
#R_gas   = const.cgs_R_gas / e * temp * 
print "=========================================================="
print "CODE UNITS:"
print "length: ", l
print "mass: ", m
print "density: ", d
print "velocity: ", u
print "time: ", time
print "pressure: ", p
print "energy: ", e
print "gravity: ", grav
print "temperature: ", temp,",that is, amu*mu/k_B = 1"
print "=========================================================="
