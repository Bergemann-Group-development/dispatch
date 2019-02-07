# -*- coding: utf-8 -*-

"""
   Coordinate transformation functions for DISPATCH2 data.
"""

import numpy as np


def rtp_to_xyz(r, theta, phi, theta_rotate=0.0, phi_rotate=0.0):
    """
    Convert 1-D spherical (r, theta, phi) coords. to 3-D Cartesian coords.
    (x, y, z) for plotting purposes.

    theta_rotate: rotate the coordinate system in the theta-direction by
                  `theta_rotate`; the default gives a r-coord that overlies
                  the -z-coord when theta=0.
    phi_rotate: rotate the coordinate system in the phi-direction by
                `phi_rotate`; the default results in a r-coord that overlies
                the x-coord when theta = pi/2.

    """

    # First establish 3-D versions of r, theta, and phi.
    rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')

    # Calculate 3-D Cartesian coordinate vectors.
    xx =  rr * np.sin(tt + theta_rotate) * np.cos(pp + phi_rotate)
    yy = -rr * np.sin(tt + theta_rotate) * np.sin(pp + phi_rotate)
    zz = -rr * np.cos(tt + theta_rotate)

    return xx, yy, zz

def zrp_to_xyz(z, r, phi, phi_rotate=0.0):
    """
    Convert 1-D cylindrical (Z, R, Phi) coords. to 3-D Cartesian coords.
    (x, y, z) for plotting purposes.

    phi_rotate: rotate the coordinate system in the phi-direction by
                'phi_rotate'; the the default gives a R-coord. that overlies
                the x-coord when phi = 0.

    """

    # First establish 3-D versions of Z, R, and Phi.
    ZZ, RR, PP = np.meshgrid(z, r, phi, indexing='ij')

    # Calculate 3-D Cartesian coordinate vectors.
    xx = RR * np.cos(PP + phi_rotate)
    yy = RR * np.sin(PP + phi_rotate)
    zz = ZZ

    return xx, yy, zz

def rtp_to_zrp(r, theta, phi, theta_rotate=0.0, phi_rotate=0.0):
    """
    Convert 1-D spherical (r, theta, phi) coords. to 3-D cylindrical coords.
    (Z, R, Phi).

    theta_rotate: rotate the coordinate system in the theta-direction by
                  `theta_rotate`; the default gives a r-coord that overlies
                  the -z-coord when theta=0.
    phi_rotate: rotate the coordinate system in the phi-direction by
                `phi_rotate`; the default results in a r-coord that overlies
                the R-coord when theta = pi/2.

    """

    # First establish 3-D versions of r, theta, and phi.
    rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')

    # Calculate 3-D cylindrical coordinate vectors.
    ZZ = -rr * np.cos(tt + theta_rotate)
    RR =  rr * abs(np.sin(tt + theta_rotate))
    PP = -pp - phi_rotate

    return ZZ, RR, PP

def zrp_to_rtp(z, r, phi, phi_rotate=0.0):
    """
    Convert 1-D cylindrical (Z, R, Phi) coords. to 3-D spherical coords.
    (r, theta, phi).

    phi_rotate: rotate the coordinate system in the phi-direction by
                'phi_rotate'; the default gives a R-coord which overlies
                the r-coord when theta=pi/2.

    """

    # First establish 3-D versions of Z, R, and Phi.
    ZZ, RR, PP = np.meshgrid(z, r, phi, indexing='ij')

    # Calculate 3-D spherical coordinate vectors.
    rr = np.sqrt(RR**2 + ZZ**2)
    tt = np.arctan2(RR, ZZ)
    pp = -PP - phi_rotate

    return rr, tt, pp

def xyz_to_zrp(x, y, z):
    """
    Convert 1-D Cartesian (x, y, z) coords. to 3-D cylindrical coords.
    (Z, R, Phi).

    The z- and Z-coords. are assumed equivalent.

    """

    # First establish 3-D versions of x, y, z.
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

    # Calculate 3-D cylindrical coordinate vectors.
    ZZ = zz
    RR = np.sqrt(xx**2 + yy**2)
    PP = np.arctan2(yy, xx)

    return ZZ, RR, PP

def xyz_to_rtp(x, y, z):
    """
    Convert 1-D Cartesian (x, y, z) coords. to 3-D spherical coords.
    (r, theta, phi).

    The z-coord. is assumed to be anti-parallel to the r-coord. when
    theta = 0.

    """

    # First establish 3-D versions of x, y, z
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

    # Calculate 3-D spherical coordinate vectors.
    rr = np.sqrt(xx**2 + yy**2 + zz**2)
    tt = np.arccos(zz / rr)
    pp = np.arccos(xx / np.sqrt(xx**2 + yy**2))

    return rr, tt, pp
