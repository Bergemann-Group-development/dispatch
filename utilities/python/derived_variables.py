# -*- coding: utf-8 -*-

"""
   Derived, physical quantities for use in analysis and plotting of DISPATCH2 results.
"""

import numpy as np
import itertools

def divergence(p, varindices, gf, lnorm=True, lcalcnorm=False, lsum=False, lzc=False):
    """
    Return the divergence of a vector field given the components.
    Input must be 3-D numpy arrays (even if one or more of the dimensions
    has length 1).
    
    Originally written by David Clarke.
    Modified by JPR for zone-centred magnetic fields and re-written in Python.

    Note: To apply this function in curvilinear coords., the geometric factors
    must be provided (which is not currently implemented)!
    
    Another note: This function is known to be relatively slow due to the nested do-loops.
    
    v1, v2, v3 =>  components of a vector field in the form of numpy arrays
    lnorm = T =>  normalise divergence
          = F =>  do not normalise divergence
    lcalcnorm = F => do not calculate the normalisation factor
              = T => calculate a normalisation constant for the divergence, but DO NOT
                     apply it, only return its value.
    lsum  = F =>  perform a reduction on the divergence, returning the
                  absolute sum ('sumd')
          = T =>  no reduction
    lzc   = F  =>  input vector components are *face*-centred (e.g., ZEUS, STAGGER)
          = T  =>  input vector components are *zone*-centred (e.g., RAMSES)

    """

    assert isinstance(varindices,list) or isinstance(varindices,tuple), "varindices must be a tuple or a list."
    assert len(varindices) == 3, "you need to supply the indices of three components of a vector field."

    TINY = 1.0e-99

    nx, ny, nz, nv = np.shape(p.data)
    vf1 = p.data[:,:,:,varindices[0]]
    vf2 = p.data[:,:,:,varindices[1]]
    vf3 = p.data[:,:,:,varindices[2]]

    divv = np.empty(vf1.shape)

    if nx > 1: ione = 1
    else: ione = 0
    if ny > 1: jone = 1
    else: jone = 0
    if nz > 1: kone = 1
    else: kone = 0

    # calculate divergence
    if not lzc:
        # face-centred inputs
        dvl1fi = 1.0 / ( gf['dvol1f'] + TINY )
        dvl2fi = 1.0 / ( gf['dvol2f'] + TINY )
        dvl3fi = 1.0 / ( gf['dvol3f'] + TINY )

        for i,j,k in itertools.product(xrange(nx-ione),xrange(ny-jone),xrange(nz-kone)):
            kp1 = k + kone
            jp1 = j + jone
            ip1 = i + ione
            divv[i,j,k] = ( ( gf['dar1c'][ip1] * vf1[ip1,j,k]
                            - gf['dar1c'][i  ] * vf1[i  ,j,k] ) * dvl1fi[i]
                          + ( gf['h32f'][jp1] * vf2[i,jp1,k]
                            - gf['h32f'][j  ] * vf2[i,j  ,k] )
                            * dvl2fi[j] * gf['dar2f'][i]
                          + ( vf3[i,j,kp1] - vf3[i,j,k] )
                            * dvl3fi[k] * gf['dar31f'][i] * gf['dar32f'][j] )
    else:
        # zone-centred inputs
        # note that, to end up with a zone-centred divergence, differences are taken
        # from -1 to +1, and thus skipping i,j,k entirely.
        dvl1ci = 1.0 / ( gf['dvol1c'] + TINY )
        dvl2ci = 1.0 / ( gf['dvol2c'] + TINY )
        dvl3ci = 1.0 / ( gf['dvol3c'] + TINY )

        for i,j,k in itertools.product(xrange(ione,nx-ione),xrange(jone,ny-jone),xrange(kone,nz-kone)):
            km1 = k - kone
            kp1 = k + kone
            jm1 = j - jone
            jp1 = j + jone
            im1 = i - ione
            ip1 = i + ione
            divv[i,j,k] = ( ( gf['dar1f'][ip1] * vf1[ip1,j,k]
                            - gf['dar1f'][im1] * vf1[im1,j,k] )
                            * ( dvl1ci[ip1] + dvl1ci[i] )
                            + ( gf['h32c'][jp1] * vf2[i,jp1,k]
                              - gf['h32c'][jm1] * vf2[i,jm1,k] )
                              * ( dvl2ci[jp1] + dvl2ci[j] ) * gf['dar2f'][i]
                            + ( vf3[i,j,kp1] - vf3[i,j,km1] )
                              * ( dvl3ci[kp1] + dvl3ci[k] )
                              * gf['dar31f'][i] * gf['dar32f'][j] )

    # calculate the reduction, if desired
    sumd = None
    if lsum: sumd = np.sum(np.abs(divv))

    # normalise the divergence, if desired
    nmax = None
    if lnorm or lcalcnorm:
        norm = np.zeros(vf1.shape)

        # Evaluate two normalising constants:
        # 1)  sumn = sum over all grid zones the ratio of (average absolute
        #            vector field) / (sum of grid zone dimensions)
        # 2)  nmax = maximum over the grid of the above ratio.
        if not lzc:
            for i,j,k in itertools.product(xrange(nx-ione),xrange(ny-jone),xrange(nz-kone)):
                kp1 = k + kone
                jp1 = j + jone
                ip1 = i + ione
                norm[i,j,k] = ( ( abs( vf1[ip1,j,k] + vf1[i,j,k] )
                                + abs( vf2[i,jp1,k] + vf2[i,j,k] )
                                + abs( vf3[i,j,kp1] + vf3[i,j,k] ) ) * 0.5
                            / ( ione * gf['dx1f'][i]
                              + jone * gf['h2c'][i] * gf['dx2f'][j]
                              + kone * gf['h31c'][i] * gf['h32c'][j] * gf['dx3f'][k] ) )
        else:
            for i,j,k in itertools.product(xrange(nx-ione),xrange(ny-jone),xrange(nz-kone)):
                norm[i,j,k] = ( ( abs( vf1[i,j,k] )
                                + abs( vf2[i,j,k] )
                                + abs( vf3[i,j,k] ) )
                            / ( ione * gf['dx1f'][i]
                              + jone * gf['h2c'][i] * gf['dx2f'][j]
                              + kone * gf['h31c'][i] * gf['h32c'][j] * gf['dx3f'][k] ) )

        sumn = np.sum(norm)
        nmax = np.nanmax(norm)

        # Apply normalising constants 'sumn' to scalar 'sumd' and 'nmax' to 'divv'.
        if lsum:
            if sumn < TINY:
                sumd = 0.0
            else:
                sumd = sumd / sumn
        if lnorm:
            if nmax < TINY:
                nmaxi = 0.0
            else:
                nmaxi = 1.0 / nmax
            divv = divv * nmaxi

    return divv, sumd, nmax

def gas_temperature (p, dindex, eindex, gf, gamma=5.0/3.0):
    """Calculate the gas temperature given the density and internal energy."""

    #import ipdb
    #ipdb.set_trace()
    e = p.data[:,:,:,eindex]
    d = p.data[:,:,:,dindex]

    return ((gamma-1.0) * e / d)

def speed (p, varidx, gf):
    """Calculate the speed = magnitude of the velocity vector."""

    s1 = p.data[:,:,:,varidx['s1']]
    s2 = p.data[:,:,:,varidx['s2']]
    s3 = p.data[:,:,:,varidx['s3']]
    d  = p.data[:,:,:,varidx['d ']]

    v1 = s1 / d
    v2 = s2 / d / gf['h2c'][:,np.newaxis,np.newaxis]
    v3 = (s3 / d / gf['h31c'][:,np.newaxis,np.newaxis] 
                 / gf['h32c'][np.newaxis,:,np.newaxis])

    # FIXME: technically incorrect for ZEUS and STAGGER data, which should be
    # averaged to zone-centre first.
    return np.sqrt(v1**2 + v2**2 + v3**2)

def gradient_scalar (p, varidx, gf, idir=0):
    """
    Calculate the gradient of a scalar (zone-centred).
    The result is a staggered, 2nd-order difference.

    """

    na = np.newaxis
    n1, n2, n3 = p.gn
    if n1 > 1: ione = 1
    else: ione = 0
    if n2 > 1: jone = 1
    else: jone = 0
    if n3 > 1: kone = 1
    else: kone = 0

    cvar = p.data[:,:,:,varidx]
    if idir == 1:
        g = np.empty_like(cvar)
        for i in xrange(1,n1):
            g[i,:,:] = (cvar[i,:,:] - cvar[i-ione,:,:]) / gf['dx1c'][i]
        return g
    elif idir == 2:
        g = np.empty_like(cvar)
        for j in xrange(1,n2):
            g[:,j,:] = (cvar[:,j,:] - cvar[:,j-jone,:]) / gf['dx2c'][j] / gf['h2c']
        return g
    elif idir == 3:
        g = np.empty_like(cvar)
        for k in xrange(1,n3):
            g[:,:,k] = (cvar[:,:,k] - cvar[:,:,k-kone]) / gf['dx3c'][k] / gf['h31c'] / gf['h32c']
        return g
    else:
        g1 = np.empty_like(cvar)
        for i in xrange(1,n1):
            g1[i,:,:] = (cvar[i,:,:] - cvar[i-ione,:,:]) / gf['dx1c'][i]
        g2 = np.empty_like(cvar)
        for j in xrange(1,n2):
            g2[:,j,:] = (cvar[:,j,:] - cvar[:,j-jone,:]) / gf['dx2c'][j] / gf['h2c']
        g3 = np.empty_like(cvar)
        for k in xrange(1,n3):
            g3[:,:,k] = (cvar[:,:,k] - cvar[:,:,k-kone]) / gf['dx3c'][k] / gf['h31c'] / gf['h32c']
        return (g1,g2,g3)
