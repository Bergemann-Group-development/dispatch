from numpy import abs, zeros, min, product
from time import time

def slopes (xx, slope_type, timeit=False):
     """
     Computes slopes for hydro and radiative transfer
     variables in 3 dimensions.
     The different slope limiters are:
     - 1 : Minmod
     - 2 : Moncen
     - 3 : Positivity preserving
     - 4 : Harmonic average
     - 5 : Van Albada & Van Leer
     - 6 : Superbee
     """
     m = xx.shape
     ndim = 3-m.count(1)
     m = list(m)
     m.append(ndim)
     pxxG = zeros(m)
     pxxD = zeros(m)

     if timeit:
         start=time()

     if(slope_type == 0):  
        return zeros(m)

     elif slope_type==1: # MINMOD
        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):
	      
                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 pxxG[i,j,k,0] = MinMod(dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    pxxG[i,j,k,1] = MinMod(dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    pxxG[i,j,k,2] = MinMod(dfmmr, -dfmml)

                 pxxD[i,j,k,:] = pxxG[i,j,k,:]

     elif slope_type==2: # MONCEN

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = 2.0 * ( xx[i,  j,  k  ]-xx[i-1,j  ,k  ] )
                 dfrmm = 2.0 * ( xx[i+1,j  ,k  ]-xx[i  ,j  ,k  ] )
                 dfcen = ( xx[i+1,j  ,k  ]-xx[i-1,j  ,k  ] ) / 2.0
                 pxxG[i,j,k,0] = MonCen(dflmm,dfrmm,dfcen)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = 2.0 * ( xx[i,  j,  k  ]-xx[i  ,j-1,k  ] )
                    dfmrm = 2.0 * ( xx[i  ,j+1,k  ]-xx[i  ,j  ,k  ] )
                    dfcen = ( xx[i  ,j+1,k  ]-xx[i  ,j-1,k  ] ) / 2.0
                    pxxG[i,j,k,1] = MonCen(dfmlm,dfmrm,dfcen)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = 2.0 * ( xx[i,  j,  k  ]-xx[i  ,j  ,k-1] )
                    dfmmr = 2.0 * ( xx[i  ,j  ,k+1]-xx[i  ,j  ,k  ] )
                    dfcen = ( xx[i  ,j  ,k+1]-xx[i  ,j  ,k-1] ) / 2.0
                    pxxG[i,j,k,2] = MonCen(dfmml,dfmmr,dfcen)

                 pxxD[i,j,k,:] = pxxG[i,j,k,:]

     elif slope_type==2: # positivity preserving

        if( ndim == 1):
            for i in range(1,m[0]-1):
               for j in range(1,m[1]-1):
                  for k in range(1,m[2]-1):
                    dflm = xx[i-1]-xx[i]
                    dfmm = xx[i  ]-xx[i]
                    dfrm = xx[i+1]-xx[i]

                    vmin = min(dflm,dfmm,dfrm)
                    vmax = max(dflm,dfmm,dfrm)

                    dfx  = 0.5*(xx[i+1]-xx[i-1])
                    dff  = abs(dfx)

                    if (dff > 0.0):
                       slop = min(1.0,min(abs(vmin),abs(vmax))/dff)
                    else:
                       slop = 1.0

                    dlim = slop

                    pxxG[i,0] = dlim*dfx
                    pxxD[i,:] = pxxG[i,:]

        elif (ndim == 2):

            for i in range(1,m[0]-1):
               for j in range(1,m[1]-1):
                  for k in range(1,m[2]-1):
                    dfll = xx[i-1,j-1]-xx[i,j]
                    dflm = xx[i-1,j  ]-xx[i,j]
                    dflr = xx[i-1,j+1]-xx[i,j]
                    dfml = xx[i  ,j-1]-xx[i,j]
                    dfmm = xx[i  ,j  ]-xx[i,j]
                    dfmr = xx[i  ,j+1]-xx[i,j]
                    dfrl = xx[i+1,j-1]-xx[i,j]
                    dfrm = xx[i+1,j  ]-xx[i,j]
                    dfrr = xx[i+1,j+1]-xx[i,j]

                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

                    dfx  = 0.5*(xx[i+1,j]-xx[i-1,j])
                    dfy  = 0.5*(xx[i,j+1]-xx[i,j-1])
                    dff  = 0.5*(abs(dfx)+abs(dfy))

                    if(dff>0.0):
                       slop = min(1.0,min(abs(vmin),abs(vmax))/dff)
                    else:
                       slop = 1.0

                    dlim = slop

                    pxxG[i,j,0] = dlim*dfx
                    pxxG[i,j,1] = dlim*dfy
                    pxxD[i,j,:] = pxxG[i,j,:]
           
        else:
            for i in range(1,m[0]-1):
               for j in range(1,m[1]-1):
                  for k in range(1,m[2]-1):

                    dflll = xx[i-1,j-1,k-1]-xx[i,j,k]
                    dflml = xx[i-1,j  ,k-1]-xx[i,j,k]
                    dflrl = xx[i-1,j+1,k-1]-xx[i,j,k]
                    dfmll = xx[i  ,j-1,k-1]-xx[i,j,k]
                    dfmml = xx[i  ,j  ,k-1]-xx[i,j,k]
                    dfmrl = xx[i  ,j+1,k-1]-xx[i,j,k]
                    dfrll = xx[i+1,j-1,k-1]-xx[i,j,k]
                    dfrml = xx[i+1,j  ,k-1]-xx[i,j,k]
                    dfrrl = xx[i+1,j+1,k-1]-xx[i,j,k]

                    dfllm = xx[i-1,j-1,k  ]-xx[i,j,k]
                    dflmm = xx[i-1,j  ,k  ]-xx[i,j,k]
                    dflrm = xx[i-1,j+1,k  ]-xx[i,j,k]
                    dfmlm = xx[i  ,j-1,k  ]-xx[i,j,k]
                    dfmmm = 0.0
                    dfmrm = xx[i  ,j+1,k  ]-xx[i,j,k]
                    dfrlm = xx[i+1,j-1,k  ]-xx[i,j,k]
                    dfrmm = xx[i+1,j  ,k  ]-xx[i,j,k]
                    dfrrm = xx[i+1,j+1,k  ]-xx[i,j,k]

                    dfllr = xx[i-1,j-1,k+1]-xx[i,j,k]
                    dflmr = xx[i-1,j  ,k+1]-xx[i,j,k]
                    dflrr = xx[i-1,j+1,k+1]-xx[i,j,k]
                    dfmlr = xx[i  ,j-1,k+1]-xx[i,j,k]
                    dfmmr = xx[i  ,j  ,k+1]-xx[i,j,k]
                    dfmrr = xx[i  ,j+1,k+1]-xx[i,j,k]
                    dfrlr = xx[i+1,j-1,k+1]-xx[i,j,k]
                    dfrmr = xx[i+1,j  ,k+1]-xx[i,j,k]
                    dfrrr = xx[i+1,j+1,k+1]-xx[i,j,k]

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, 
                               dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, 
                               dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, 
                               dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, 
                               dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = 0.5*(xx[i+1,j,k]-xx[i-1,j,k])
                    dfy  = 0.5*(xx[i,j+1,k]-xx[i,j-1,k])
                    dfz  = 0.5*(xx[i,j,k+1]-xx[i,j,k-1])
                    dff  = 0.5*(abs(dfx)+abs(dfy)+abs(dfz))

                    if (dff > 0.0):
                       dlim = min(1.0,min(abs(vmin),abs(vmax))/dff)
                    else:
                       dlim = 1.0

                    pxxG[i,j,k,0] = dlim*dfx

                    if (ndim > 1):  # slopes along y direction
                       pxxG[i,j,k,1] = dlim*dfy

                    if (ndim > 2):  # slopes along z direction
                       pxxG[i,j,k,2] = dlim*dfz

                    pxxD[i,j,k,:] = pxxG[i,j,k,:]

     elif slope_type==4: # HARMONIC AVERAGE

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 pxxG[i,j,k,0] = MoyHarm(dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    pxxG[i,j,k,1] = MoyHarm(dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    pxxG[i,j,k,2] = MoyHarm(dfmmr, -dfmml)

                 pxxD[i,j,k,:] = pxxG[i,j,k,:]

     elif slope_type==5: # Van Albada Van Leer

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 pxxG[i,j,k,0] = VAVL(dfrmm, -dflmm )

                 if(ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    pxxG[i,j,k,1] = VAVL(dfmrm, -dfmlm)

                 if(ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    pxxG[i,j,k,2] = VAVL(dfmmr, -dfmml)

                 pxxD[i,j,k,:] = pxxG[i,j,k,:]

     elif slope_type==6: # Superbee

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 pxxD[i,j,k,0] = MinMod(dfrmm, -2.0*dflmm)
                 pxxG[i,j,k,0] = MinMod(2.0*dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    pxxD[i,j,k,1] = MinMod(dfmrm, -2.0*dfmlm)
                    pxxG[i,j,k,1] = MinMod(2.0*dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    pxxD[i,j,k,2] = MinMod(dfmmr, -2.0*dfmml)
                    pxxG[i,j,k,2] = MinMod(2.0*dfmmr, -dfmml)

     if timeit:
         musec = (time()-start)/product(m[0:3])*1e6
         print('{:.2f} microsec per cell'.format(musec))

     return pxxD, pxxG

def SIGN(a,b):
    if b>0:
        return a
    else:
        return -a

def MinMod(pa,pb):
    """ Minmod averaging function """
    MMod = 0.5*min([abs(pa),abs(pb)])*(SIGN(1.0,pa) + SIGN(1.0,pb))
    return MMod

def MonCen (pa,pb,pc):
    """ Monotonised centred averaging function """
    if ((pa*pb < 0.0) or (pa*pc < 0.0) or (pb*pc < 0.0)):
        mcen = 0.0    
    elif ((abs(pa) < abs(pb)) and (abs(pa) < abs(pc))):
        mcen = pa
    elif ((abs(pb) < abs(pa)) and (abs(pb) < abs(pc))):
        mcen = pb
    else:
        mcen = pc
    return mcen

def MoyHarm(pa,pb):
    """ Harmonic averaging function """
    if (pa*pb > 0.0):
       MMod = 2.*pa*pb/(pa+pb)
    else:
       MMod = 0.
    return MMod

def VAVL(pa,pb):
    """ Van Albada & Van Leer averaging function """
    if(pa*pb > 0.0):
       e = 1.0e-8    
       num = ( pa*pa + e*e )*pb + ( pb*pb + e*e )*pa
       den = ( pa*pa +  pb*pb + 2.0*e*e )
       vv  = num/den
    else:
       vv  = 0.0    
    return vv
