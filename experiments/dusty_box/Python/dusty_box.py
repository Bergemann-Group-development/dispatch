import numpy as np

def dusty_box( regime , t_end , vx_gas , vx_dust , rho_gas , rho_dust , K , a ):

    dv0 = vx_gas - vx_dust
    rhoi = rho_gas**(-1) + rho_dust**(-1)
    if t_end > 0.0:
        # calculate results for a range of times
        t = np.linspace(0.0,t_end,200)
    else:
        # by choosing a negative `t_end`, one can calculate the solution at
        # a single time given by `abs(t_end)`.
        t = abs(t_end)

    vs  = (rho_gas*vx_gas+rho_dust*vx_dust)/(rho_gas+rho_dust)

    if regime == 'linear':
        dv  = dv0 * np.exp( -K*rhoi*t )
    elif regime == 'quadratic':
        den = 1.0 + abs(dv0)*K*rhoi*t
        dv  = dv0/den
    elif regime == 'power':
        den = 1.0 + a*(abs(dv0)**a)*K*rhoi*t
        dv  = dv0/(den**(1.0/a))
    elif regime == 'third_order':
        e   = np.exp( -K*rhoi*t )
        den = 1.0+a*(dv0**2)*( 1.0-(e**2) )
        dv  = dv0*e/np.sqrt(den)
    elif regime == 'mixed':
        e   = np.sign(dv0)
        s   = np.sinh(K*rhoi*t)
        c   = np.cosh(K*rhoi*t)
        num = s+np.sqrt(1.0+a*(dv0**2))*c
        den = c+np.sqrt(1.0+a*(dv0**2))*s
        dv  = (e/np.sqrt(a))*np.sqrt( ( (num/den)**2 ) - 1.0 )

    dvg = vs + rho_dust*dv/(rho_gas+rho_dust)
    dvd = vs - rho_gas *dv/(rho_gas+rho_dust)
    q   = np.array([t,dv,dvg,dvd])
    return q
