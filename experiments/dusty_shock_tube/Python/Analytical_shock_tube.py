import numpy as np


#test = 1
def shock_tube( q_l , q_r , gamma , t ):
    #if test == 1:
    # q_l = [ rho_l , u_l , P_l ]
    # q_r = [ rho_r , u_r , P_r ]
    
    # Test-values
    #q_l      = [1.0,0.0,1.0]            # Initial conditions on LHS; [rho_l,u_l,P_l]
    #q_r      = [0.125,0.0,0.1]          # Initial conditions on RHS; [rho_r,u_r,P_r]
    #gamma = 1.4
    #t     = 0.2
    
    # Allocate initial parameters
    rho_l = q_l[0]; rho_r = q_r[0]
    u_l   = q_l[1]; u_r   = q_r[1]
    P_l   = q_l[2]; P_r   = q_r[2]
    
    c_l   = np.sqrt(gamma*P_l/rho_l)
    c_r   = np.sqrt(gamma*P_r/rho_r)
    
    P     = np.linspace(P_r,P_l,num=500)
    
    x     = np.linspace(0.0,1.0,num=500)
    x_0   = 0.5
    
    if P_l > P_r: #For a rarefaction wave to the left and shock wave to the right
        # Compute intermediate values of rho, u, P around contact discontinuity
        phi_l  = u_l + 2.0*c_l/(gamma-1.0) * ( 1.0 - (P/P_l)**((gamma-1.0)/(2.0*gamma)) )
        
        beta   = (gamma+1.0)/(gamma-1.0)
        phi_r  = u_r - 2.0*c_r/np.sqrt(2.0*gamma*(gamma-1.0)) * ( (1.0-P/P_r)/np.sqrt(1.0+beta*P/P_r) )
        
        idx_s  = np.argmin(abs(phi_l-phi_r))
        P_s    = P[idx_s]
        u_s    = phi_l[idx_s]
        
        rho_sl = (P_s/P_l)**(1.0/gamma) * rho_l
        rho_sr = ((1.0+beta*P_s/P_r)/(P_s/P_r+beta)) * rho_r
        
        c_sl   = np.sqrt(gamma*P_s/rho_sl)
        c_sr   = np.sqrt(gamma*P_s/rho_sr)
        
        # Compute values of rho, u, P in rarefaction wave
        u      = ( (gamma-1.0)*u_l+2.0*(c_l+(x-x_0)/t) )/( gamma+1.0 )
        rho    = ( ( rho_l**gamma * (u-(x-x_0)/t)**2.0 )/( gamma*P_l ) )**(1.0/(gamma-1.0))
        P      = ( (P_l)/(rho_l**gamma) ) * (rho)**gamma
        
        # Compute all x-values and corresponding indices for all wave transitions
        x_1    = x_0 + (u_l-c_l)*t                               # Left edge of rarefaction wave
        idx_1  = np.argmin(abs(x-x_1))
        x_2    = x_0 + (u_s-c_sl)*t                              # Right edge of rarefaction wave
        idx_2  = np.argmin(abs(x-x_2))
        x_3    = x_0 + u_s*t                                     # Location of contact discontinuity
        idx_3  = np.argmin(abs(x-x_3))
        x_4    = x_0 + (rho_sr*u_s-rho_r*u_r)/(rho_sr-rho_r) * t # Location of shock wave
        idx_4  =  np.argmin(abs(x-x_4))
        
        # Insert correct values of rho, u, P and in the different regions
        rho[0:idx_1] = rho_l;   rho[idx_2:idx_3] = rho_sl;   rho[idx_3:idx_4] = rho_sr;   rho[idx_4:] = rho_r
        u[0:idx_1]   = u_l;     u[idx_2:idx_4]   = u_s;      u[idx_4:]        = u_r
        P[0:idx_1]   = P_l;     P[idx_2:idx_4]   = P_s;      P[idx_4:]        = P_r
        s = ( ( np.log(P)-np.log(rho)*gamma )*rho )/( gamma-1.0 )
        p = u*rho
    else: #For a shock wave to the left and rarefaction wave to the right
        # Compute intermediate values of rho, u, P around contact discontinuity
        phi_r  = u_r - 2.0*c_r/(gamma-1.0) * ( 1.0 - (P/P_r)**((gamma-1.0)/(2.0*gamma)) )
        
        beta   = (gamma+1.0)/(gamma-1.0)
        phi_l  = u_l + 2.0*c_l/np.sqrt(2.0*gamma*(gamma-1.0)) * ( (1.0-P/P_l)/np.sqrt(1.0+beta*P/P_l) )
        
        idx_s  = np.argmin(abs(phi_l-phi_r))
        P_s    = P[idx_s]
        u_s    = phi_l[idx_s]
        
        rho_sr = (P_s/P_r)**(1.0/gamma) * rho_r
        rho_sl = ((1.0+beta*P_s/P_l)/(P_s/P_l+beta)) * rho_l
        
        c_sl   = np.sqrt(gamma*P_s/rho_sl)
        c_sr   = np.sqrt(gamma*P_s/rho_sr)
        
        # Compute values of rho, u, P in rarefaction wave
        u      = ( (gamma-1.0)*u_r-2.0*(c_r-x/t) )/( gamma+1.0 )
        rho    = ( ( rho_r**gamma * (x/t-u)**2.0 )/( gamma*P_r ) )**(1.0/(gamma-1.0))
        P      = ( (P_r)/(rho_r**gamma) ) * (rho)**gamma
        
        # Compute all x-values and corresponding indices for all wave transitions
        x_1    = x_0 + (rho_l*u_l-rho_sl*u_s)/(rho_l-rho_sl) * t # Location of shock wave
        idx_1  = np.argmin(abs(x-x_1))
        x_2    = x_0 + u_s*t                                     # Location of contact discontinuity
        idx_2  = np.argmin(abs(x-x_2))
        x_3    = x_0 + (u_s+c_sr)*t                              # Left edge of rarefaction wave
        idx_3  = np.argmin(abs(x-x_3))
        x_4    = x_0 + (u_r+c_r)*t                               # Right edge of rarefaction wave
        idx_4  =  np.argmin(abs(x-x_4))
        
        # Insert correct values of rho, u, P and in the different regions
        rho[0:idx_1] = rho_l;   rho[idx_1:idx_2] = rho_sl;   rho[idx_2:idx_3] = rho_sr;   rho[idx_4:] = rho_r
        u[0:idx_1]   = u_l;     u[idx_1:idx_3]   = u_s;      u[idx_4:]        = u_r
        P[0:idx_1]   = P_l;     P[idx_1:idx_3]   = P_s;      P[idx_4:]        = P_r
        s = ( ( np.log(P)-np.log(rho)*gamma )*rho )/( gamma-1.0 )
        p = u*rho
    q = np.array([rho,s,p,P,u,x])
    return q