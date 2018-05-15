import numpy as np
#import matplotlib.pyplot as pl


#test = 1
def shock_tube( q_l , q_r , gamma , t ):
    #if test == 1:
    # q_l = [ rho_l , u_l , P_l ]
    # q_r = [ rho_r , u_r , P_r ]
    
    # Test-values
    #q_l   = np.array([ 0.5 , 0.3 , 0.0 ])                                        # Initial conditions on LHS; [rho_l,u_l,P_l]
    #q_r   = np.array([ 0.5 ,-0.3 , 0.0 ])                                        # Initial conditions on RHS; [rho_r,u_r,P_r]
    #gamma = 1.4
    #t     = 0.2
    
    # Redefine initial parameters
    rho_l   = q_l[0]; rho_r   = q_r[0]
    u_l     = q_l[1]; u_r     = q_r[1]
    P_l     = q_l[2]; P_r     = q_r[2]
    rhotr_l = 0.0;    rhotr_r = 0.0
    if np.size(q_l) == 4:
        rhotr_l = q_l[3]
    if np.size(q_r) == 4:
        rhotr_r = q_r[3]
    
    
    # Calculate sound speeds
    c_l   = np.sqrt(gamma*P_l/rho_l)
    c_r   = np.sqrt(gamma*P_r/rho_r)
    
    # Define beta-parameter
    beta  = (gamma+1.0)/(gamma-1.0)
    
    # Define pressure values used for iterative solution
    P_n   = 500
    P_i   = np.linspace(0.0,10.0,num=P_n)
    
    x_n       = 1000
    x_start   = 0.0
    x_stop    = 1.0
    x_ghost   = (x_stop-x_start)/5
    x_0       = 0.5
    x         = np.linspace(x_start-x_ghost,x_stop+x_ghost,num=x_n)
    idx_start = np.argmin(abs(x-x_start))
    idx_stop  = np.argmin(abs(x-x_stop))
    
    rho   = np.zeros(x_n)
    u     = np.zeros(x_n)
    P     = np.zeros(x_n)
    rhotr = np.zeros(x_n)
    
    # Construct Hugoniot loci and integral curves given the initial conditions
    Hug_l = u_l + 2.0*c_l/np.sqrt(2.0*gamma*(gamma-1.0)) * ( (1.0-P_i/P_l)/np.sqrt(1.0+beta*P_i/P_l) )
    Hug_r = u_r - 2.0*c_r/np.sqrt(2.0*gamma*(gamma-1.0)) * ( (1.0-P_i/P_r)/np.sqrt(1.0+beta*P_i/P_r) )
    int_l = u_l + 2.0*c_l/(gamma-1.0) * ( 1.0 - (P_i/P_l)**((gamma-1.0)/(2.0*gamma)) )
    int_r = u_r - 2.0*c_r/(gamma-1.0) * ( 1.0 - (P_i/P_r)**((gamma-1.0)/(2.0*gamma)) )
    Phi_l = np.array([Hug_l,int_l])
    Phi_r = np.array([Hug_r,int_r])
    
    # Find intersections between Phi_l and Phi_r
    for l in range(2):
        for r in range(0,2):
            
            idx_s  = np.argmin(abs(Phi_l[l,:]-Phi_r[r,:]))
            P_temp = P_i[idx_s]
            
            # Check that the intersection is physical
            if P_temp*(-1.0)**l >= P_l*(-1.0)**l and P_temp*(-1.0)**r >= P_r*(-1.0)**r:
                
                # Save physically allowed intermediate values of P and u
                P_s = P_temp
                u_s = Phi_l[l,idx_s]
                
                # Location of 2-wave contact discontinuity
                x_2   = x_0 + u_s*t
                idx_2 = np.argmin(abs(x-x_2))
                
                # Passive tracer
                rhotr[0:idx_2] = rhotr_l; rhotr[idx_2:]    = rhotr_r
                
                # Calculate specifics of 1-wave
                if l == 0: # If 1-wave is a 1-shock
                    
                    # Intermediate value of rho
                    rho_sl = ((1.0+beta*P_s/P_l)/(P_s/P_l+beta)) * rho_l
                    
                    # Location of 1-shock
                    x_1    = x_0 + (rho_l*u_l-rho_sl*u_s)/(rho_l-rho_sl) * t
                    idx_1  = np.argmin(abs(x-x_1))
                    
                    # Solution between LHS and contact discontinuity
                    rho[0:idx_1] = rho_l; rho[idx_1:idx_2] = rho_sl
                    u[0:idx_1]   = u_l;   u[idx_1:idx_2]   = u_s
                    P[0:idx_1]   = P_l;   P[idx_1:idx_2]   = P_s
                    
                elif l == 1: # If 1-wave is a rarefaction wave
                    
                    # Intermediate value of rho
                    rho_sl = (P_s/P_l)**(1.0/gamma) * rho_l
                    
                    # Intermediate value of sound speed
                    c_sl   = np.sqrt(gamma*P_s/rho_sl)
                    
                    # Compute values of rho, u, P in 1-rarefaction
                    u_rare   = ( (gamma-1.0)*u_l+2.0*(c_l+(x-x_0)/t) )/( gamma+1.0 )
                    rho_rare = ( ( rho_l**gamma * (u_rare-(x-x_0)/t)**2.0 )/( gamma*P_l ) )**(1.0/(gamma-1.0))
                    P_rare   = ( (P_l)/(rho_l**gamma) ) * (rho_rare)**gamma
                    
                    # Location of 1-rarefaction
                    x_11   = x_0 + (u_l-c_l)*t     # Left edge of rarefaction wave
                    idx_11 = np.argmin(abs(x-x_11))
                    x_12   = x_0 + (u_s-c_sl)*t    # Right edge of rarefaction wave
                    idx_12 = np.argmin(abs(x-x_12))
                    
                    # Solution between LHS and contact discontinuity
                    rho[0:idx_11] = rho_l; rho[idx_11:idx_12] = rho_rare[idx_11:idx_12]; rho[idx_12:idx_2] = rho_sl
                    u[0:idx_11]   = u_l;   u[idx_11:idx_12]   = u_rare[idx_11:idx_12];   u[idx_12:idx_2]   = u_s
                    P[0:idx_11]   = P_l;   P[idx_11:idx_12]   = P_rare[idx_11:idx_12];   P[idx_12:idx_2]   = P_s
                
                # Calculate specifics of 3-wave
                if r == 0: # If 3-wave is a 3-shock
                
                    # Intermediate value of rho
                    rho_sr = ((1.0+beta*P_s/P_r)/(P_s/P_r+beta)) * rho_r
                    
                    # Location of 3-shock
                    x_3    = x_0 + (rho_sr*u_s-rho_r*u_r)/(rho_sr-rho_r) * t
                    idx_3  = np.argmin(abs(x-x_3))
                    
                    # Solution between RHS and contact discontinuity
                    rho[idx_2:idx_3] = rho_sr; rho[idx_3:] = rho_r
                    u[idx_2:idx_3]   = u_s;    u[idx_3:]   = u_r
                    P[idx_2:idx_3]   = P_s;    P[idx_3:]   = P_r
                    
                elif r == 1: # If 3-wave is a rarefaction wave
                    
                    # Intermediate value of rho
                    rho_sr = (P_s/P_r)**(1.0/gamma) * rho_r
                    
                    # Intermediate value of sound speed
                    c_sr   = np.sqrt(gamma*P_s/rho_sr)
                    
                    # Compute values of rho, u, P in 3-rarefaction
                    u_rare   = ( (gamma-1.0)*u_r-2.0*(c_r-(x-x_0)/t) )/( gamma+1.0 )
                    rho_rare = ( ( rho_r**gamma * ((x-x_0)/t-u_rare)**2.0 )/( gamma*P_r ) )**(1.0/(gamma-1.0))
                    P_rare   = ( (P_r)/(rho_r**gamma) ) * (rho_rare)**gamma
                    
                    # Location of 3-rarefaction
                    x_31   = x_0 + (u_s+c_sr)*t      # Left edge of rarefaction wave
                    idx_31 = np.argmin(abs(x-x_31))
                    x_32   = x_0 + (u_r+c_r)*t       # Right edge of rarefaction wave
                    idx_32 =  np.argmin(abs(x-x_32))
                    
                    # Solution between RHS and contact discontinuity
                    rho[idx_2:idx_31] = rho_sr; rho[idx_31:idx_32] = rho_rare[idx_31:idx_32]; rho[idx_32:] = rho_r
                    u[idx_2:idx_31]   = u_s;    u[idx_31:idx_32]   = u_rare[idx_31:idx_32];   u[idx_32:]   = u_r
                    P[idx_2:idx_31]   = P_s;    P[idx_31:idx_32]   = P_rare[idx_31:idx_32];   P[idx_32:]   = P_r
        
    
    s = ( ( np.log(P)-np.log(rho)*gamma )*rho )/( gamma-1.0 )
    p = u*rho
    q = np.array([rho,s,p,P,u,rhotr,x])
    q = q[:,idx_start:idx_stop]
    return q