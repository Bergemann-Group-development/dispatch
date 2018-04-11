import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
from read_patches import *
import dispatch_utils as du
from dusty_wave import dusty_wave

pl.cla()
pl.clf()
pl.close('all')

datadir  = '../data'          # Data directory
run      = '.'
#iout    = [0,5,10,15,20]           # Indices for output times 
iout     = [0,5,10]                    # Indices for output times 
colr     = ['b','c','g','y','r']    # Colorscheme for different output times
#iout    = [0,1,2]                   # Indices for output times
#colr    = ['m','g','y']                # Colorscheme for different output times
verbose  = False
advect   = True
analytic = False
gamma    = 1.4                      # Value of adiabatic exponent
q_l      = [ [1.0,0.0,1.0,0.5] , [1.0,0.0,0.0,0.5] ]
                                    # Initial conditions on LHS; [rho_l,u_l,P_l]
q_r      = [ [0.125,0.0,0.1,0.5] , [0.125,0.0,0.0,0.5] ]
                                    # Initial conditions on RHS; [rho_r,u_r,P_r]
npatches = 5                        # Number of patches in output data
nv       = 6                        # Number of variables per component in this particular case
nvs      = 13                       # Number of variables per component in generic case



# data structure: p.data[x,y,z,params]
# rho  s  p_x  p_y  p_z  (B_x)  (B_y)  (B_z) (rho_t)
#  0   1   2    3    4     5      6      7    5 / 8

# additional data: cycles for second parameters
# P  v_x  v_y  v_z  
# 9  10   11   12

#===============================================================================

# data structure for advect, two components, no mhd: p.data[x,y,z,params]
# rho  s  p_x  p_y  p_z  r_t 
#  0   1   2    3    4    5

# additional data: cycles for second parameters
# P  v_x  v_y  v_z  
# 9  10   11   12


plparams = np.array([10,0,10+nvs,0+nvs])
nplots   = len(plparams)


# Reading the files
#files = list_info_files(datadir)      # list of all metadata files
#pp    = []
pl.figure()


for j in range(nplots):
    pl.subplot(np.ceil(np.sqrt(nplots)),np.ceil(np.sqrt(nplots)),j+1)
    colridx = 0
    xmin    = 100
    xmax    = 0
    ymin    = 100
    ymax    = 0
    qidx    = -1
    pl.xlabel(r'$x$ in code units')
    for i in iout:
        first_patch = True
        for k in range(npatches):
            p = du.patch(id=k+1,iout=i,run=run,data=datadir)
            li = p.li[0] - 1
            ui = p.ui[0]
            
            if plparams[j]==0 or plparams[j]==0+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,0+nv*c2]
                qidx   = 0
                pl.ylabel(r'$\rho$ in code units')
            elif plparams[j]==1 or plparams[j]==1+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,1+nv*c2]
                qidx   = 1
                pl.ylabel(r'$s$ in code units')
            elif plparams[j]==2 or plparams[j]==2+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,2+nv*c2]
                qidx   = 2
                pl.ylabel(r'$p_x$ in code units')
            elif plparams[j]==3 or plparams[j]==3+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,3+nv*c2]
                pl.ylabel(r'$p_y$ in code units')
            elif plparams[j]==4 or plparams[j]==4+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,4+nv*c2]
                pl.ylabel(r'$p_z$ in code units')
            elif plparams[j]==5 or plparams[j]==5+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,5+nv*c2]
                if advect:
                    pl.ylabel(r'$\rho_{tracer}$ in code units')
                    qidx = 5
                else:
                    pl.ylabel(r'$B_x$ in code units')
            elif plparams[j]==6 or plparams[j]==6+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,6+nv*c2]
                pl.ylabel(r'$B_y$ in code units')
            elif plparams[j]==7 or plparams[j]==7+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,7+nv*c2]
                pl.ylabel(r'$B_z$ in code units')
            elif plparams[j]==8 or plparams[j]==8+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,8+nv*c2]
                pl.ylabel(r'$\rho_{tracer}$ in code units')
                qidx = 5
            elif plparams[j]==9 or plparams[j]==9+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = np.exp( p.data[:,0,0,1+nv*c2]*(gamma-1.0)/p.data[:,0,0,0+nv*c2] + np.log(p.data[:,0,0,0+nv*c2])*gamma ) \
                         * ( np.sum(abs(p.data[:,0,0,1+nv*c2])) > 0.0 )
                qidx   = 3
                pl.ylabel(r'$P$ in code units')
            elif plparams[j]==10 or plparams[j]==10+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,2+nv*c2]/p.data[:,0,0,0+nv*c2]
                qidx   = 4
                pl.ylabel(r'$v_x$ in code units')
            elif plparams[j]==11 or plparams[j]==11+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,3+nv*c2]/p.data[:,0,0,0+nv*c2]
                pl.ylabel(r'$v_y$ in code units')
            elif plparams[j]==12 or plparams[j]==12+nvs:
                c2     = ( plparams[j] >= nvs )
                pldata = p.data[:,0,0,4+nv*c2]/p.data[:,0,0,0+nv*c2]
                pl.ylabel(r'$v_z$ in code units')
            else:
                pldata = p.data[:,0,0,0]
            
            
            
            
            if first_patch:
                pl.plot( p.x[li:ui] , pldata , colr[colridx] , label = r'$t$ = '+np.str(round(p.time,4)) )
                first_patch = False
            else:
                pl.plot( [p.x[li] , p.x[li]] , [-100.0,100.0] , 'k:' , lw=0.7)
                pl.plot( p.x[li:ui] , pldata , colr[colridx] )
            
            ymin = min([ymin,pldata.min()])
            ymax = max([ymax,pldata.max()])
            xmin = min([xmin,p.x[li]])
            xmax = max([xmax,p.x[ui]])
        
        #q  = [ rho , s , p_x , P , v_x , rhotr , x ]
        if analytic:
            q   = shock_tube(q_l[c2],q_r[c2],gamma,p.time)
            if qidx != -1:
                pl.plot( q[6,:] , q[qidx,:] , '--'+colr[colridx] )
        colridx += 1
    
    pl.axis( xmin = xmin , xmax = xmax , ymin = ymin-(ymax*1.1-ymax) , ymax = 0+ymax*1.1)
    pl.legend()



pl.show()