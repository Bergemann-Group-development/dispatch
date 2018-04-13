import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as pl
#from read_patches import *
import dispatch_utils as du
import scipy.interpolate as sc

pl.cla()
pl.clf()
pl.close('all')

datadir  = '../data'          # Data directory
run      = '.'
#iout    = [0,5,10,15,20]           # Indices for output times 
iout     = [5]                     # Indices for output times 
colr     = ['b','c','g','y','r']    # Colorscheme for different output times
#iout    = [0,1,2]                  # Indices for output times
#colr    = ['m','g','y']            # Colorscheme for different output times
verbose      = False
advect       = True
show_patches = True
gamma    = 5.0/3.0
patches  = [3,3,3]                  # Patch structure in output data
npatches = np.prod(patches)         # Number of patches in output data
nv       = 6                        # Number of variables per component in this particular case
nvs      = 13                       # Number of variables per component in generic case
plot_dim = '2D'                     # Dimension 



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


plparams = np.array([0,0+nvs,2,12])
nplots   = len(plparams)


# Reading the files
#files = list_info_files(datadir)      # list of all metadata files
#pp    = []
pl.figure()


if plot_dim == '2D':
    iout = [iout[0], iout[0]]
for j in range(nplots):
    pl.subplot(np.ceil(np.sqrt(nplots)),np.ceil(np.sqrt(nplots)),j+1)
    colridx =  0
    xmin    =  100
    xmax    = -100
    ymin    =  100
    ymax    = -100
    zmin    =  100
    zmax    = -100
    Vmin    =  100
    Vmax    = -100
    first_time = True
    dx      = 0.0
    dy      = 0.0
    pl.xlabel(r'$x$ in code units')
    for i in iout:
        k   =  0
        for iz in range(patches[2]):
            for iy in range(patches[1]):
                for ix in range(patches[0]):                   
                    p  = du.patch(id=k+1,iout=i,run=run,data=datadir)
                    li = p.li[:] - 1
                    ui = p.ui[:]
                    
                    if plparams[j]==0 or plparams[j]==0+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,0+nv*c2]
                        unit   = r'$\rho$ in code units'
                    elif plparams[j]==1 or plparams[j]==1+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,1+nv*c2]
                        unit   = r'$s$ in code units'
                    elif plparams[j]==2 or plparams[j]==2+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,2+nv*c2]
                        unit   = r'$p_x$ in code units'
                        dx     = p.dx[0]/2.0
                    elif plparams[j]==3 or plparams[j]==3+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,3+nv*c2]
                        unit   = r'$p_y$ in code units'
                    elif plparams[j]==4 or plparams[j]==4+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,4+nv*c2]
                        unit   = r'$p_z$ in code units'
                    elif plparams[j]==5 or plparams[j]==5+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,5+nv*c2]
                        if advect:
                            unit = r'$\rho_{tracer}$ in code units'
                        else:
                            unit = r'$B_x$ in code units'
                        dx     = p.dx[0]/2.0
                    elif plparams[j]==6 or plparams[j]==6+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,6+nv*c2]
                        unit   = r'$B_y$ in code units'
                    elif plparams[j]==7 or plparams[j]==7+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,7+nv*c2]
                        unit   = r'$B_z$ in code units'
                    elif plparams[j]==8 or plparams[j]==8+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,:,8+nv*c2]
                        unit   = r'$\rho_{tracer}$ in code units'
                    elif plparams[j]==9 or plparams[j]==9+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = np.exp( p.data[:,:,:,1+nv*c2]*(gamma-1.0)/p.data[:,:,:,0+nv*c2] + np.log(p.data[:,:,:,0+nv*c2])*gamma ) \
                                * ( np.sum(abs(p.data[:,:,:,1+nv*c2])) > 0.0 )
                        unit   = r'$P$ in code units'
                    elif plparams[j]==10 or plparams[j]==10+nvs:
                        c2     = ( plparams[j] >= nvs )
                        dx     = p.dx[0]/2.0
                        pldata = p.data[0:-1,:,:,0+nv*c2]/2.0 + p.data[1:,:,:,0+nv*c2]/2.0
                        pldata = np.insert(pldata, 0, 1, axis=0)
                        pldata = p.data[:,:,:,2+nv*c2]/pldata
                        unit   = r'$v_x$ in code units'
                    elif plparams[j]==11 or plparams[j]==11+nvs:
                        c2     = ( plparams[j] >= nvs )
                        dy     = p.dx[1]/2.0
                        pldata = p.data[:,0:-1,:,0+nv*c2]/2.0 + p.data[:,1:,:,0+nv*c2]/2.0
                        pldata = np.insert(pldata, 0, 1, axis=1)
                        pldata = p.data[:,:,:,3+nv*c2]/p.data[:,:,:,0+nv*c2]
                        unit   = r'$v_y$ in code units'
                    elif plparams[j]==12 or plparams[j]==12+nvs:
                        c2     = ( plparams[j] >= nvs )
                        pldata = p.data[:,:,0:-1,0+nv*c2]/2.0 + p.data[:,:,1:,0+nv*c2]/2.0
                        pldata = np.insert(pldata, 0, 1, axis=2)
                        pldata = p.data[:,:,:,4+nv*c2]/p.data[:,:,:,0+nv*c2]
                        unit   = r'$v_z$ in code units'
                    else:
                        pldata = p.data[:,:,:,0]
                        unit   = r'$\rho$ in code units'
                    
                    
                    xmin  = min([xmin,p.x[li[0]+1]])
                    xmax  = max([xmax,p.x[ui[0]]])
                    Vminp = float(np.min(pldata)); Vmin = min(Vmin,Vminp)
                    Vmaxp = float(np.max(pldata)); Vmax = max(Vmax,Vmaxp)
                    print p.id, p.time, p.dt, p.x[p.li[0]]
                    if plot_dim=='2D' and first_time==False:
                        if iz+1 == (patches[2]+1)/2:
                            V    = np.linspace(Vmin,Vmax,100)
                            X,Y = np.meshgrid(p.x[li[0]:ui[0]]-dx,p.y[li[1]:ui[1]]-dy)
                            pl.contourf(X,Y,pldata[li[0]:ui[0],li[1]:ui[1],(ui[2]-li[2])/2+li[2]].transpose(),V)
                            if (ix+1 > 1) and show_patches:
                                pl.plot( [p.x[li[0]+1] , p.x[li[0]+1]] , [p.y[li[1]],p.y[ui[1]-1]] , 'k:' , lw=0.7)
                            if (iy+1 > 1) and show_patches:
                                pl.plot( [p.x[li[0]] , p.x[ui[0]-1]] , [p.y[li[1]+1],p.y[li[1]+1]] , 'k:' , lw=0.7)
                                
                            pl.ylabel(r'$y$ in code units')
                            ymin = min([ymin,p.y[li[1]]])
                            ymax = max([ymax,p.y[ui[1]-1]])
                            pl.title(unit)
                    elif plot_dim=='1D':
                        if iz+1 == (patches[2]+1)/2 and iy+1 == (patches[1]+1)/2:
                            pl.plot( p.x[li[0]:ui[0]]-dx , pldata[li[0]:ui[0] , (ui[1]-li[1])/2+li[1],(ui[2]-li[2])/2+li[2]].transpose() , 'b' )
                            pl.plot( p.x[li[0]:ui[0]]-dx , pldata[li[0]:ui[0] , (ui[1]-li[1])/2+li[1],(ui[2]-li[2])/2+li[2]].transpose() , 'bo' )
                            if (k+1 > 1) and show_patches:
                                pl.plot( [p.x[li[0]+1] , p.x[li[0]+1]] , [-100.0,1000.0] , 'k:' , lw=0.7)
                            
                            ymin = min([ymin,pldata.min()])
                            ymax = max([ymax,pldata.max()])
                            pl.ylabel(unit)
                    
                    k = k+1
        first_time = False
    
    pl.axis( xmin=xmin , xmax=xmax , ymin=ymin, ymax=ymax )
    if plot_dim == '2D':
        pl.colorbar()
    elif plot_dim == '1D':
        pl.axis( xmin = xmin , xmax = xmax , ymin = ymin-(ymax*1.1-ymax) , ymax = 0+ymax*1.1)

    
    



pl.show()
