import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
import os

pl.cla()
pl.clf()
pl.close('all')

datadir  = 'data'          # Data directory
outdir   = 'snap_shots'    # Output directory containing snap shots
format   = '.pdf'
run      = '.'
iend     = 2
colr     = ['b','c','g','y','r']    # Colorscheme for different output times
#iout    = [0,1,2]                  # Indices for output times
#colr    = ['m','g','y']            # Colorscheme for different output times
verbose      = False
advect       = True
show_patches = True
gamma    = 1.0001 #5.0/3.0
Mstar    = 1.0
aplanet  = 1.0
sigma0   = 0.002 * Mstar/(np.pi*aplanet**2)
patches  = [1,5,1]                  # Patch structure in output data
npatches = np.prod(patches)         # Number of patches in output data
nv       = 5                        # Number of variables per component in this particular case
nvs      = 14                       # Number of variables per component in generic case
plot_dim = '2D'                     # Dimension 



# data structure: p.data[x,y,z,params]
# rho  s  p_x  p_y  L_z  (B_x)  (B_y)  (B_z) (rho_t)
#  0   1   2    3    4     5      6      7    5 / 8

# additional data: cycles for second parameters
# P  v_x  v_y  v_z  p_z
# 9  10   11   12   13

#===============================================================================

# data structure for advect, two components, no mhd: p.data[x,y,z,params]
# rho  s  p_x  p_y  L_z  r_t 
#  0   1   2    3    4    5

# additional data: cycles for second parameters
# P  v_x  v_y  v_z  p_z
# 9  10   11   12   13


plparams = np.array([0])
nplots   = len(plparams)

# Determine min and max values
xmin    =  100
xmax    = -100
ymin    =  100
ymax    = -100
zmin    =  100
zmax    = -100
Vmin = np.zeros(nplots)+100.; Vmax = np.zeros(nplots)-100.
first_time = True

# Create output directory if it is not present
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Start figure generation
for ii in range(0,2): # The first cycle determines min and max values, the second cycle creates the figures
    for i in range(0,iend+1):
        if first_time == False:
            fig=pl.figure(figsize=[10,6])
        k = 0
        for j in range(nplots):
            if first_time == False:
                pl.subplot(np.ceil(np.sqrt(nplots)),np.ceil(np.sqrt(nplots)),j+1)
                pl.xlabel(r'$r$ in code units')
            colridx =  0
            for iz in range(patches[2]):
                for iy in range(patches[1]):
                    for ix in range(patches[0]):
                        p  = du.patch(id=k+1,iout=i,run=run,data=datadir)
                        li = p.li[:]
                        ui = p.ui[:]
                        
                        R,Z,PHI = np.meshgrid(p.y,p.x,p.z)
                        
                        if plparams[j]==0 or plparams[j]==0+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = np.log10(p.data[:,:,:,0+nv*c2]/sigma0)
                            unit   = r'log($\rho / \rho _0$) in code units'
                        elif plparams[j]==1 or plparams[j]==1+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,1+nv*c2]
                            unit   = r'$s$ in code units'
                        elif plparams[j]==2 or plparams[j]==2+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,2+nv*c2]
                            unit   = r'$p_z$ in code units'
                        elif plparams[j]==3 or plparams[j]==3+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,3+nv*c2]
                            unit   = r'$p_r$ in code units'
                        elif plparams[j]==4 or plparams[j]==4+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,4+nv*c2]
                            unit   = r'$L_{\phi}$ in code units'
                        elif plparams[j]==5 or plparams[j]==5+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,10+nv*c2]
                            if advect:
                                unit = r'$\rho_{tracer}$ in code units'
                            else:
                                unit = r'$B_z$ in code units'
                        elif plparams[j]==6 or plparams[j]==6+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,6+nv*c2]
                            unit   = r'$B_r$ in code units'
                        elif plparams[j]==7 or plparams[j]==7+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,7+nv*c2]
                            unit   = r'$B_{\phi}$ in code units'
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
                            pldata = p.data[:,:,:,2+nv*c2]/p.data[:,:,:,0+nv*c2]
                            unit   = r'$v_z$ in code units'
                        elif plparams[j]==11 or plparams[j]==11+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,3+nv*c2]/p.data[:,:,:,0+nv*c2]
                            unit   = r'$v_r$ in code units'
                        elif plparams[j]==12 or plparams[j]==12+nvs:
                            c2     = ( plparams[j] >= nvs )
                            pldata = p.data[:,:,:,4+nv*c2]/R/p.data[:,:,:,0+nv*c2]
                            unit   = r'$v_{\phi}$ in code units'
        		elif plparams[j]==13 or plparams[j]==13+nvs:
 			    c2     = ( plparams[j] >= nvs )
 			    pldata = p.data[:,:,:,4+nv*c2]/R
 			    unit   = r'$p_{\phi}$ in code units'
                        else:
                            pldata = p.data[:,:,:,0]
                            unit   = r'$\rho$ in code units'
                        
                        if first_time == True:
                            ymin  = min([ymin,p.y[li[1]]])
                            ymax  = max([ymax,p.y[ui[1]-1]])
                            Vminp = float(np.min(pldata)); Vmin[j] = min(Vmin[j],Vminp)
                            Vmaxp = float(np.max(pldata)); Vmax[j] = max(Vmax[j],Vmaxp)
                        if plot_dim=='2D' and first_time == False:
                            if ix+1 == (patches[0]+1)/2:
                                V    = np.linspace(Vmin[j],Vmax[j],100)
                                Y,Z = np.meshgrid(p.y[li[1]:ui[1]+1],p.z[li[2]:ui[2]+1])
                                pl.contourf(Y,Z,pldata[(ui[0]-li[0])/2+li[0],li[1]:ui[1]+1,li[2]:ui[2]+1].transpose(),V)
                                if (iy+1 > 1) and show_patches:
                                    pl.plot( [p.y[li[1]] , p.y[li[1]]] , [p.z[li[2]],p.z[ui[2]-1]] , 'k:' , lw=0.7)
                                if (iz+1 > 1) and show_patches:
                                    pl.plot( [p.y[li[1]] , p.y[ui[1]-1]] , [p.z[li[2]],p.z[li[2]]] , 'k:' , lw=0.7)
                                    
                                pl.ylabel(r'$\phi$ in code units')
                                zmin = min([zmin,p.z[li[2]]])
                                zmax = max([zmax,p.z[ui[2]-1]])
                        elif plot_dim=='1D' and first_time == False:
                            if ix+1 == (patches[0]+1)/2 and iz+1 == (patches[2]+1)/2:
                                pl.plot( p.y[li[1]:ui[1]+1] , pldata[(ui[0]-li[0])/2+li[0] , li[1]:ui[1]+1, (ui[2]-li[2])/2+li[2]].transpose() , colr[colridx] )
                                pl.plot( p.y[li[1]:ui[1]+1] , pldata[(ui[0]-li[0])/2+li[0] , li[1]:ui[1]+1, (ui[2]-li[2])/2+li[2]].transpose() , 'o'+colr[colridx], label = r'$t$ = '+np.str(round(p.time,4)))
                                if (k+1 > 1) and show_patches:
                                    pl.plot( [p.y[li[1]] , p.y[li[1]]] , [-100.0,1000.0] , 'k:' , lw=0.7)
                                
                                zmin = min([zmin,pldata.min()])
                                zmax = max([zmax,pldata.max()])
                                pl.ylabel(unit)
                        
                        k = k+1
            colridx = colridx + 1
            
            if first_time == False:
                if plot_dim == '2D':
                    pl.axis( xmin=ymin , xmax=ymax , ymin=float(zmin), ymax=float(zmax) )
                    cbar = pl.colorbar()
                    cbar.ax.set_ylabel(unit)
                elif plot_dim == '1D':
                    pl.axis( xmin = ymin , xmax = ymax , ymin = Vmin , ymax = Vmax)
                    pl.legend()
        if first_time == False:
            fig.savefig(outdir+'/snap_'+np.str(i)+format)
            print 'Saved snap shot '+np.str(i)+' in '+outdir
    first_time = False
    

#pl.show()
print 'Done.'
