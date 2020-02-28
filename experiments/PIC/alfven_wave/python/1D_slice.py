# -*- coding: utf-8 -*-
'''
    Show the time evolution in a sequence of panels
'''
import numpy as np
import matplotlib.pyplot as pl
import dispatch
import os
from sys import argv

def one_dim(pp,iv,vmin,vmax,varname):
    num_patches = len(pp.patches)
    n=pp.patches[0].n
    ng=pp.patches[0].ng
    box=pp.patches[0].box
    x_center=int(box[0]*n[0]/2 -1)
    y_center=int(box[1]*n[1]/2 -1)
    pos_arr = np.zeros([int(box[0]),int(box[1]),int(box[2])])
    patch_arr = np.zeros([np.shape(pos_arr)[0]*n[0], np.shape(pos_arr)[1]*n[1], np.shape(pos_arr)[2]*n[2]])
    for ip in range(num_patches):
        p = pp.patches[ip]
        position = [int(p.position[0]-0.5), int(p.position[1]-0.5), int(p.position[2]-0.5)]
        pos_arr[position[0], position[1], position[2]] = p.id
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        no_gz=data[ng[0]:-ng[0],ng[1]:-ng[1],ng[2]:-ng[2]]
        patch_arr[position[0]*n[0]:position[0]*n[0]+n[0],position[1]*n[1]:position[1]*n[1]+n[1],position[2]*n[2]:position[2]*n[2]+n[2]] = data #no_gz#.transpose()
    final = patch_arr[x_center,y_center,:].transpose()
    #pl.imshow(final, cmap='jet')#,vmin=vmin, vmax=vmax)
    pl.plot(final,'.')
    pl.title(varname); 
    pl.plot([31.5,31.5],[-10,10],'k-')
    pl.plot([63.5,63.5],[-10,10],'k-')
    pl.ylim([min(final)-0.2*min(final),max(final)+0.2*max(final)])
    #pl.xlim([int(n[0])-5,2*int(n[0])+5])
    #pl.colorbar()

def one_dim_slices(data,run,save_plots=False,ylim0=0,ylim1=3):
    "Plot d, e, px, py for a slice through y"
    # Make a list of the snapshots
    ss = []
    start = 0
    end = sum(os.path.isdir(data+'/'+run+'/'+i) for i in os.listdir(data+'/'+run)) - 1
    step = 1
    for i in range(start,end,step):
        s = dispatch.snapshot(i,run=run,data=data,copy=True)
        ss.append(s)



    pl.figure(1, figsize=(11,7))
    for it in range(start,end,step):
        print('it=',it)
        pl.clf()
        pl.suptitle('iout={}                            time={:5.1f}'.format(it,ss[it].time))
        pl.subplot(2,2,1); one_dim(ss[it], 0, 0., 1.1 ,'d')#,ylim0,ylim1)
        pl.subplot(2,2,2); one_dim(ss[it], 7, 0.00, 0.005, 'bz')#,ylim0,ylim1)
        pl.subplot(2,2,3); one_dim(ss[it], 2, .008, .0, 'py')#,ylim0,ylim1)
        pl.subplot(2,2,4); one_dim(ss[it], 3, -.20, +.20, 'pz')#,ylim0,ylim1)
        if it==0:
            pl.tight_layout()
        pl.draw()
        pl.pause(0.001)
        if save_plots == True:
            pl.savefig("figures/{}/{:03d}_1D".format(run,it))

def run():
    data = '../data'
    run = 'input27'
    a = one_dim_slices(data,run,True)

if __name__=='__main__':
    run()
