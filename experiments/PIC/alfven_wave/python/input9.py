# -*- coding: utf-8 -*-
'''
    Show the time evolution in a sequence of panels
'''

import numpy as np
import matplotlib.pyplot as pl
import dispatch
import os
#import sys
#sys.path.insert(0,'../../../utilities/python')
#import dispatch_utils as du
#import plot_utils as pu
#reload(du)

data = '../data'
run = 'input9_mhd'

def im3(pp,iv,vmin,vmax,varname): #,ylim0=0,ylim1=3):
    #pl.xlim(0,1)#; pl.ylim(ylim0,ylim1)
    num_patches = len(pp.patches)
    n=pp.patches[0].n
    g=pp.patches[0].ng
    print(g)
    m=int(n[1]/2 -1)
    pos_arr = np.zeros([int(pp.patches[0].box[0]),int(pp.patches[0].box[1]),int(pp.patches[0].box[2])])
    #print(np.shape(pos_arr)[0])
    patch_arr = np.zeros([np.shape(pos_arr)[0]*n[0], np.shape(pos_arr)[1]*n[1], np.shape(pos_arr)[2]*n[2]])
    for ip in range(num_patches):
        #print(pp.patches[ip].position[0])
        #print(pp.patches[ip].position[1])
        #print(pp.patches[ip].position[2])
        p = pp.patches[ip]
        #print(p.position[0], type(p.position[0]))
        position = [int(p.position[0]-0.5), int(p.position[1]-0.5), int(p.position[2]-0.5)]
        pos_arr[position[0], position[1], position[2]] = p.id
        #print(p.id)
        #print(p.box[0])
        #print(p.box[1])
        #print(p.box[2])
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        #no_gz=np.sum(data[g[0]:-g[0],g[1]:-g[1],g[2]:-g[2]],axis=1)
        no_gz=data[g[0]:-g[0],g[1]:-g[1],g[2]:-g[2]]
        patch_arr[position[0]*n[0]:position[0]*n[0]+n[0],position[1]*n[1]:position[1]*n[1]+n[1],position[2]*n[2]:position[2]*n[2]+n[2]] = no_gz#.transpose()
    #print(pos_arr)
    #print(np.shape(patch_arr.transpose()))
    #final = np.sum(patch_arr,axis=1).transpose()#.reshape(n[0]*3,n[0]*3)
    final = patch_arr[:,15,:].transpose()
    pl.imshow(final, cmap='jet')
    pl.title(varname); 
    pl.colorbar()

def y_slices(datadir,save_plots='False',ylim0=0,ylim1=3):
    "Plot d, e, px, py for a slice through y"
    # Make a list of the snapshots
    ss = []
    start = 0
    end = sum(os.path.isdir(data+'/'+run+'/'+i) for i in os.listdir(data+'/'+run+'/')) - 3 
    #print(end)
    step = 1
    for i in range(start,end,step):
        s = dispatch.snapshot(i,run=run,data=data)
        ss.append(s)



#    pl.ion()
    pl.figure(1, figsize=(11,7))
#    pl.show()
    for it in range(start,end,step):
        print('it=',it)
        pl.clf()
        pl.suptitle('iout={}                            time={:5.1f}'.format(it,ss[it].time))
        #print "e, max,min,mean", np.max(pp[1].data[:,:,:,0]), np.min(pp[1].data[:,:,:,0]), np.mean(pp[1].data[:,:,:,0])
        #for ip in range(num_patches):
        #    pp = ss[it].patches[ip]
        #    #density = pp.data[0]
        pl.subplot(2,2,1); im3(ss[it], 0, 470, 530 ,'d')#,ylim0,ylim1)
        #pl.subplot(1,4,1); im3(pp, 1, 0.01, 0.08 ,'e' )
        pl.subplot(2,2,2); im3(ss[it], 7, 0.00, 0.08, 'bz')#,ylim0,ylim1)
        pl.subplot(2,2,3); im3(ss[it], 1, -.20, +.20, 'px')#,ylim0,ylim1)
        pl.subplot(2,2,4); im3(ss[it], 2, -.20, +.20, 'py')#,ylim0,ylim1)
        #pl.subplot(1,4,4); im3(pp, 4, -.20, +.20, 'pz')
        if it==0:
            pl.tight_layout()
        pl.draw()
        pl.pause(0.001)
        if save_plots == 'True':
            #pl.savefig('images/border/border_xz_plan'+str(iout))
            pl.savefig("{}/{:03d}".format(run,it))

#def border_images(datadir,save_plots='False',ylim0=0,ylim1=3):
#    pl.ion()
#    pl.figure(1, figsize=(11,7))
#    pl.show()
#    for iout in range(11):
#        pl.clf()
#        pp=du.patches(iout,datadir)    
#        pl.suptitle('iout={}                            time={:5.1f}'\
#        .format(iout,pp[0].time))
#        #print "e, max,min,mean", np.max(pp[1].data[:,:,:,0]), np.min(pp[1].data[:,:,:,0]), np.mean(pp[1].data[:,:,:,0])
#        pl.subplot(2,2,1); im3(pp, 0,  470,  530,  'd', ylim0,ylim1)
#        pl.subplot(2,2,2); im3(pp, 7, 0.00, 0.08, 'bz', ylim0,ylim1)
#        pl.subplot(2,2,3); im3(pp, 2, -.20, +.20, 'px', ylim0,ylim1)
#        pl.subplot(2,2,4); im3(pp, 3, -.20, +.20, 'py', ylim0,ylim1)
#        if iout==0:
#            pl.tight_layout()
#        pl.draw()
#        pl.pause(0.001)
#        if save_plots == 'True':
#            pl.savefig('images/boundary/'+str(iout))

#x = border_images('input3','False',1.5,2.5)
a = y_slices(run,'True')

