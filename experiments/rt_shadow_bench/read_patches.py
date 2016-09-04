import mayavi.mlab as ml
from mayavi.tools import pipeline as pl
import itertools
import numpy as np
from dispatch_data import *

def draw_cube_borders(fig,p):
    """Draw an outline of a Cartesian patch."""

    # flatten the tuple of corners
    llc = [x-y for x,y in zip(p.pos,p.dx)]
    cc = [x for y in p.active_corners for x in y]
    vertices = set(itertools.permutations(cc,3))

    for v,w in itertools.product(vertices,vertices):
        x0, y0, z0 = v
        x1, y1, z1 = w
        s = sum(bool(x) for x in [x0 == x1, y0==y1, z0==z1])
        if s == 2:
            ml.plot3d([x0+llc[0],x1+llc[0]], [y0+llc[1],y1+llc[1]], [z0+llc[1],z1+llc[1]],
                      figure=fig, tube_radius=None, line_width=0.5)

    print p.id

if __name__ == "__main__":

    datadir = 'data/'
    prefix = ''
    dumps = [0]
    ThisPatch = StaggerMHDPatch

    m = 0 # counter
    for iout in dumps:
        files = fnmatch.filter(os.listdir(datadir), prefix + "*_{0:04d}.txt".format(iout))

        patches = []
        times = []
        for f in files:
            patches.append(ThisPatch(datadir+f.strip('.txt'),verbose=False,read_derivs=False))

        for p in patches:
            if p.task_type == 'immersed_boundary_t':
                print "Skipping ",p.id
                continue
            if p.id == 1:
                root = p

            print p.filename, p.time
            times.append(p.time)

        m += 1
    #fig = ml.figure('3d', size=(700,700), fgcolor=(1,1,1), bgcolor=(0.5,0.5,0.5))
    #ml.clf(fig)

    #for p in patches:
    #    llc = [x-y for x,y in zip(p.pos,p.dx)]
    #    var = p.data[:,:,:,p.varidx['d']].squeeze()
    #    x_i = dispatch_grid(p)
    #    x, y, z = np.meshgrid(x_i[0]+llc[0],x_i[1]+llc[1],x_i[2]+llc[2],indexing='ij')
    #    pts = ml.points3d(x,y,z,var,colormap='jet')
    #    draw_cube_borders(fig,p)
