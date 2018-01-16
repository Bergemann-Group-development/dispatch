# -*- coding: utf-8 -*-

"""
   Support functions for reading or calculating DISPATCH2 grid and geometry information.
"""

import numpy as np
from abc import abstractproperty, abstractmethod, ABCMeta
from curvilinear_transforms import *
from dispatch_data import Patch
from legacy_reader import mesh_types

class Grid(dict):
    """An abstract grid class based on `dict()`."""

    __metaclass__= ABCMeta

class GeometricFactors(Grid):
    """Calculate and store the geometric factors used by curvilinear grids."""

    def __init__(self):
        """Constructor."""

        # Define geometric factors with the same notation as in `mesh_mod` ("c"
        # for zone-centred and "f" for face-centred).
        self['h2c'] = None
        self['h2f'] = None
        self['h31c'] = None
        self['h31f'] = None
        self['h32c'] = None
        self['h32f'] = None
        self['dx1c'] = None
        self['dx1f'] = None
        self['dx2c'] = None
        self['dx2f'] = None
        self['dx3c'] = None
        self['dx3f'] = None
        self['dvol1c'] = None
        self['dvol1f'] = None
        self['dvol2c'] = None
        self['dvol2f'] = None
        self['dvol3c'] = None
        self['dvol3f'] = None
        self['dar1c'] = None
        self['dar1f'] = None
        self['dar2c'] = None
        self['dar2f'] = None
        self['dar31c'] = None
        self['dar31f'] = None
        self['dar32c'] = None
        self['dar32f'] = None

    def init_grid(self, p, x_i):
        """Initialise geometric factors based on coord. type."""

        if p.mesh_type == 'Cartesian': self.init_Cartesian(p, x_i)
        elif p.mesh_type == 'cylindrical': self.init_cylindrical(p, x_i)
        elif p.mesh_type == 'spherical': self.init_spherical(p, x_i)

    def init_Cartesian(self, p, x_i):
        """Initialise geometric factors for a Cartesian coord. system."""

        n1, n2, n3 = p.nwrite

        # 1-direction
        self['h2c'] = np.ones(n1)
        self['h2f'] = np.ones(n1)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        self['h32c'] = np.ones(n2)
        self['h32f'] = self['h32c'].view()

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.dx[0]
        self['dx1f'] = np.ones(n1) * p.dx[0]
        self['dx2c'] = np.ones(n2) * p.dx[1]
        self['dx2f'] = np.ones(n2) * p.dx[1]
        self['dx3c'] = np.ones(n3) * p.dx[2]
        self['dx3f'] = np.ones(n3) * p.dx[2]

        # volume elements
        self['dvol1c'] = np.ones(n1) * p.dx[0]
        self['dvol1f'] = np.ones(n1) * p.dx[0]
        self['dvol2c'] = np.ones(n2) * p.dx[1]
        self['dvol2f'] = np.ones(n2) * p.dx[1]
        self['dvol3c'] = np.ones(n3) * p.dx[2]
        self['dvol3f'] = np.ones(n3) * p.dx[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.dx[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.dx[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.dx[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.dx[0] / self['dvol1f']
        self['dar32c'] = p.dx[1] / self['dvol2c']
        self['dar32f'] = p.dx[1] / self['dvol2f']

    def init_cylindrical(self, p, x_i):
        """Initialise geometric factors for a cylindrical coord. system."""

        n1, n2, n3 = p.nwrite

        # 1-direction
        self['h2c'] = np.ones(n1)
        self['h2f'] = np.ones(n1)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        pos_c = np.array(x_i[1])
        pos_f = np.array(x_i[1] - 0.5 * p.dx[1])
        self['h32c'] = abs(pos_c)
        self['h32f'] = abs(pos_f)

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.dx[0]
        self['dx1f'] = np.ones(n1) * p.dx[0]
        self['dx2c'] = np.ones(n2) * p.dx[1]
        self['dx2f'] = np.ones(n2) * p.dx[1]
        self['dx3c'] = np.ones(n3) * p.dx[2]
        self['dx3f'] = np.ones(n3) * p.dx[2]

        # volume elements
        self['dvol1c'] = np.ones(n1) * p.dx[0]
        self['dvol1f'] = np.ones(n1) * p.dx[0]
        self['dvol2c'] = np.empty_like(pos_c)
        self['dvol2f'] = np.empty_like(pos_f)
        for j in xrange(len(pos_c)-1):
            self['dvol2c'][j] = 0.5 * abs( self['h32f'][j+1] * pos_f[j+1]
                                         - self['h32f'][j  ] * pos_f[j  ] )
            self['dvol2f'][j] = 0.5 * abs( self['h32c'][j  ] * pos_c[j  ]
                                         - self['h32c'][j-1] * pos_c[j-1] )
        self['dvol3c'] = np.ones(n3) * p.dx[2]
        self['dvol3f'] = np.ones(n3) * p.dx[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.dx[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.dx[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.dx[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.dx[0] / self['dvol1f']
        self['dar32c'] = p.dx[1] / self['dvol2c']
        self['dar32f'] = p.dx[1] / self['dvol2f']

    def init_spherical(self, p, x_i):
        """Initialise geometric factors for a spherical coord. system."""

        n1, n2, n3 = p.nwrite

        # 1-direction
        rpos_c = np.array(x_i[0])
        rpos_f = np.array(x_i[0] - 0.5 * p.dx[0])
        self['h2c'] = abs(rpos_c)
        self['h2f'] = abs(rpos_f)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        tpos_c = np.array(x_i[1])
        tpos_f = np.array(x_i[1] - 0.5 * p.dx[1])
        self['h32c'] = abs(np.sin(tpos_c))
        self['h32f'] = abs(np.sin(tpos_f))

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.dx[0]
        self['dx1f'] = np.ones(n1) * p.dx[0]
        self['dx2c'] = np.ones(n2) * p.dx[1]
        self['dx2f'] = np.ones(n2) * p.dx[1]
        self['dx3c'] = np.ones(n3) * p.dx[2]
        self['dx3f'] = np.ones(n3) * p.dx[2]

        # volume elements
        self['dvol1c'] = np.empty_like(rpos_c)
        self['dvol1f'] = np.empty_like(rpos_f)
        for i in xrange(len(rpos_c)-1):
            self['dvol1c'] = (self['h2f'][i+1] * self['h31f'][i+1] * rpos_f[i+1] / 3.0
                            - self['h2f'][i  ] * self['h31f'][i  ] * rpos_f[i  ] / 3.0)
            self['dvol1f'] = (self['h2c'][i  ] * self['h31c'][i  ] * rpos_c[i  ] / 3.0
                           -  self['h2c'][i-1] * self['h31c'][i-1] * rpos_c[i-1] / 3.0)
        self['dvol2c'] = np.empty_like(tpos_c)
        self['dvol2f'] = np.empty_like(tpos_f)
        for j in xrange(len(tpos_c)-1):
            self['dvol2c'] = -np.cos(tpos_c[j+1]) + np.cos(tpos_c[j  ])
            self['dvol2f'] = -np.cos(tpos_f[j  ]) + np.cos(tpos_f[j-1])
        self['dvol3c'] = np.ones(n3) * p.dx[2]
        self['dvol3f'] = np.ones(n3) * p.dx[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.dx[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.dx[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.dx[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.dx[0] / self['dvol1f']
        self['dar32c'] = p.dx[1] / self['dvol2c']
        self['dar32f'] = p.dx[1] / self['dvol2f']

def dispatch_grid(p, truecoords=True):
    """
    Given grid information for a patch, calculate the three *cell-centred*
    coordinate arrays in the *coord. system of the patch*.
    To convert to Cartesian coords., see `plotting_grid`.

    """

    if np.any(p.n != p.nwrite):
        if truecoords:
            x1 = p.corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
            x2 = p.corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
            x3 = p.corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]
        else:
            x1 = p.llc[0] - p.nghost[0] * p.dx[0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
            x2 = p.llc[1] - p.nghost[1] * p.dx[1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
            x3 = p.llc[2] - p.nghost[2] * p.dx[2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]
    else:
        if truecoords:
            x1 = p.active_corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
            x2 = p.active_corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
            x3 = p.active_corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]
        else:
            x1 = p.llc[0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
            x2 = p.llc[1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
            x3 = p.llc[2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]

    if p.ioformat == 2:
        x1 = x1 - p.dx[0] * 0.5
        x2 = x2 - p.dx[1] * 0.5
        x3 = x3 - p.dx[2] * 0.5

    return (x1, x2, x3)

def plotting_grid(p, x_i, slicedir, sliceloc, sliceeps=1.0e-2, truecoords=False,
                  phi_rotate=0.0, theta_rotate=0.0):
    """
    Given slice information, return the appropriate coord. vectors for plotting,
    including the conversion to Cartesian coords. (required for plotting with
    matplotlib).
    
    If the `truecoords` keyword is set True, then the conversion to Cartesian
    coords. is skipped.
    
    The following curvilinear coord. systems are known:
      (x,y,z), (Z,R,Phi), (r,theta,phi)
    The `phi_rotate` and `theta_rotate` arguments allow for patch coords. which
    are rotated relative to the convention of theta = phi = 0. They are silently
    passed to the coord. transformation functions.

    """

    assert (slicedir > 1) or (slicedir < 3), "slice direction must be 1, 2, or 3."

    if p.mesh_type == 'Cartesian' or truecoords:
        x, y, z = np.meshgrid(x_i[0], x_i[1], x_i[2], indexing='ij')
    elif p.mesh_type == 'cylindrical':
        x, y, z = zrp_to_xyz(x_i[0], x_i[1], x_i[2], phi_rotate=phi_rotate)
    elif p.mesh_type == 'spherical':
        x, y, z = rtp_to_xyz(x_i[0], x_i[1], x_i[2], phi_rotate=phi_rotate, theta_rotate=theta_rotate)

    # adjust coords. using the patch position (which is in Cartesian coords.).
    #x += p.pos[0]
    #y += p.pos[1]
    #z += p.pos[2]

    # assign coords. actually needed for plotting; cyclical.
    if slicedir == 1:
        xp, yp = y,z
    elif slicedir == 2:
        xp, yp = z,x
    elif slicedir == 3:
        xp, yp = x,y

    # account for symmetry (i.e., less than 3-D).
    sliceloc_act = sliceloc
    if p.mesh_type == 'Cartesian':
        if slicedir == 1 and len(x_i[0]) == 1: sliceloc_act = x_i[0]
        if slicedir == 2 and len(x_i[1]) == 1: sliceloc_act = x_i[1]
        if slicedir == 3 and len(x_i[2]) == 1: sliceloc_act = x_i[2]
    if p.mesh_type == 'cylindrical':
        if slicedir == 3 and len(x_i[0]) == 1: sliceloc_act = x_i[0]
    if p.mesh_type == 'spherical':
        if slicedir == 3 and len(x_i[1]) == 1: sliceloc_act = np.cos(x_i[1])

    # finally, determine the (nearest) slice location for the current patch.
    if slicedir == 1:
        diff = x - sliceloc_act
    elif slicedir == 2:
        diff = y - sliceloc_act
    elif slicedir == 3:
        diff = z - sliceloc_act
    # create mask
    m = np.ma.masked_outside(diff,-sliceeps,sliceeps)
    xpcut = np.ma.array(xp,mask=m.mask)
    ypcut = np.ma.array(yp,mask=m.mask)

    return xpcut, ypcut, m.mask

def patch_grid_edges(p, slicedir, scaling=1.0, truecoords=True):
    """
    Using the corners of the patch, create a matplotlib patch that draws
    the grid edges.

    """

    from matplotlib.patches import Rectangle, Wedge

    assert isinstance(p, Patch), "Input must be an instance of a DISPATCH `Patch`!"

    if truecoords:
        if slicedir == 1:
            llx, urx = p.active_corners[0][1], p.active_corners[1][1]
            lly, ury = p.active_corners[0][2], p.active_corners[1][2]
        elif slicedir == 2:
            llx, urx = p.active_corners[0][0], p.active_corners[1][0]
            lly, ury = p.active_corners[0][2], p.active_corners[1][2]
        elif slicedir == 3:
            llx, urx = p.active_corners[0][0], p.active_corners[1][0]
            lly, ury = p.active_corners[0][1], p.active_corners[1][1]
    else:
        if slicedir == 1:
            llx, urx = p.llc[1], p.llc[1] + p.size[1]
            lly, ury = p.llc[2], p.llc[2] + p.size[2]
        elif slicedir == 2:
            llx, urx = p.llc[0], p.llc[0] + p.size[0]
            lly, ury = p.llc[2], p.llc[2] + p.size[2]
        elif slicedir == 3:
            llx, urx = p.llc[0], p.llc[0] + p.size[0]
            lly, ury = p.llc[1], p.llc[1] + p.size[1]

    if p.mesh_type == mesh_types[1] or truecoords: # Cartesian
        llx *= scaling
        lly *= scaling
        urx *= scaling
        ury *= scaling
        edges = Rectangle((llx,lly),width=(urx-llx),height=(ury-lly),fc='none',lw=1.0,
                          ec='k',linestyle='dashed',zorder=0.01/min(p.dx)+0.1)
    elif p.mesh_type == mesh_types[2]: # spherical
        thetamin = round(p.active_corners[0][2]*180.0/np.pi,3)
        thetamax = round(p.active_corners[1][2]*180.0/np.pi,3)
        llx, urx = p.active_corners[0][0], p.active_corners[1][0]
        llx *= scaling
        lly *= scaling
        urx *= scaling
        ury *= scaling
        edges = Wedge((p.pos[0]*scaling,p.pos[1]*scaling),urx,thetamin,thetamax,width=(urx-llx),fc='none',
                      lw=1.0,ec='k',linestyle='dotted',zorder=0.01/min(p.dx)+0.1)
    elif p.mesh_type == mesh_types[3]: # cylindrical
        phimin = round(p.active_corners[0][2]*180.0/np.pi,3)
        phimax = round(p.active_corners[1][2]*180.0/np.pi,3)
        llx, urx = p.active_corners[0][1], p.active_corners[1][1]
        llx *= scaling
        lly *= scaling
        urx *= scaling
        ury *= scaling
        edges = Wedge((p.pos[0]*scaling,p.pos[1]*scaling),urx,phimin,phimax,width=(urx-llx),fc='none',
                      lw=1.0,ec='k',linestyle='dotted',zorder=0.01/min(p.dx)+0.1)

    return edges
