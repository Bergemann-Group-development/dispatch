# -*- coding: utf-8 -*-

"""
   Support functions for reading or calculating DISPATCH2 grid and geometry information.

"""

import numpy as np

class GeometricFactors(dict):
    """Calculate and store the geometric factors used by curvilinear grids."""

    def __init__(self, patch):
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
        # initialize the grid
        self.init_grid(patch)

    def init_grid(self, p):
        """Initialise geometric factors based on coord. type."""

        if p.mesh_type == 'Cartesian': self.init_Cartesian(p)
        elif p.mesh_type == 'cylindrical': self.init_cylindrical(p)
        elif p.mesh_type == 'spherical': self.init_spherical(p)

    def init_Cartesian(self, p):
        """Initialise geometric factors for a Cartesian coord. system."""

        n1, n2, n3 = p.ncell

        # 1-direction
        self['h2c'] = np.ones(n1)
        self['h2f'] = np.ones(n1)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        self['h32c'] = np.ones(n2)
        self['h32f'] = self['h32c'].view()

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.ds[0]
        self['dx1f'] = np.ones(n1) * p.ds[0]
        self['dx2c'] = np.ones(n2) * p.ds[1]
        self['dx2f'] = np.ones(n2) * p.ds[1]
        self['dx3c'] = np.ones(n3) * p.ds[2]
        self['dx3f'] = np.ones(n3) * p.ds[2]

        # volume elements
        self['dvol1c'] = np.ones(n1) * p.ds[0]
        self['dvol1f'] = np.ones(n1) * p.ds[0]
        self['dvol2c'] = np.ones(n2) * p.ds[1]
        self['dvol2f'] = np.ones(n2) * p.ds[1]
        self['dvol3c'] = np.ones(n3) * p.ds[2]
        self['dvol3f'] = np.ones(n3) * p.ds[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.ds[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.ds[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.ds[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.ds[0] / self['dvol1f']
        self['dar32c'] = p.ds[1] / self['dvol2c']
        self['dar32f'] = p.ds[1] / self['dvol2f']

    def init_cylindrical(self, p):
        """Initialise geometric factors for a cylindrical coord. system."""

        n1, n2, n3 = p.ncell

        # 1-direction
        self['h2c'] = np.ones(n1)
        self['h2f'] = np.ones(n1)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        pos_c = np.array(p.y )
        pos_f = np.array(p.ys)
        self['h32c'] = abs(pos_c)
        self['h32f'] = abs(pos_f)

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.ds[0]
        self['dx1f'] = np.ones(n1) * p.ds[0]
        self['dx2c'] = np.ones(n2) * p.ds[1]
        self['dx2f'] = np.ones(n2) * p.ds[1]
        self['dx3c'] = np.ones(n3) * p.ds[2]
        self['dx3f'] = np.ones(n3) * p.ds[2]

        # volume elements
        self['dvol1c'] = np.ones(n1) * p.ds[0]
        self['dvol1f'] = np.ones(n1) * p.ds[0]
        self['dvol2c'] = np.empty_like(pos_c)
        self['dvol2f'] = np.empty_like(pos_f)
        for j in xrange(len(pos_c)-1):
            self['dvol2c'][j] = 0.5 * abs( self['h32f'][j+1] * pos_f[j+1]
                                         - self['h32f'][j  ] * pos_f[j  ] )
            self['dvol2f'][j] = 0.5 * abs( self['h32c'][j  ] * pos_c[j  ]
                                         - self['h32c'][j-1] * pos_c[j-1] )
        self['dvol3c'] = np.ones(n3) * p.ds[2]
        self['dvol3f'] = np.ones(n3) * p.ds[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.ds[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.ds[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.ds[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.ds[0] / self['dvol1f']
        self['dar32c'] = p.ds[1] / self['dvol2c']
        self['dar32f'] = p.ds[1] / self['dvol2f']

    def init_spherical(self, p):
        """Initialise geometric factors for a spherical coord. system."""

        n1, n2, n3 = p.ncell

        # 1-direction
        rpos_c = np.array(p.x)
        rpos_f = np.array(p.xs)
        self['h2c'] = abs(rpos_c)
        self['h2f'] = abs(rpos_f)
        self['h31c'] = self['h2c'].view()
        self['h31f'] = self['h2f'].view()

        # 2-direction
        tpos_c = np.array(p.y)
        tpos_f = np.array(p.ys)
        self['h32c'] = abs(np.sin(tpos_c))
        self['h32f'] = abs(np.sin(tpos_f))

        # linear size elements
        self['dx1c'] = np.ones(n1) * p.ds[0]
        self['dx1f'] = np.ones(n1) * p.ds[0]
        self['dx2c'] = np.ones(n2) * p.ds[1]
        self['dx2f'] = np.ones(n2) * p.ds[1]
        self['dx3c'] = np.ones(n3) * p.ds[2]
        self['dx3f'] = np.ones(n3) * p.ds[2]

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
        self['dvol3c'] = np.ones(n3) * p.ds[2]
        self['dvol3f'] = np.ones(n3) * p.ds[2]

        # area elements
        self['dar1c'] = self['h31c'] * self['h2c']
        self['dar1f'] = self['h31f'] * self['h2f']
        self['dar2c'] = self['h31f'] * p.ds[0] / self['dvol1c']
        self['dar2f'] = self['h31c'] * p.ds[0] / self['dvol1f']
        self['dar31c'] = self['h2f'] * p.ds[0] / self['dvol1c']
        self['dar31f'] = self['h2c'] * p.ds[0] / self['dvol1f']
        self['dar32c'] = p.ds[1] / self['dvol2c']
        self['dar32f'] = p.ds[1] / self['dvol2f']
