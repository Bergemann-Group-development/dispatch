# Pythn 2/3 compatibility
from __future__ import print_function
import sys
if sys.version_info[0]!=2:
    from importlib import reload

import numpy as np

"""
   Support functions for reading legacy format DISPATCH2 data.
"""

mesh_types = {1: 'Cartesian',
              2: 'spherical',
              3: 'cylindrical',}
mesh_types_r = {v: k for k, v in mesh_types.items()}

def get_mesh_type(current_mesh_type):
    """
    A quick and hard-wired way to turn the integer mesh type index into a
    string. The following must mirror the `mesh_types` variable in the 
    header of `mesh_mod.f90`.
    """
    assert isinstance(current_mesh_type, int), "expecting an integer for `mesh_type`"
    try:
        mesh_string = mesh_types[current_mesh_type]
    except:
        raise Exception("Unrecognised mesh type!")
    return mesh_string
              
def info (self):
    print("        patch id:", self.id)
    print("    no. of steps:", self.istep)
    print("         coords.:", self.mesh_type)
    print("      full patch:", self.m)
    print("   unigrid patch:", self.n)
    print("     total cells:", self.gn)
    print("   lower indices:", self.li)
    print("   upper indices:", self.ui)
    print("no. of variables:", self.nvar)
    print("     ghost cells:", self.nghost)
    print("            time:", self.time)
    print("        position:", self.pos)
    print("         delta s:", self.ds)
    print("        velocity:", self.vel)
    print("          origin:", self.origin)
    print("          l.l.c.:", self.llc)
    print("     vol. centre:", self.centre)
    print("           level:", self.level)
    print("         quality:", self.quality)
    print("       task type:", self.kind)

def read_version0 (self, fp, line, verbose):
    """
    An concrete method for reading information for a single patch.
    The following must mirror what is in `patch_mod.f90::output_legacy`.
    """
    self.kind    = 'ramses'
    self.eos     = 'ideal'
    self.opacity = 'none'
    self.gamma   = 1.4
    # line 1: time, pos(1:3), ds(1:3), vel(1:3), origin(1:3), llc(1:3), cent(1:3)
    self.time    = float(line[0])
    self.level   = int  (line[19])
    self.quality = float(line[20])
    self.pos     = np.array(line[ 1: 4],dtype=np.float64)
    self.ds      = np.array(line[ 4: 7],dtype=np.float64)
    self.vel     = np.array(line[ 7:10],dtype=np.float64)
    self.origin  = np.array(line[10:13],dtype=np.float64)
    self.llc     = np.array(line[13:16],dtype=np.float64)
    self.centre  = np.array(line[16:19],dtype=np.float64)
    # line 2: m(1:3), n(1:3), li(1:3), ui(1:3), gn(1:3), nghost(1:3), nvar
    line = fp.readline().split()
    self.m      = np.array(line[ 0: 3],dtype=np.int32)
    self.nwrite =  (int(line[3]),int(line[4]), int(line[5]))
    self.li     = np.array(line[ 6: 9],dtype=np.int32)-1
    self.ui     = np.array(line[ 9:12],dtype=np.int32)-1
    self.n      = np.array(line[12:15],dtype=np.int32)
    self.gn     = np.array(line[15:18],dtype=np.int32)
    self.nghost = np.array(line[18:21],dtype=np.int32)
    self.nvar   = int(line[21])
    # line 3: id, iout, istep, mesh_type
    line = fp.readline().split()
    self.id    = int(line[0])
    self.iout  = int(line[1])
    self.istep = int(line[2])
    self.dt    = float(line[3])
    self.mesh_type = get_mesh_type(int(line[4]))
    if (verbose):
        info (self)
  
def read_version1 (self, fp, verbose):
    """
    An concrete method for reading information for a single patch.
    The following must mirror what is in `patch_mod.f90::output_legacy`.
    """
    # lines 1-3: solver, eos, opacity
    self.kind    = fp.readline().strip()
    self.eos     = fp.readline().strip()
    self.opacity = fp.readline().strip()
    # line 4: time, dtime, pos(1:3), ds(1:3), vel(1:3), origin(1:3), llc(1:3), cent(1:3)
    line = fp.readline().split()
    self.time,self.dt       = np.array(line[ 0: 2],dtype=np.float64)
    self.pos                = np.array(line[ 2: 5],dtype=np.float64)
    self.ds                 = np.array(line[ 5: 8],dtype=np.float64)
    self.vel                = np.array(line[ 8:11],dtype=np.float64)
    self.origin             = np.array(line[11:14],dtype=np.float64)
    self.llc                = np.array(line[14:17],dtype=np.float64)
    self.centre             = np.array(line[17:20],dtype=np.float64)
    self.quality,self.gamma = np.array(line[20:22],dtype=np.float64)
    # line 5: m(1:3), n(1:3), li(1:3), ui(1:3), gn(1:3), nghost(1:3), nvar
    line = fp.readline().split()
    self.m      = np.array(line[ 0: 3],dtype=np.int32)
    self.nwrite =  (int(line[3]),int(line[4]), int(line[5]))
    self.li     = np.array(line[ 6: 9],dtype=np.int32)-1
    self.ui     = np.array(line[ 9:12],dtype=np.int32)-1
    self.n      = np.array(line[12:15],dtype=np.int32)
    self.gn     = np.array(line[15:18],dtype=np.int32)
    self.nghost = np.array(line[18:21],dtype=np.int32)
    self.nvar   = np.int  (line[21   ])
    self.level  = np.int  (line[22   ])
    # line 3: id, iout, istep, mesh_type
    line = fp.readline().split()
    self.id    = np.int(line[0])
    self.iout  = np.int(line[1])
    self.istep = np.int(line[2])
    self.mesh_type = get_mesh_type(int(line[3]))

def single_patch (self, fp, fd=None, verbose=False):
    if fd:
        self.ioformat=3
    else:
        line=fp.readline()
        chars=line[0:2]
        if chars=='  ':
            self.ioformat=0
            line=line.split()
            read_version0 (self, fp, line, verbose)
        else:
            self.ioformat=np.int(chars)
            read_version1 (self, fp, verbose)
    derived (self, verbose)

def derived(self, verbose=False):
    self.dx   = self.ds
    self.corners        = self._determine_corners()
    self.active_corners = self._determine_corners(active_only=True)
    self.corners = np.array(self.corners)
    self.active_corners = np.array(self.active_corners)
    self.size = self.ds*self.n
    self.x=self.corners[0][0]+(0.5+np.arange(self.gn[0]))*self.ds[0]
    self.y=self.corners[0][1]+(0.5+np.arange(self.gn[1]))*self.ds[1]
    self.z=self.corners[0][2]+(0.5+np.arange(self.gn[2]))*self.ds[2]
    if self.ioformat==2 or self.ioformat==4:
        self.size = self.ds*(self.n-1)
        self.active_corners[1] = self.active_corners[1] - self.ds
        self.corners[1]        = self.corners[1]        - self.ds
        self.x = self.x - 0.5*self.ds[0]
        self.y = self.y - 0.5*self.ds[1]
        self.z = self.z - 0.5*self.ds[2]
    self.xs = self.x - 0.5*self.ds[0]
    self.ys = self.y - 0.5*self.ds[1]
    self.zs = self.z - 0.5*self.ds[2]
    a = self.active_corners
    self.extent=np.array([[a[0,1],a[1,1],a[0,2],a[1,2]],\
                          [a[0,0],a[1,0],a[0,2],a[1,2]],\
                          [a[0,0],a[1,0],a[0,1],a[1,1]]])
    if (verbose):
        info (self)
