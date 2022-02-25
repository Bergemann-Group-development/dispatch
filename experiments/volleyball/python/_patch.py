import copy as copy
import numpy as np
import pylab

class patch():
    def __init__(self,point,verbose=0,pos=np.array([0.5,0.5,0.5]),size=np.array([1.0,1.0,1.0]),n_cell=20,ng=2,mem=False):
        self.verbose = verbose
        self.position= np.zeros(3); self.position[:]= point.xyz()
        self.size    = np.zeros(3); self.size[:]    = point.size
        self.point   = copy.deepcopy(point)
        self.point.unit_vectors()
        self.n       = np.zeros(3); self.n[:]       = n_cell
        self.ng      = np.zeros(3); self.ng[:]      = ng
        self.n_cell  = n_cell
        self.new_mesh()
        if mem:
          self.add_mem()

    def new_mesh(self):
        self.lb  = np.zeros(3,dtype='int'); self.lb[:]  = 0
        self.li  = np.zeros(3,dtype='int'); self.li[:]  = self.lb + self.ng
        self.lo  = np.zeros(3,dtype='int'); self.lo[:]  = self.li - 1
        self.ui  = np.zeros(3,dtype='int'); self.ui[:]  = self.li + self.n - 1
        self.uo  = np.zeros(3,dtype='int'); self.uo[:]  = self.ui + 1
        self.ub = np.zeros(3, dtype='int')
        self.ub[:] = self.ui + self.ng
        self.gn  = self.ub
        self.r   = np.zeros((3,self.gn[0]+1))
        self.o   = (self.li+self.ui)*0.5
        self.d   = np.zeros(3)
        self.d   = self.size/self.n
        self.llc = self.position - 0.5*self.size
        self.urc = self.position + 0.5*self.size
        self.centre = self.llc + 0.5*self.size
        for i in range(self.gn[0]+1):
            self.r[:,i]=(i+1-0.5-self.ng)*self.d+self.llc-self.centre

    def add_mem(self):
      """ add pseudo memory array. For vectors trasnformation testing. """
      self.mem=np.zeros((self.gn[0]+1,self.gn[1]+1,self.gn[2]+1,3))

    def cell2pos(self,p):
      """ convert cell index into position in the local reference frame """
      return [self.r[0,p[0]],self.r[1,p[1]],self.r[2,p[2]]]


    def pos2cell(self,p):
      """ convert position in the local reference frame into cell index """
      pos = p/self.d + self.o + 0.0001
      return pos.astype(int)

    def loc2glob(self,p):
      """ convert local position into global cartesian coordinates """
      return np.matmul(self.point.e, p)+self.position

    def glob2loc(self, p):
      """ convert global cartesian coordinates into local position """
      return np.matmul(p-self.position, self.point.e)

    def cell2pos_glob(self,p):
      """ return cell position in global cartesian coordinates """
      pos = [self.r[0, p[0]], self.r[1, p[1]], self.r[2, p[2]]]
      return self.loc2glob(pos)

    def to_my_coords(self, nbor, pos):
      """ convert point position in nbor frame of reference to local coordinates """
      dp = nbor.position-self.position
      npos = np.matmul(np.matmul(nbor.point.e, pos), self.point.e)+np.matmul(dp, self.point.e)
      return npos


    def to_my_coords2(self, nbor, pos):
      """ convert point position in nbor frame of reference to local coordinates """
      e1=self.point.e
      e2=nbor.point.e.transpose()
      dp = nbor.position-self.position
      e=np.matmul(e2,e1)
      npos = np.matmul(pos, e)+np.matmul(dp, e1)
      return npos

    def indices(self, tl, tu, target, d, e):
      """ transformation from target index space to source index space. reuses some
          things from to_my_coords, so probably this can be optimized."""
      ps=np.zeros((self.gn[0]+1,self.gn[1]+1,self.gn[2]+1,3))
      c=1/self.d
      o=self.o
      for i in range(tl[0],tu[0]+1):
        for j in range(tl[1], tu[1]+1):
          for k in range(tl[2], tu[2]+1):
            pos=[target.r[0,i],target.r[1,j],target.r[2,k]]
            ps[i,j,k,:] = (np.matmul(e,pos)-d)*c+o
      return ps




    def erot(self, nbor):
      """ Return rotation matrix between two patches """
      e1 = self.point.e
      e2 = nbor.point.e.transpose()
      e = np.matmul(e2, e1)
      return e

    def get_position(self,source):
      """ get source position and orientation in target coordinates """
      dp = source.position-self.position
      pos=np.matmul(dp,self.point.e)
      e  = self.erot(source)
      return pos, e

    def all_corners(self,guards=False):
      if guards:
          l = self.lb
          u = self.ub
      else:
          l = self.li
          u = self.ui
      pl = [self.r[0, l[0]], self.r[1, l[1]], self.r[2, l[2]]]-0.5*self.d
      pu = [self.r[0, u[0]], self.r[1, u[1]], self.r[2, u[2]]]+0.5*self.d

      p1 = [pl[0], pl[1], pl[2]]
      p2 = [pu[0], pl[1], pl[2]]
      p3 = [pl[0], pu[1], pl[2]]
      p4 = [pu[0], pu[1], pl[2]]
      p5 = [pl[0], pl[1], pu[2]]
      p6 = [pu[0], pl[1], pu[2]]
      p7 = [pl[0], pu[1], pu[2]]
      p8 = [pu[2], pu[2], pu[2]]
      return [p1,p2,p3,p4,p5,p6,p7,p8]

    def from_source(self, lt, ut, src, filled=True):
          # convert to per unit mass here
          # get the distance and rotation matrix from src to self
      d, e = self.get_position(src)
      #Interpolation coordinates and weights p[:,:,:,0:2]
      p = src.indices(lt, ut, self,d, e)
      l = src.li
      u = src.ui
      m = np.shape(self.mem)
      nv = m[3]
      tp = np.zeros((m[0], m[1], m[2], m[3]))
      sp = src.mem

      for i in range(lt[0], ut[0]+1):
          for j in range(lt[1], ut[1]+1):
              for k in range(lt[2], ut[2]+1):
                  for iv in range(nv):
                      j1 = p[i, j, k, 0]
                      j2 = p[i, j, k, 1]
                      j3 = p[i, j, k, 2]
                      j1 = min(max(j1, l[0]), u[0]-1)
                      j2 = min(max(j2, l[1]), u[1]-1)
                      j3 = min(max(j3, l[2]), u[2]-1)
                      p1 = p[i, j, k, 0] - j1
                      p2 = p[i, j, k, 1] - j2
                      p3 = p[i, j, k, 2] - j3
                      q1 = 1-p1
                      q2 = 1-p2
                      q3 = 1-p3
                      tp[i, j, k, iv] = \
                          q3*(q2*(q1*sp[j1, j2, j3,     iv] + p1*sp[j1+1, j2, j3,    iv]) +\
                              p2*(q1*sp[j1, j2+1, j3,   iv] + p1*sp[j1+1, j2+1, j3,  iv])) + \
                          p3*(q2*(q1*sp[j1, j2, j3+1,   iv] + p1*sp[j1+1, j2, j3+1,  iv]) +\
                              p2*(q1*sp[j1, j2+1, j3+1, iv] + p1*sp[j1+1, j2+1, j3+1,iv]))

      outside = np.zeros(3, dtype='bool')
      for i in range(lt[0], ut[0]+1):
          outside[0] = i < self.li[0] or i > self.ui[0]
          for j in range(lt[1], ut[1]+1):
              outside[1] = j < self.li[1] or j > self.ui[1]
              for k in range(lt[2], ut[2]+1):
                  outside[2] = k < self.li[2] or k > self.ui[2]
                  if any(outside):
                      v = tp[i, j, k, :]
                      v = np.matmul(v,e)
                      self.mem[i, j, k, :] = v
                      self.filled[i,j,k]   = src.point.id



    def plot_patch(self, ax, guards=False, rotate=False, color='k', alpha=0.8, linestyle='solid'):
      [p1,p2,p3,p4,p5,p6,p7,p8] = self.all_corners(guards)
      # rotate the patch
      if rotate:
          e = self.point.e
          #p1=matmul(p1,e)
          #p2=matmul(p2,e)
          #p3=matmul(p3,e)
          #p4=matmul(p4,e)
          #p5=matmul(p5,e)
          #p6=matmul(p6,e)
          #p7=matmul(p7,e)
          #p8=matmul(p8,e)
          p1 = np.matmul(e, p1)
          p2 = np.matmul(e, p2)
          p3 = np.matmul(e, p3)
          p4 = np.matmul(e, p4)
          p5 = np.matmul(e, p5)
          p6 = np.matmul(e, p6)
          p7 = np.matmul(e, p7)
          p8 = np.matmul(e, p8)
      p1 = p1+self.position
      p2 = p2+self.position
      p3 = p3+self.position
      p4 = p4+self.position
      p5 = p5+self.position
      p6 = p6+self.position
      p7 = p7+self.position
      p8 = p8+self.position
      # bottom
      ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p2[0], p4[0]], [p2[1], p4[1]], [p2[2], p4[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p3[0], p4[0]], [p3[1], p4[1]], [p3[2], p4[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p1[0], p3[0]], [p1[1], p3[1]], [p1[2], p3[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      # top
      ax.plot([p5[0], p6[0]], [p5[1], p6[1]], [p5[2], p6[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p6[0], p8[0]], [p6[1], p8[1]], [p6[2], p8[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p7[0], p8[0]], [p7[1], p8[1]], [p7[2], p8[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p5[0], p7[0]], [p5[1], p7[1]], [p5[2], p7[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      # walls
      ax.plot([p1[0], p5[0]], [p1[1], p5[1]], [p1[2], p5[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p2[0], p6[0]], [p2[1], p6[1]], [p2[2], p6[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p3[0], p7[0]], [p3[1], p7[1]], [p3[2], p7[2]],
              color=color, alpha=alpha, linestyle=linestyle)
      ax.plot([p4[0], p8[0]], [p4[1], p8[1]], [p4[2], p8[2]],
              color=color, alpha=alpha, linestyle=linestyle)
