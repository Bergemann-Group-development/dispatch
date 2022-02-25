import copy as copy
import numpy as np
import pylab
import matplotlib.pyplot as pl
import _sphere_points as sp

class volleyball():
  def __init__(self,n_wide=12,verbose=0,r_min=0.7, r_max=0.75):
      self.verbose=0
      self.r_min=r_min
      self.r_max=r_max
      f=1-np.arctan(np.sqrt(0.5))/n_wide/2
      self.radius=r_max*f
      self.ntheta=n_wide
  #
  def symmetries(self,point, points):
      for i in (0,1):  # loop over front and back
          for j in (0,1,2): # loop over the three front faces: East, North, Up
              new=copy.deepcopy(point)
              point.id+=1
              points.append(new)
              point.permute()
          point.reflect()
  #
  def set_radius(self,radius):
      self.radius=radius
  #
  def generate_layer(self,symmetric=True):
      self.theta0=np.arctan(np.sqrt(0.5))
      self.dtheta=self.theta0/self.ntheta
      point=sp.point_on_sphere(self.dtheta)
      point.r = self.radius
      point.size= self.radius*self.dtheta/(1-self.dtheta)
      point.r = point.r - point.size/2
      save=copy.deepcopy(point)
      points=[]
      point.go_north(0.5*self.dtheta)
      point.go_east(0.5*self.dtheta)
      while(point.is_inside(self.theta0,self.dtheta)):
          #point.r = self.radius #- point.size/2
          save = copy.deepcopy(point)
          if symmetric:
              while (point.is_inside(self.theta0,self.dtheta)):
                  self.symmetries(point, points)
                  point.east_west()
                  self.symmetries(point, points)
                  point.north_south()
                  self.symmetries(point, points)
                  point.east_west()
                  self.symmetries(point, points)
                  point.north_south()
                  # back to start
                  point.go_east(self.dtheta)
          else:
              while (point.is_inside(self.theta0,self.dtheta)):
                  new=copy.deepcopy(point)
                  points.append(new)
                  point.east_west()
                  new=copy.deepcopy(point)
                  points.append(new)
                  point.north_south()
                  new=copy.deepcopy(point)
                  points.append(new)
                  point.east_west()
                  new=copy.deepcopy(point)
                  points.append(new)
                  point.north_south()
                  # back to start
                  point.go_east(self.dtheta)
          save.id=point.id
          point=save
          # next row
          point.go_north(self.dtheta)
          save=copy.deepcopy(point)
      sz=point.size
      step=sz*(1-self.dtheta)*(1-self.dtheta**2)
      new_rad=self.radius-step
      return points, new_rad

  def generate_nbors(self):
    self.theta0 = np.arctan(np.sqrt(0.5))
    self.dtheta = self.theta0/self.ntheta
    point = sp.point_on_sphere(self.dtheta)
    radius = self.radius
    # generate a set of nbors outside
    point.r = radius
    point.size = self.radius*self.dtheta/(1-self.dtheta)
    point.r = point.r - point.size/2
    save = copy.deepcopy(point)
    points = []
    point.go_north(0.5*self.dtheta)
    point.go_east(0.5*self.dtheta)
    new = copy.deepcopy(point)
    points.append(new)
    point.id+=1
    self.generate_around_point(point,points)
    save.id=point.id

    # generate a set of nbors on the same radius
    sz = point.size
    step = sz*(1-self.dtheta)*(1-self.dtheta**2)
    radius = radius-step
    point = copy.deepcopy(save)
    point.r = radius
    point.size = radius*self.dtheta/(1-self.dtheta)
    point.r = point.r - point.size/2
    save = copy.deepcopy(point)
    point.go_north(0.5*self.dtheta)
    point.go_east(0.5*self.dtheta)
    new = copy.deepcopy(point)
    points.insert(0,new)
    point.id+=1
    self.generate_around_point(point, points)
    save.id=point.id

    # generate a set of nbors inside
    sz = point.size
    step = sz*(1-self.dtheta)*(1-self.dtheta**2)
    radius = radius-step
    point = copy.deepcopy(save)
    point.r = radius
    point.size = radius*self.dtheta/(1-self.dtheta)
    point.r = point.r - point.size/2
    save = copy.deepcopy(point)
    point.go_north(0.5*self.dtheta)
    point.go_east(0.5*self.dtheta)
    new = copy.deepcopy(point)
    points.append(new)
    point.id+=1
    self.generate_around_point(point, points)

    return points

  def generate_pt_nbors(self,point):
    self.theta0 = np.arctan(np.sqrt(0.5))
    self.dtheta = self.theta0/self.ntheta
    radius = self.radius
    # generate a set of nbors outside
    save = copy.deepcopy(point)
    points = []
    new = copy.deepcopy(point)
    points.append(new)
    point.id += 1
    self.generate_around_point(point, points)
    save.id = point.id

    # generate a set of nbors on the same radius
    sz = point.size
    step = sz*(1-self.dtheta)*(1-self.dtheta**2)
    radius = radius-step
    point = copy.deepcopy(save)
    point.r = radius
    point.size = radius*self.dtheta/(1-self.dtheta)
    point.r = point.r - point.size/2
    save = copy.deepcopy(point)
    point.go_north(0.5*self.dtheta)
    point.go_east(0.5*self.dtheta)
    new = copy.deepcopy(point)
    points.insert(0, new)
    point.id += 1
    self.generate_around_point(point, points)
    save.id = point.id

    # generate a set of nbors inside
    sz = point.size
    step = sz*(1-self.dtheta)*(1-self.dtheta**2)
    radius = radius-step
    point = copy.deepcopy(save)
    point.r = radius
    point.size = radius*self.dtheta/(1-self.dtheta)
    point.r = point.r - point.size/2
    save = copy.deepcopy(point)
    point.go_north(0.5*self.dtheta)
    point.go_east(0.5*self.dtheta)
    new = copy.deepcopy(point)
    points.append(new)
    point.id += 1
    self.generate_around_point(point, points)

    return points

  def generate_around_point(self,point, points):
    new = copy.deepcopy(point)
    point.id+=1
    new.go_west(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_west(self.dtheta)
    new.go_north(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_west(self.dtheta)
    new.go_south(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_north(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_south(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_east(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_east(self.dtheta)
    new.go_north(self.dtheta)
    points.append(new)
    #
    new = copy.deepcopy(point)
    point.id += 1
    new.go_east(self.dtheta)
    new.go_south(self.dtheta)
    points.append(new)

  def generate_nbor(self, point):
    self.theta0 = np.arctan(np.sqrt(0.5))
    self.dtheta = self.theta0/self.ntheta
    new = copy.deepcopy(point)
    new.go_west(self.dtheta)
    new.id+=1
    return new
