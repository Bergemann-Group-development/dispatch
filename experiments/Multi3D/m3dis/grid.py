import os
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy import integrate, stats, optimize
import matplotlib.pyplot as plt
from .auxfuncs import plt_kwargs
from .fast_trans import fast_trans
import copy
from scipy.spatial import Delaunay
import re
from .CoG2 import CoG
from .auxfuncs import *

#=====================================================================
class grid(CoG):
  def __init__(self, runs):
    cogs = []
    points = []
    atmos_files = np.unique([run.nml['atmos_params']['atmos_file'] for run in runs])
    for atmos_file in atmos_files:
      cog = CoG([run for run in runs if run.nml['atmos_params']['atmos_file'] == atmos_file])

      sname = re.findall(r't\d{2,4}g\d{2}[mp]\d{2}', atmos_file)
      mname = re.findall(r'[sp]\d{4}.+g.+m.+t.+z.+\.mod', atmos_file)

      if len(sname) > 0:
        cog.star = stagger_name_to_params(sname[0])
      elif len(mname) > 0:
        cog.star = marcs_name_to_params(mname[0])
      else:
        print(atmos_file)
        raise ValueError('Atmos name could not be converted into stellar parameters')
    
      cog.label = 'Teff={:.0f}, logg={:.1f}, feh={:+.1f}'.format(*cog.star)
      points.append(cog.star)
      cogs.append(cog)

    self.runs = np.array(runs)
    self.cogs = np.array(cogs)
    self.points = np.array(points)
    self.delaunay = Delaunay(points)

  #=====================================================================
  def get_weights(self, query_point, verbose=True):
    # Find the simplex containing the query point
    simplex = self.delaunay.find_simplex(query_point)

    if simplex == -1:
        if verbose: print("Query point is outside the convex hull.")
        return None, None

    # Get the vertices of the simplex
    vertices = self.delaunay.simplices[simplex]

    # Compute barycentric coordinates
    simplex_points = self.points[vertices]  # Coordinates of the simplex vertices
    A = np.hstack((simplex_points, np.ones((len(simplex_points), 1))))  # Add 1 for affine coordinates
    b = np.append(query_point, 1)  # Query point in affine coordinates
    weights = np.linalg.solve(A.T, b)  # Solve for barycentric coordinates
    mask = weights > 1e-10  # Mask for vertices with non-zero weight
    weights = weights[mask]
    cogs = self.cogs[vertices[mask]]
    return weights, cogs

  #=====================================================================
  def crop(self, star, **kwargs):
    weights, runs = self.get_weights(star)
    xy = np.array([run.crop(**kwargs) for run in runs])
    xx = xy[0, 0]
    yy = np.sum(xy[:,1].T * weights, axis=1)
    return xx, yy