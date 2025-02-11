import os
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy import integrate, stats, optimize
import matplotlib.pyplot as plt
from .m3dis import m3dis
from .auxfuncs import plt_kwargs
from .fast_trans import fast_trans
import copy
from scipy.spatial import Delaunay
import re
from .CoG2 import CoG
from .auxfuncs import *
import multiprocessing as mp
from tqdm import tqdm
import plotly.graph_objs as go

def read_run(folder):
  return m3dis.read(folder, read_atmos=False, read_atom=False)

#=====================================================================
class grid(CoG):
  def __init__(self, runs, params=None):
    if type(runs[0]) == str:
      with mp.Pool() as pool:
        runs = list(tqdm(pool.imap(read_run, runs), total=len(runs)))

    cogs = []
    points = []
    atmos_files = np.unique([run.nml['atmos_params']['atmos_file'] for run in runs])
    for i, atmos_file in enumerate(atmos_files):
      cog = CoG([run for run in runs if run.nml['atmos_params']['atmos_file'] == atmos_file])

      if params is None:
        sname = re.findall(r't\d{2,4}g\d{2}[mp]\d{2}', atmos_file)
        mname = re.findall(r'[sp]\d{4}.+g.+m.+t.+z.+\.mod', atmos_file)

        if len(sname) > 0:
          cog.star = stagger_name_to_params(sname[0])
        elif len(mname) > 0:
          cog.star = marcs_name_to_params(mname[0])
        else:
          print(atmos_file)
          raise ValueError('Atmos name could not be converted into stellar parameters')
      else:
        cog.star = params[i]

      cog.label = 'Teff={:.0f}, logg={:.1f}, feh={:+.1f}'.format(*cog.star[:3])
      if len(cog.star) > 3:
        cog.label += ', vmic={:.1f}'.format(cog.star[3])

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
  
  #=====================================================================
  def get_CoG(self, star):
    weights, cogs = self.get_weights(star, verbose=True)
    if cogs is None: return None
    if len(cogs) == 1: return cogs[0]
    if weights is None: return None
    nabnds = np.array([len(cog.runs) for cog in cogs]) # number of abnds in each CoG
    if not np.all(nabnds == nabnds[0]):
      print('Warning: Not all CoGs have the same number of runs.')
    nabnds = nabnds[0]

    runs = []
    for i in range(nabnds):
      runs.append(m3dis.combine([cog.runs[i] for cog in cogs], weights=weights))

    abnds = np.dot(weights, [cog.abnds for cog in cogs])

    return CoG(runs, abnds=abnds)
  
  #=====================================================================
  def plot_nodes(self, size=5, color='k', cmap='Viridis', opacity=0.8, **kwargs):
    teff = self.points[:,0]
    logg = self.points[:,1]
    feh = self.points[:,2]

    # Create a 3D scatter plot
    fig = go.Figure(data=[
      go.Scatter3d(
        x=teff,
        y=feh,
        z=logg,
        mode='markers',
        marker=dict(
            size=size,
            color=color,  # Set color to the 'feh' values
            colorscale=cmap,  # Choose a colorscale
            opacity=opacity,
            **kwargs
        )),
    ])

    # Set labels
    fig.update_layout(
        scene=dict(
            xaxis_title='Teff',
            yaxis_title='Fe/H',
            zaxis=dict(
                title='Logg',
                autorange='reversed'  # Invert the z-axis
            )
        ),
        margin=dict(l=0, r=0, b=0, t=0)  # Reduce padding around the plot
    )

    # Show plot
    fig.show()
