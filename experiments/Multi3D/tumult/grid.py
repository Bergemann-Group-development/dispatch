import os
import numpy as np
import pandas as pd
from scipy import integrate, stats, optimize
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from .tumult import tumult
from .fast_trans import fast_trans
import copy
from scipy.spatial import Delaunay, ConvexHull, Voronoi
import re
from .CoG import CoG
from .auxfuncs import *
import multiprocessing as mp
from tqdm import tqdm
import plotly.graph_objs as go
from functools import partial

#=====================================================================
def read_run(folder):
  return tumult.read(folder, read_atmos=False, read_atom=False)

#=====================================================================
def read_list(folders):
  with mp.Pool() as pool:
    run_list = list(tqdm(pool.imap(read_run, folders), total=len(folders)))
    return np.array(run_list)

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
    self.points[:,0] /= 1e3
    self.delaunay = Delaunay(self.points)

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
  def get_CoG(self, star, method='Delaunay'):
    if method=='Delaunay':
      weights, cogs = self.get_weights(star, verbose=True)
    elif method=='NaturalNeighbor':
      weights, nbors = compute_nn_weights(self.points, star)
      cogs = self.cogs[nbors]

    if cogs is None: return None
    if len(cogs) == 1: return cogs[0]
    if weights is None: return None
    nabnds = np.array([len(cog.runs) for cog in cogs]) # number of abnds in each CoG
    if not np.all(nabnds == nabnds[0]):
      print('Warning: Not all CoGs have the same number of runs.')
    nabnds = nabnds[0]

    runs = []
    for i in range(nabnds):
      runs.append(tumult.combine([cog.runs[i] for cog in cogs], weights=weights))

    abnds = np.dot(weights, [cog.abnds for cog in cogs])

    return CoG(runs, abnds=abnds)
  
  #=====================================================================
  def plot_nodes(self, fig=None, size=5, color='black', cmap='Viridis', opacity=0.8, **kwargs):
    teff = self.points[:,0]
    logg = self.points[:,1]
    feh = self.points[:,2]

    if fig is None:
      fig = go.Figure()

    # Create a 3D scatter plot
    fig.add_trace(go.Scatter3d(
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
    )

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

    return fig
  
  #=====================================================================
  def plot_point(self, teff, logg, feh, fig=None, **kwargs):
    if fig is None:
      fig = go.Figure()

    fig.add_trace(go.Scatter3d(
      x=[teff],
      y=[feh],
      z=[logg],
      **kwargs,
      ),
    )

    return fig


#=====================================================================
def is_point_in_array(point, points, tol=1e-9):
    point = np.asarray(point)
    points = np.asarray(points)
    # Ensure point shape matches array dimensions
    if point.shape != (points[0].shape):
        raise ValueError(f"Point must have shape ({points[0].shape},), but got {point.shape}")
    
    # Find matching rows
    matches = np.all(np.isclose(points, point, atol=tol), axis=1)
    if np.any(matches):
        return np.where(matches)[0][0]
    else:
        return -1

#=====================================================================
def get_nbors(vor, query_index):
    """ Find neighbors of the query point in the Voronoi diagram. """
    if query_index < 0: query_index += len(vor.points)
    mask = np.any(vor.ridge_points == query_index, axis=1)
    point_ridges = vor.ridge_points[mask]
    nbors = np.sum(point_ridges, axis=1) - query_index
    ridge_vertices = [vor.ridge_vertices[i] for i in np.where(mask)[0]]
    return nbors, ridge_vertices

#=====================================================================
def point_in_convex_hull(point, hull):
    # Use the hull equations to check if the point lies inside
    return np.all(np.dot(hull.equations[:, :-1], point) + hull.equations[:, -1] <= 0)

#=====================================================================  
def compute_enclosed_measure(points):
    points = np.asarray(points)
    dim = points.shape[1]  # Number of columns determines dimension
    if dim == 1:  # 1D: Compute length
        return np.max(points) - np.min(points)
    elif dim >= 2:  # 2D, 3D, or higher: Use ConvexHull
        hull = ConvexHull(points)
        return hull.volume if dim > 2 else hull.area  # Area for 2D, Volume for 3D+
    else:
        raise ValueError("Invalid input: The input points must have at least 1 dimension.")

#=====================================================================
def compute_nn_weights(points, query, draw=False):
    idx = is_point_in_array(query, points)
    if idx >= 0:
        return np.array([1]), np.array([idx])
    
    hull = ConvexHull(points)
    if not point_in_convex_hull(query, hull):
        return np.array([0]), np.array([-1])
    
    vor1 = Voronoi(points)
    vor2 = Voronoi(np.vstack([points, query]))
    query_region = vor2.regions[vor2.point_region[-1]]
    query_hull = ConvexHull(vor2.vertices[query_region])
    nbors, ridge_verts = get_nbors(vor2, -1)
    areas = np.zeros(len(nbors))

    if draw:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for n, nbor in enumerate(nbors):
        region = vor1.regions[vor1.point_region[nbor]]
        vertices = vor2.vertices[ridge_verts[n]]
        removed_vertices = []
        for i in region:
            if i >= 0:
                point = vor1.vertices[i]
                if point_in_convex_hull(point, query_hull):
                    removed_vertices.append(point)

        
        vertices = np.vstack([vertices, removed_vertices])
        areas[n] = compute_enclosed_measure(vertices)
        if draw:
            hull = ConvexHull(vertices)
            vertices = vertices[hull.vertices]
            plt.gca().add_patch(Polygon(vertices, color=colors[n], alpha=0.5))

    weights = areas / np.sum(areas)
    return weights, nbors


#=====================================================================
def compute_delaunay_weights(points, query, delaunay=None):
    idx = is_point_in_array(query, points)
    if idx >= 0:
        return np.array([1]), np.array([idx])
    
    hull = ConvexHull(points)
    if not point_in_convex_hull(query, hull):
        return np.array([0]), np.array([-1])
    
    if delaunay is None:
      delaunay = Delaunay(points)

    simplex = delaunay.find_simplex(query)

    if simplex == -1:
        return np.array([0]), -1

    vertices = delaunay.simplices[simplex]
    simplex_points = points[vertices]
    
    A = np.hstack((simplex_points, np.ones((len(simplex_points), 1))))
    b = np.append(query, 1)
    
    weights = np.linalg.solve(A.T, b)
    mask = weights > 1e-10
    weights = weights[mask]
    
    return weights, vertices[mask]

def keep_k_closest(points, query, k=100):
    """Keep only the k closest points to query_point"""
    # Calculate distances
    distances = np.linalg.norm(points - query, axis=1)
    
    # Get indices of k smallest distances (O(n) operation)
    k_closest_indices = np.argpartition(distances, k)[:k]
    
    return points[k_closest_indices], k_closest_indices