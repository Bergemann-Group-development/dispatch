from .auxfuncs import pchip_mono
from .spectrum import SpectrumClass
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, ConvexHull, Voronoi
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

class SpectrumGrid(object):
  """
  A class to hold a grid of spectra.
  """
  def __init__(self, spectra, params, normalize=True, param_names=None):
    self.spectra = spectra
    self.params = np.array(params)
    self.npoints = len(self.spectra)
    if self.npoints != self.params.shape[0]:
      raise ValueError("Number of spectra and parameters must have the same first dimension.")
    
    if self.params.ndim == 1:
      self.ndim = 1
      self.param_names = [param_names] or ['param_0']
    else:
      self.ndim = self.params.shape[-1]
      self.param_names = param_names or [f'param_{i}' for i in range(self.ndim)]

    if normalize:
       def normalize_params(p):
          return (p - np.min(params, axis=0)) / np.ptp(params, axis=0)
       
       # establish normalization routine
       self.normalize_params = normalize_params
       self.points = self.normalize_params(self.params)
    else:
        self.points = self.params
    
    if self.points.ndim > 1 and self.points.shape[-1] > 1:
      self.delaunay = Delaunay(self.points)
    else:
      self.delaunay = None
      self.points = self.points.ravel()
      order = np.argsort(self.points)
      self.points = self.points[order]
      self.spectra = [self.spectra[i] for i in order]

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
    spectra = [self.spectra[i] for i in vertices[mask]]
    return weights, spectra
  
  #=====================================================================
  def crop(self, *args, method='Delaunay', **kwargs):
    if args:
      p = np.array(args)
    else:
      p = np.array([kwargs.pop(name) for name in self.param_names])

    if hasattr(self, 'normalize_params'):
      p = self.normalize_params(p)

    if self.delaunay is not None:
      if method=='Delaunay':
        weights, spectra = self.get_weights(p)
      elif method=='NaturalNeighbor':
        weights, spectra = compute_nn_weights(self.points, p)

      spectra = [s.crop(**kwargs) for s in spectra]
      yy = np.array([spec.ydata for spec in spectra])
      yy = np.sum(yy * weights[:, None], axis=0)
    else:
      i = np.searchsorted(self.points, p)[0]
      slc = slice(max(i-2,0), min(i+2,self.npoints))
      spectra = [s.crop(**kwargs) for s in self.spectra[slc]]
      yy = np.array([spec.ydata for spec in spectra])
      func = pchip_mono(self.points[slc], yy, axis=0)
      yy = func(p).squeeze()

    spec = SpectrumClass()
    spec.wave = spectra[0].wave
    data = {spectra[0].data.columns[0]:yy}
    spec.data = pd.DataFrame(data)
    return spec
    
  #=================================================================
  def fit(self, oxx, oyy, bounds=None, p0=None, **kwargs):
    #-----------------------------------------------------------------
    bounds = bounds or {}
    p0 = p0 or {}

    p = {}
    b = {}
    for i, name in enumerate(self.param_names):
      vals = self.params[:,i] if self.ndim > 1 else self.params
      p[name] = np.mean(vals)
      b[name] = (np.min(vals), np.max(vals))

    p.update(p0)
    b.update(bounds)
    return SpectrumClass.fit(self, oxx, oyy, bounds=b, p0=p, **kwargs)
  
  #=====================================================================
  def plot_nodes(self, fig=None, size=5, color='black', cmap='Viridis', opacity=0.8, **kwargs):
    import plotly.graph_objs as go
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
    import plotly.graph_objs as go
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