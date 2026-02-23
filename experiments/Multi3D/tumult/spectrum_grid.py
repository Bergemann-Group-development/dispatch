from .auxfuncs import pchip_mono
from .spectrum import Spectrum
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, ConvexHull, Voronoi
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

ensure_list = lambda x: x if isinstance(x, (list, np.ndarray, pd.Series)) else [x]

class SpectrumGrid(object):
    """
    A class to hold a grid of spectra.
    """
    #=====================================================================
    def __init__(self, spectra, normalize=True, **kwargs):
        self.spectra = ensure_list(spectra)
        self.normalize = normalize
        processed_kwargs = {k: ensure_list(v) for k, v in kwargs.items()}
        self.params = pd.DataFrame(processed_kwargs)
        self._build_grid()

    def _build_grid(self):
        """Builds the interpolation structure based on current dimensionality."""
        self.npoints = len(self.spectra)
        
        # If no params or only one point, we can't build a grid
        if self.params.empty or self.npoints <= 1:
            self.ndim = 0
            self.points = None
            self.delaunay = None
            if self.npoints == 1:
                self.reference = self.spectra[0]
            return

        self.ndim = self.params.shape[1]
        self.param_names = self.params.columns.tolist()

        if self.normalize:
            # Pandas peak-to-peak equivalent
            p_min = self.params.min()
            p_ptp = self.params.max() - p_min
            p_ptp = p_ptp.replace(0, 1.0) # Avoid division by zero
            
            self._p_min, self._p_ptp = p_min, p_ptp
            # Normalize for triangulation
            self.points = ((self.params - self._p_min) / self._p_ptp).values
            
            # Helper for crop method
            self.normalize_params = lambda p_arr: (p_arr - self._p_min.values) / self._p_ptp.values
        else:
            self.points = self.params.values

        if self.ndim > 1:
            try:
                self.delaunay = Delaunay(self.points)
            except Exception as e:
                print(f"Delaunay build failed: {e}")
                self.delaunay = None
        else:
            self.delaunay = None
            # 1D logic: Ensure sorted order for interpolation
            self.points = self.points.ravel()
            order = np.argsort(self.points)
            self.points = self.points[order]
            self.params = self.params.iloc[order].reset_index(drop=True)
            self.spectra = [self.spectra[i] for i in order]

    def add_data(self, new_spectra, tol=1e-8, **kwargs):
        """
        Adds spectra with dynamic parameter expansion.
        Example: grid.add_data(spec, x=10, y=20)
        """
        new_spectra = ensure_list(new_spectra)
        
        # Convert scalars to lists so DataFrame constructor doesn't complain
        processed_kwargs = {k: ensure_list(v) for k, v in kwargs.items()}
        new_params = pd.DataFrame(processed_kwargs)

        if len(new_params) != len(new_spectra):
            raise ValueError(f"Got {len(new_spectra)} spectra but {len(new_params)} param rows.")

        # Logic for first-time parameter addition (Original point = Origin)
        if self.params.empty and self.npoints == 1:
            # Create a zero-filled row for the existing spectrum
            self.params = pd.DataFrame(np.zeros((1, new_params.shape[1])), 
                                      columns=new_params.columns)

        # Merge and fill missing dimensions with 0
        combined = pd.concat([self.params, new_params], axis=0, ignore_index=True, sort=False).fillna(0)

        # Duplicate check on physical units
        if not self.params.empty:
            existing_count = len(self.params)
            unique_indices = []
            
            for i in range(len(new_params)):
                new_row_val = combined.iloc[existing_count + i].values
                existing_vals = combined.iloc[:existing_count].values
                dist = np.linalg.norm(existing_vals - new_row_val, axis=1)
                
                if np.any(dist < tol):
                    print(f"Point {new_row_val} is a duplicate. Skipping.")
                else:
                    unique_indices.append(i)
            
            if not unique_indices:
                return

            # Keep only the unique new spectra and the old data
            new_spectra = [new_spectra[i] for i in unique_indices]
            final_rows = list(range(existing_count)) + [existing_count + i for i in unique_indices]
            self.params = combined.iloc[final_rows].reset_index(drop=True)
        else:
            self.params = combined

        self.spectra.extend(new_spectra)
        self._build_grid()
  
    #=====================================================================
    def crop(self, *args, method='Delaunay', **kwargs):
        if self.params.empty:
            return self.spectra[0].crop(*args, **kwargs)

        p = args or [kwargs.pop(name) for name in self.param_names]
        p = np.array(p)

        if hasattr(self, 'normalize_params'):
            p = self.normalize_params(p)

        if self.delaunay is not None:
            if method=='Delaunay':
                weights, spectra = self.get_weights(p)
            elif method=='NaturalNeighbor':
                weights, spectra = compute_nn_weights(self.points, p)
            else:
                raise ValueError(f"Unknown interpolation method: {method}")

            spec = Spectrum.combine([s.crop(**kwargs) for s in spectra], weights=weights)
        else:
            i = np.searchsorted(self.points, p)[0]
            slc = slice(max(i-2,0), min(i+2,self.npoints))
            spectra = [s.crop(**kwargs) for s in self.spectra[slc]]
            yy = np.array([spec.ydata for spec in spectra])
            yy = pchip_mono(self.points[slc], yy, axis=0)(p)[0]
            spec = Spectrum(spectra[0].xdata, yy, xunit=spectra[0].xunit, yunit=spectra[0].yunit)

        return spec

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
    
    #=================================================================
    def fit(self, spec, bounds=None, p0=None, **kwargs):
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
        return Spectrum.fit(self, spec, bounds=b, p0=p, **kwargs)
  
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