import threading
from .utils import pchip_mono
from .spectrum import Spectrum
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, ConvexHull, Voronoi, HalfspaceIntersection, cKDTree
from scipy.optimize import linprog

ensure_list = lambda x: x if isinstance(x, (list, np.ndarray, pd.Series)) else [x]

class SpectrumGrid(object):
    """
    A class to hold a grid of spectra.

    ``separable`` names parameter axes that form a full tensor-product lattice with the
    remaining ("scatter") axes — every scatter node must carry a spectrum at every
    combination of the separable axes' values (validated at build time, loud ValueError
    otherwise). Separable axes are interpolated per-pixel with monotone-cubic
    ``pchip_mono`` (shape-preserving, no overshoot — the right tool for smooth monotone
    responses like flux vs. abundance or microturbulence), while the scatter axes keep
    the linear Delaunay / natural-neighbour weights. The names are whatever the caller
    used as grid kwargs; nothing here attaches meaning to specific axis names.

    ``separable=None`` (default) auto-flags the axis of a single-axis grid, preserving
    the historical behaviour where 1-D grids interpolate with pchip; pass
    ``separable=[]`` to force linear interpolation there too.
    """
    #=====================================================================
    def __init__(self, spectra, normalize=True, separable=None, **kwargs):
        self.spectra = ensure_list(spectra)
        self.normalize = normalize
        self.separable = list(separable) if separable is not None else None
        processed_kwargs = {k: ensure_list(v) for k, v in kwargs.items()}
        self.params = pd.DataFrame(processed_kwargs)
        self.p0 = {}
        # node weights at p0 fixed by the caller (the synthesis node selection); a crop
        # without a query point uses them instead of interpolating afresh
        self.interp_weights = None
        self._build_grid()

    def _build_grid(self):
        """Builds the interpolation structure based on current dimensionality."""
        self.npoints = len(self.spectra)
        self._sep_names, self._sep_axes = [], {}
        self._sep_lattice, self._scatter_names, self._scatter_points = None, [], None

        # If no params or only one point, we can't build a grid
        if self.params.empty or self.npoints <= 1:
            self.ndim = 0
            self.points = None
            if self.npoints == 1:
                self.reference = self.spectra[0]
            return

        self.param_names = self.params.columns.tolist()

        if self.normalize:
            points_full, self._p_min, self._p_ptp = normalize_points(self.params.values)
        else:
            points_full = self.params.values
            self._p_min = np.zeros(self.params.shape[1])
            self._p_ptp = np.ones(self.params.shape[1])

        # Drop uniform columns — they collapse the point cloud and break Delaunay
        varying = np.ptp(points_full, axis=0) > 0
        self.points = points_full[:, varying]
        self.normalize_params = lambda p_arr: ((p_arr - self._p_min) / self._p_ptp)[varying]
        self.ndim = self.points.shape[1]
        self.p0 = {name: self.params[name].median() for name in self.param_names}

        varying_names = [n for n, v in zip(self.param_names, varying) if v]
        sep_names = self._resolve_separable(varying_names)

        # No triangulation here: crop triangulates the node subset each query needs
        # (get_interp_weights), so the grid only keeps the points.
        if self.ndim == 1:
            # 1D logic: Ensure sorted order for interpolation
            self.points = self.points.ravel()
            order = np.argsort(self.points)
            self.points = self.points[order]
            self.params = self.params.iloc[order].reset_index(drop=True)
            self.spectra = [self.spectra[i] for i in order]

        if sep_names:
            self._build_separable(sep_names, varying_names)

    #=====================================================================
    def _resolve_separable(self, varying_names):
        """The separable axes actually in play: the declared names restricted to the
        varying columns (uniform axes are pinned constants anyway), or — with no
        declaration — the single axis of a 1-D grid (the historical pchip behaviour)."""
        if self.separable is None:
            return list(varying_names) if len(varying_names) == 1 else []
        unknown = [n for n in self.separable if n not in self.param_names]
        if unknown:
            raise ValueError(
                f"separable axes {unknown} are not grid parameters {self.param_names}")
        return [n for n in self.separable if n in varying_names]

    #=====================================================================
    def _build_separable(self, sep_names, varying_names):
        """Index the spectra as a tensor product: scatter node x separable lattice.

        Fails loudly when the points do not form a complete lattice (every scatter node
        must carry every combination of the separable axes' values, exactly once)."""
        scatter_names = [n for n in varying_names if n not in sep_names]
        axes = {n: np.unique(self.params[n].to_numpy(dtype=float)) for n in sep_names}
        shape = tuple(axes[n].size for n in sep_names)

        if scatter_names:
            scat = self.params[scatter_names].to_numpy(dtype=float)
            scatter_points, inverse = np.unique(scat, axis=0, return_inverse=True)
            inverse = np.asarray(inverse).ravel()
        else:
            scatter_points = np.zeros((1, 0))
            inverse = np.zeros(self.npoints, dtype=int)

        n_lattice = int(np.prod(shape))
        if len(scatter_points) * n_lattice != self.npoints:
            raise ValueError(
                f"separable axes {sep_names} do not form a tensor product: "
                f"{len(scatter_points)} scatter nodes x {n_lattice} lattice points "
                f"!= {self.npoints} spectra")

        lattice = np.full((len(scatter_points),) + shape, -1, dtype=int)
        for i in range(self.npoints):
            idx = tuple(int(np.searchsorted(axes[n], float(self.params[n].iloc[i])))
                        for n in sep_names)
            if lattice[(inverse[i],) + idx] >= 0:
                raise ValueError(
                    f"duplicate grid point: scatter node {scatter_points[inverse[i]]} "
                    f"already has a spectrum at {dict(zip(sep_names, idx))}")
            lattice[(inverse[i],) + idx] = i

        self._sep_names = sep_names
        self._sep_axes = axes
        self._sep_lattice = lattice
        self._scatter_names = scatter_names
        self._scatter_points = scatter_points

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

        # Duplicate check in parameter units, against the stored points and within the
        # batch (the first of several coincident new points is kept)
        existing_count = len(self.params)
        values = combined.to_numpy(dtype=float)
        new_vals = values[existing_count:]
        duplicate = np.zeros(len(new_vals), dtype=bool)
        if existing_count:
            dist, _ = cKDTree(values[:existing_count]).query(new_vals, k=1)
            duplicate |= dist < tol
        for i, j in cKDTree(new_vals).query_pairs(tol):
            duplicate[max(i, j)] = True
        for row in new_vals[duplicate]:
            print(f"Point {row} is a duplicate. Skipping.")
        unique_indices = np.flatnonzero(~duplicate).tolist()
        if not unique_indices:
            return

        # Keep only the unique new spectra and the old data
        new_spectra = [new_spectra[i] for i in unique_indices]
        final_rows = list(range(existing_count)) + [existing_count + i for i in unique_indices]
        self.params = combined.iloc[final_rows].reset_index(drop=True)

        self.spectra.extend(new_spectra)
        self._build_grid()

    def find(self, tol=1e-8, **params):
        """The spectrum stored at ``params`` (within ``tol``, in parameter units), or None.

        Every grid parameter must be given; unknown names match nothing.
        """
        if self.params.empty or set(params) != set(self.params.columns):
            return None
        names = list(self.params.columns)
        values = self.params[names].to_numpy(dtype=float)
        query = np.array([float(params[n]) for n in names])
        dist = np.linalg.norm(values - query, axis=1)
        i = int(np.argmin(dist))
        return self.spectra[i] if dist[i] < tol else None
  
    #=====================================================================
    def crop(self, *args, method='Delaunay', **kwargs):
        if self.interp_weights is not None and not args \
                and not set(self.params.columns) & set(kwargs):
            return Spectrum.combine([s.crop(**kwargs) for s in self.spectra],
                                    weights=self.interp_weights)
        if self.params.empty:
            return self.spectra[0].crop(*args, **kwargs)

        p = args or [kwargs.pop(name, self.p0[name]) for name in self.param_names]
        query = {name: float(v) for name, v in zip(self.param_names, p)}

        if self._sep_names:
            return self._crop_separable(query, method, kwargs)

        p = np.array([query[name] for name in self.param_names])
        if hasattr(self, 'normalize_params'):
            p = self.normalize_params(p)

        # shared node selection: drops on-grid axes, exact NN / barycentric weights,
        # fast local stencil + tail truncation for NaturalNeighbor; handles 1-D linearly
        points = self.points if self.points.ndim == 2 else self.points[:, None]
        weights, indices = get_interp_weights(points, p, method)

        spectra = [self.spectra[i] for i in indices]
        return Spectrum.combine([s.crop(**kwargs) for s in spectra], weights=weights)

    #=====================================================================
    def _crop_separable(self, query, method, kwargs):
        """Tensor-product interpolation: linear weights across the scatter nodes,
        per-pixel monotone-cubic (pchip) collapse along the separable axes."""
        if self._scatter_names:
            q = np.array([query[n] for n in self._scatter_names], dtype=float)
            weights, sel = get_interp_weights(self._scatter_points, q, method)
        else:
            weights, sel = np.array([1.0]), np.array([0])

        collapsed = [self._collapse_node(int(si), query, kwargs) for si in sel]
        if not self._scatter_names:
            return collapsed[0]
        return Spectrum.combine(collapsed, weights=weights)

    #=====================================================================
    def _collapse_node(self, node, query, kwargs):
        """Collapse one scatter node's separable lattice to a single spectrum, one axis
        at a time, each through the same local 4-point pchip stencil the 1-D grid always
        used (pchip is node-local, so the window equals the full-axis interpolant within
        the bracketing interval; beyond the axis ends it degrades to the 2-point line)."""
        sub = self._sep_lattice[node]
        windows = []
        for n in self._sep_names:
            vals = self._sep_axes[n]
            i = int(np.searchsorted(vals, query[n]))
            windows.append(slice(max(i - 2, 0), min(i + 2, vals.size)))
        sub = sub[tuple(windows)]

        flat = sub.ravel()
        specs = [self.spectra[j].crop(**kwargs) for j in flat]
        yy = np.array([s.ydata for s in specs])
        yy = yy.reshape(sub.shape + yy.shape[1:])
        for w, n in zip(windows, self._sep_names):
            yy = pchip_mono(self._sep_axes[n][w], yy, axis=0)([query[n]])[0]

        ref = specs[0]
        return Spectrum(ref.xdata, yy, xunit=ref.xunit, yunit=ref.yunit)
    
    #=================================================================
    def fit(self, spec, bounds=None, p0=None, **kwargs):
        #-----------------------------------------------------------------
        bounds = bounds or {}
        p0 = p0 or {}

        p = {}
        b = {}
        for name in getattr(self, "param_names", []):   # empty for a single-node grid
            vals = self.params[name]
            p[name] = np.median(vals)
            b[name] = (np.min(vals), np.max(vals))

        p.update(p0)
        b.update(bounds)
        return Spectrum.fit(self, spec, bounds=b, p0=p, **kwargs)


#=====================================================================
def normalize_points(points):
    """Normalize rows to [0,1] per column. Returns (points_norm, p_min, p_ptp)."""
    points = np.asarray(points)
    p_min = points.min(axis=0)
    p_ptp = points.max(axis=0) - p_min
    p_ptp[p_ptp == 0] = 1.0
    return (points - p_min) / p_ptp, p_min, p_ptp

#=====================================================================
def get_delaunay_weights(points, query_point, delaunay=None, verbose=True):
    """Returns (weights, indices) for the enclosing simplex, or (None, None)."""
    if delaunay is None:
        delaunay = Delaunay(points)
    simplex = delaunay.find_simplex(query_point)
    if simplex == -1:
        return None, None
    vertices = delaunay.simplices[simplex]
    A = np.hstack((points[vertices], np.ones((len(vertices), 1))))
    b = np.append(query_point, 1)
    weights = np.linalg.solve(A.T, b)
    mask = weights > 1e-10
    weights = weights[mask]
    vertices = vertices[mask]
    if vertices is None:
        if verbose: print("Query point is outside the convex hull.")
        return None, None
    else:
        return weights, vertices

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
def _nn_halfspaces(site, others):
    """Bisector half-spaces {x: ||x-site|| <= ||x-o||} for each ``o`` in ``others``.

    Returned as rows ``[a_0..a_{d-1}, b]`` with the convention ``a·x + b <= 0`` (the
    half-space of points at least as close to ``site`` as to ``o``).
    """
    others = np.atleast_2d(np.asarray(others, dtype=float))
    site = np.asarray(site, dtype=float)
    a = others - site                                       # normal points toward ``o``
    b = -0.5 * (np.einsum("ij,ij->i", others, others) - site @ site)
    return np.column_stack([a, b])


def _chebyshev_center(halfspaces):
    """Interior point and inscribed-ball radius of ``{x: A x + b <= 0}``.

    Radius ``<= 0`` (or a ``None`` point) means the region is empty or lower-dimensional,
    i.e. the candidate neighbour steals no volume.
    """
    A, b = halfspaces[:, :-1], halfspaces[:, -1]
    d = A.shape[1]
    norms = np.linalg.norm(A, axis=1)
    c = np.zeros(d + 1); c[-1] = -1.0                       # maximise the radius r
    res = linprog(c, A_ub=np.column_stack([A, norms]), b_ub=-b,
                  bounds=[(None, None)] * d + [(0, None)], method="highs")
    if not res.success:
        return None, 0.0
    return res.x[:d], float(res.x[-1])


def _polytope_volume(halfspaces, interior):
    """``d``-volume and vertices of the convex polytope ``{A x + b <= 0}``."""
    try:
        verts = HalfspaceIntersection(halfspaces, np.asarray(interior, float)).intersections
        return float(ConvexHull(verts).volume), verts
    except Exception:
        return 0.0, None


def get_nn_weights(points, query, draw=False, hull=None, vor1=None):
    """Sibson natural-neighbour weights/indices for ``query`` within ``points``.

    Computed exactly from the Sibson definition: the weight of neighbour ``n`` is the
    fraction of the query's *inserted* Voronoi cell ``V_new(q)`` stolen from ``n``'s old
    cell, ``vol(V_new(q) ∩ V_old(n)) / vol(V_new(q))``. Within ``V_new(q)`` the only
    competitors are ``q``'s own Voronoi neighbours, so both regions are convex polytopes
    of bisector half-spaces among ``q`` and those candidates, and each stolen volume is a
    half-space intersection. This is robust on degenerate (regular/cocircular) grids,
    where the previous area-stealing heuristic over-counted next-ring nodes and lost
    linear precision; here ``sum(w_i * p_i) == q`` to machine precision.

    ``hull`` (ConvexHull of ``points``) is reused across queries when supplied. ``vor1``
    is accepted for backward compatibility and no longer used.
    """
    if draw:
        # only the debug scatter path needs matplotlib; deferred so importing this
        # module doesn't pay its ~250ms startup cost
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon

    points = np.asarray(points, dtype=float)
    query = np.asarray(query, dtype=float)

    idx = is_point_in_array(query, points)
    if idx >= 0:
        return np.array([1.0]), np.array([idx])

    if hull is None:
        hull = ConvexHull(points)
    if not point_in_convex_hull(query, hull):
        return np.array([0.0]), np.array([-1])

    # candidate natural neighbours: the Voronoi neighbours of the inserted query point
    vor2 = Voronoi(np.vstack([points, query]))
    nbors = np.unique(get_nbors(vor2, -1)[0])
    cand = points[nbors]
    scale = float(np.linalg.norm(np.ptp(points, axis=0))) or 1.0
    tol = 1e-9 * scale

    hq = _nn_halfspaces(query, cand)                        # V_new(q): q beats every candidate
    weights, keep, polys = [], [], []
    for k, n in enumerate(nbors):
        others = np.delete(cand, k, axis=0)                 # n is the nearest *old* site
        hs = np.vstack([hq, _nn_halfspaces(points[n], others)]) if len(others) else hq
        interior, radius = _chebyshev_center(hs)
        if interior is None or radius < tol:
            continue
        vol, verts = _polytope_volume(hs, interior)
        if vol > 0:
            weights.append(vol); keep.append(int(n)); polys.append(verts)

    weights = np.asarray(weights, dtype=float)
    if weights.sum() <= 0:
        return np.array([0.0]), np.array([-1])
    weights /= weights.sum()

    if draw and points.shape[1] == 2:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, verts in enumerate(polys):
            v = verts[ConvexHull(verts).vertices]
            plt.gca().add_patch(Polygon(v, color=colors[i % len(colors)], alpha=0.5))

    return weights, np.asarray(keep)


#=====================================================================
# Node selection: pick the interpolation stencil (weights + indices) for a query.
#=====================================================================

# natural-neighbour weights have a long tail of negligible nodes; keep only the highest-
# weight ones up to this cumulative fraction (the rest barely affect the result but each
# costs a full synthesis). Delaunay weights (a single simplex) are never truncated.
_NN_WEIGHT_KEEP = 0.995

# a query coordinate within this fraction of an axis' span of a grid value is treated as
# lying exactly on the grid (so that dimension is dropped from the interpolation).
_ON_GRID_TOL = 1e-4

# cache of per-subgrid structures, keyed by the subgrid's normalised-point bytes. The
# subgrid (and hence its triangulation) repeats across queries that pin the same axes to
# the same grid values, e.g. scanning Teff/logg at a fixed [Fe/H] and vmic.
_SUBGRID_CACHE = {}


def _subgrid(points):
    """Cached normalisation + lazily-built Delaunay / NN structures for a node set.

    ``points`` is a raw (n, d) array; returns a dict with ``pn`` (normalised), ``pmin``,
    ``pptp`` and slots ``delaunay`` / ``nn`` filled on first use by the method that needs them.
    """
    key = (points.shape, points.tobytes())
    cached = _SUBGRID_CACHE.get(key)
    if cached is None:
        pn, pmin, pptp = normalize_points(points)
        if len(_SUBGRID_CACHE) >= 8:
            _SUBGRID_CACHE.clear()
        cached = _SUBGRID_CACHE[key] = {'pn': pn, 'pmin': pmin, 'pptp': pptp,
                                        'delaunay': None, 'nn': None}
    return cached


def _nn_select(points_norm, query_norm, grid_nn, max_nb=100, max_sub=2000, max_ring=4):
    """Natural-neighbour node indices for the query — fast local path, or None to defer.

    The query's natural neighbours are determined entirely by its local Voronoi cell, so:
    take the nearest grid node, grow its Voronoi-adjacency neighbourhood ring by ring, and
    build a small local Voronoi from it. Once the query's natural neighbours are all
    *interior* to that neighbourhood (none on its convex hull), enlarging it can't change
    them — the result is exact. Returns ``None`` (caller does the full grid computation)
    for a degenerate, out-of-local-hull, or huge-stencil (sparse-region) query.
    """
    vor1, kdtree, _ = grid_nn
    rp = vor1.ridge_points
    host = int(kdtree.query(query_norm)[1])
    nodes, frontier = {host}, {host}
    for r in range(1, max_ring + 1):
        nxt = set()
        for n in frontier:
            m = np.any(rp == n, axis=1)
            nxt.update(int(x) for x in (rp[m].sum(1) - n))
        nxt -= nodes
        nodes |= nxt
        frontier = nxt
        if r < 2:
            continue
        if len(nodes) > max_sub:
            return None
        loc = np.fromiter(sorted(nodes), dtype=int)
        try:
            w, li = get_nn_weights(points_norm[loc], query_norm)
        except Exception:
            continue   # degenerate local subset — grow another ring
        if li[0] < 0 or len(li) > max_nb:
            return None
        hull_v = {int(v) for v in ConvexHull(points_norm[loc]).vertices}
        if not any(int(i) in hull_v for i in li):
            return w, loc[li]   # (weights, indices); all natural neighbours interior -> exact
    return None


def _truncate_weights(weights, indices, keep):
    """Keep the highest-weight nodes whose cumulative weight reaches ``keep``, renormalised.

    Returns ``(weights, indices)`` for the kept nodes (at least one). Dropping the
    negligible tail removes syntheses that barely change the interpolated spectrum.
    """
    weights = np.asarray(weights, dtype=float)
    indices = np.asarray(indices)
    order = np.argsort(weights)[::-1]                      # high -> low
    csum = np.cumsum(weights[order])
    n = min(max(int(np.searchsorted(csum, keep * csum[-1])) + 1, 1), len(order))
    sel = order[:n]
    w = weights[sel]
    return w / w.sum(), indices[sel]


_INTERP_LOCK = threading.RLock()


def get_interp_weights(points, query, method='Delaunay', keep=_NN_WEIGHT_KEEP):
    """Thread-safe access to cached interpolation geometry."""
    with _INTERP_LOCK:
        return _get_interp_weights(points, query, method, keep)


def _get_interp_weights(points, query, method='Delaunay', keep=_NN_WEIGHT_KEEP):
    """(weights, indices) of the interpolation nodes for ``query`` in raw grid ``points``.

    The fundamental node-selection entry point: it pins on-grid axes (dimension reduction),
    dispatches to Delaunay or natural-neighbour weights, uses the fast local stencil for
    NN, and truncates the negligible NN tail.

    On a regular grid, a query coordinate that lands exactly on a grid value makes that a
    redundant dimension — interpolating in it wastes nodes (averaging values that cancel)
    and makes local triangulations degenerate. So pin such axes to their grid value,
    restrict to the matching nodes, and interpolate only in the remaining (active) dimensions.

    ``method`` is ``'Delaunay'`` or ``'NaturalNeighbor'``; ``keep`` is the cumulative-weight
    fraction retained for NN (ignored by Delaunay). Indices index into ``points``.
    """
    points = np.asarray(points, dtype=float)
    n, ndim = points.shape
    query = np.asarray(query, dtype=float)

    active, pinned = [], []
    for d in range(ndim):
        vals = np.unique(points[:, d])
        nearest = vals[int(np.argmin(np.abs(vals - query[d])))]
        span = float(np.ptp(vals)) or 1.0
        if abs(nearest - query[d]) <= _ON_GRID_TOL * span:
            pinned.append((d, nearest))
        else:
            active.append(d)

    mask = np.ones(n, dtype=bool)
    for d, v in pinned:
        mask &= np.isclose(points[:, d], v)
    sub = np.flatnonzero(mask)

    if not active:                       # query coincides with a grid node
        return np.array([1.0]), sub[:1]

    outside = "Target is outside the convex hull of available atmospheres."
    P = points[np.ix_(sub, active)]

    if len(active) == 1:                 # single off-grid axis -> 1-D linear interp
        xs = P[:, 0]                     # (Delaunay/Voronoi need >=2-D, so handle 1-D here)
        order = np.argsort(xs)
        xs_s = xs[order]
        qx = float(query[active][0])
        j = max(1, min(int(np.searchsorted(xs_s, qx)), len(xs_s) - 1))
        x0v, x1v = xs_s[j - 1], xs_s[j]
        if np.isclose(x1v, x0v):
            return np.array([1.0]), sub[order[[j - 1]]]
        t = float(np.clip((qx - x0v) / (x1v - x0v), 0.0, 1.0))
        return np.array([1.0 - t, t]), sub[order[[j - 1, j]]]

    c = _subgrid(P)
    pn, qn = c['pn'], (query[active] - c['pmin']) / c['pptp']

    if method == 'Delaunay':
        if c['delaunay'] is None:
            c['delaunay'] = Delaunay(pn)
        w, li = get_delaunay_weights(pn, qn, c['delaunay'])
        if li is None:
            raise ValueError(outside)
        return w, sub[li]

    if method == 'NaturalNeighbor':
        if c['nn'] is None:
            c['nn'] = (Voronoi(pn), cKDTree(pn), ConvexHull(pn))
        sel = _nn_select(pn, qn, c['nn'])
        if sel is None:                  # local shortcut declined -> full subgrid
            vor1, _kd, hull = c['nn']
            w, li = get_nn_weights(pn, qn, hull=hull, vor1=vor1)
        else:
            w, li = sel
        if len(li) == 0 or li[0] < 0:
            raise ValueError(outside)
        w, li = _truncate_weights(w, li, keep)
        return w, sub[li]

    raise ValueError(f"Unknown method: {method!r}. Use 'Delaunay' or 'NaturalNeighbor'.")