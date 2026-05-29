from .ServerManager import update_source_data, clear_data
from .auxfuncs import *
from .widgets import HTMLDataFrame
from .StudioAtmos import StudioAtmos
from .synth import *
import ipyformkit as ifk
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, ConvexHull, Voronoi

class StudioGrid(object):
    icon = 'border-all'

    def __init__(self, label, data, files, format):
        super().__init__()
        self.df = pd.DataFrame(data)
        self.grid = self.normalise(self.df)
        self.files = np.array(files)
        self.format = format
        self.label = label

    def init_plot_form(self, studio):
        columns = tuple(self.df.columns)
        fields = {('x-axis', 'y-axis'): (columns, columns[1:])}
        plot_form = ifk.Form(fields)
        self.xaxis_wid = plot_form.inputs_dict['x-axis']
        self.yaxis_wid = plot_form.inputs_dict['y-axis']
        if len(columns) > 0:
            self.yaxis_wid.value = columns[1]

        @log_callback_errors
        def on_axis_change(change):
            self.yaxis_wid.options = [o for o in columns if o != self.xaxis_wid.value]
            self.draw_grid_data(studio)

        self.xaxis_wid.observe(on_axis_change, names='value')
        self.yaxis_wid.observe(on_axis_change, names='value')
        return plot_form

    def init_interp_form(self, studio):
        keys = tuple(self.df.iloc[0].index)
        vals = tuple(self.df.iloc[0].values)
        self.interp_keys = keys

        if 'vmic' not in [k.lower() for k in keys]:
            keys = ('vmic',) + keys
            vals = (0.0,) + vals
        if 'feh' not in [k.lower() for k in keys]:
            keys = ('FeH',) + keys
            vals = (0.0,) + vals

        keys = chunk_tuple(keys, chunk_size=3)
        vals = chunk_tuple(vals, chunk_size=3)
        interp_fields = {k:v for k,v in zip(keys, vals)}
        interp_fields['interpolate'] = None
        atmos_interp_form = ifk.Form(interp_fields)
        interp_btn = atmos_interp_form.inputs_dict['interpolate']

        @log_callback_errors
        def on_val_change(change):
            xval = self.xaxis_wid.value
            yval = self.yaxis_wid.value
            xwid = atmos_interp_form.inputs_dict[xval]
            ywid = atmos_interp_form.inputs_dict[yval]
            self.indices, self.wts = self.get_nbors()
            selection = self.df.loc[self.indices][[xval, yval]].drop_duplicates()
            indices = self.plot_df.reset_index().merge(selection)['index'].tolist()
            @update_source_data(studio.atmos_interp_source)
            def update_atmos_interp(data):
                data['x'] = [xwid.value]
                data['y'] = [ywid.value]
                studio.atmos_scatter_source.selected.indices = indices
            update_atmos_interp()

        for wid in [v for v in atmos_interp_form.inputs_dict.values() if v is not interp_btn]:
            wid.observe(on_val_change, names='value')

        self.xaxis_wid.observe(on_val_change, names='value')
        self.yaxis_wid.observe(on_val_change, names='value')

        @log_callback_errors
        def on_interp_click(btn):
            studio.add_object('Atmos', self.interpolate_atmosphere)

        interp_btn.on_click(on_interp_click)
        return atmos_interp_form

    def normalise(self, val):
        normed = (val - self.df.min()) / (self.df.max() - self.df.min())
        return normed.values

    def chosen(self, studio):
        if not hasattr(self, 'plot_form'):
            self.plot_form = self.init_plot_form(studio)

        if not hasattr(self, 'atmos_interp_form'):
            self.atmos_interp_form = self.init_interp_form(studio)

        self.draw_grid_data(studio)
        instructions = {'tabs': (self.atmos_interp_form, self.plot_form, self.info),
                        'titles': ('Interpolate atmosphere', 'Plot Settings', 'Info')}
        return instructions
    
    @property
    def info(self):
        df = self.df.copy()
        df['file'] = [f.name for f in self.files]
        return HTMLDataFrame(df)
    
    def draw_grid_data(self, studio):
        clear_data(studio.atmos_line_source)
        xval = self.xaxis_wid.value
        yval = self.yaxis_wid.value
        self.plot_df = self.df[[xval, yval]].drop_duplicates().reset_index(drop=True)

        @update_source_data(studio.atmos_scatter_source)
        def update_atmos_scatter(data):
            data['x'] = self.plot_df[xval].values
            data['y'] = self.plot_df[yval].values
            studio.atmos_fig.xaxis.axis_label = xval
            studio.atmos_fig.yaxis.axis_label = yval
        update_atmos_scatter()

    def get_nbors(self):
        values = self.atmos_interp_form.check_and_return_values()
        query = [values[k] for k in self.interp_keys]
        query = self.normalise(query)

        subgrid, subidx = keep_k_closest(self.grid, query)
        flat = np.all(subgrid == subgrid[0], axis=0)
        subgrid = subgrid[:, ~flat]
        query = query[~flat]
        wts, indices = compute_nn_weights(subgrid, query)

        if any(indices < 0):
            logging.info(f'Query outside of Convex Hull!')
            return np.array([]), np.array([])
        
        mask = wts > 1e-2*np.max(wts)
        wts = wts[mask]
        wts = wts / np.sum(wts)
        indices = indices[mask]
        indices = subidx[indices]
        return indices, wts

    def interpolate_atmosphere(self, studio, name, runner=None):
        values = self.atmos_interp_form.check_and_return_values()
        query = [values[k] for k in self.interp_keys]
        query = self.normalise(query)

        runs = []
        for i, file in enumerate(self.files[self.indices]):
            settings = {'atmos_params':
                            {'atmos_file':str(file),
                            'atmos_format':self.format,},
                        'linelist_params':
                            {'line_lists':''},}
            if 'vmic' in self.interp_keys:
                settings['atmos_params']['vmic'] = values['vmic']
            if 'FeH' in self.interp_keys:
                settings['atmos_params']['feh'] = values['FeH']

            run_name = f'{name}_{i}' 
            run = run_synth(studio, run_name, settings, runner=runner)
            if run is None: return
            runs.append(run)
        
        comb = tumult.combine(runs, weights=self.wts)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as atmos_file:
            comb.to_dispatch1d_atmos(atmos_file.name, vmic=values['vmic'], FeH=values['FeH'])
            settings = {'atmos_params':
                            {'atmos_file':atmos_file.name,
                            'atmos_format':'Text',},
                        'linelist_params':
                            {'line_lists':''},}
            run = run_synth(studio, name, settings, runner=runner)
            if run:
                run.nbors = runs
                run.nbors_wts = self.wts
                new_atmos = StudioAtmos.from_run(run)
                return new_atmos
            

def keep_k_closest(points, query, k=100):
    """Keep only the k closest points to query_point"""
    # Calculate distances
    distances = np.linalg.norm(points - query, axis=1)
    
    # Get indices of k smallest distances (O(n) operation)
    k_closest_indices = np.argpartition(distances, k)[:k]
    
    return points[k_closest_indices], k_closest_indices


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

