from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .llist_annotations import get_llist_table
from .fit_form import fit_form, set_fit_form
from .rerun_form import rerun_form, set_rerun_form
from .llist_annotations import plot_line_annotations
from . import ServerManager
import tumult
import numpy as np
import tempfile
import logging
import copy
from tumult.auxfuncs import update_dict


def run_synth(studio, name, settings, line_mask=None, default=None, llist=None, folder=None,
              runner=None):
    logging.info(f'Executing run {name}')

    if folder:
        folder.mkdir(parents=True, exist_ok=True)
    else:
        folder = studio.synth_folder

    input_name = folder / name
    log_file = 'test.log'

    update_dict(settings, {'composition_params':{'abund_file':saved_settings['abundance file']},})

    if line_mask:
        lam = line_mask.get_wavelength_grid(daa=0.01)

    if not default:
        default = module_dir / 'default.nml'

    process = None
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as llist_file:
        if llist is not None: llist.write_to_tumult_format(llist_file.name)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as lam_file:
            if line_mask is not None:
                np.savetxt(lam_file.name, lam)
                settings['spectrum_params'] = {'lam_file':lam_file.name}

            if llist is not None:
                settings['linelist_params'] = {'line_lists':llist_file.name}

            ofolder, process = tumult.launch(input_name, nml_settings=settings, default=default,
                                    log_file=log_file, threads=2, overwrite=True, verbose=0,
                                    runner=runner)

    if process and process.returncode == 0:
        logging.info(f'Reading run from {ofolder}')
        run = tumult.read(ofolder, read_llist=True)
        run.studio = studio
        if line_mask:
            run = StudioSynth.from_run(run, source=studio.synth_source, line_mask=line_mask)

        return run
    else:
        return None


def monitor_synth(name, log_file, progress):
    folder = pathlib.Path('data') / name
    n, nmax = tumult.monitor(folder, log_file=log_file)
    old_nmax = progress.max
    old_n = progress.value
    new_nmax = nmax+1

    # we always want some progress to be visible
    nmin = np.ceil(new_nmax*0.01)
    # do not drop below previous progress
    nprev = int(old_n*new_nmax/old_nmax)
    new_n = max(n, nprev, nmin)

    progress.value = new_n
    progress.max = new_nmax


class StudioSynth(tumult.tumult):
    @classmethod
    #=================================================================
    def from_run(cls, run, source=None, **kwargs):
        # Copy all attributes from run
        b = cls.__new__(cls)
        b.__dict__.update(run.__dict__)
        # Initialize any new attributes
        for key, value in kwargs.items():
            setattr(b, key, value)

        if source:
            b.source = source
            b._init_source()

        return b
    
    @property
    #=================================================================
    def label(self):
        return self.folder.name
    
    #=================================================================
    def delete(self):
        pass
    
    #=================================================================
    def _init_source(self):
        rid= id(self)

        @update_source_data(self.source)
        def delete_synth(data):
            i = data['run_ids'].index(rid)
            for key in data:
                data[key].pop(i)
        
        @update_source_data(self.source)
        def hide_synth(data):
            i = data['run_ids'].index(rid)
            data['alphas'][i] = 0.0

        @update_source_data(self.source)
        def show_synth(data):
            i = data['run_ids'].index(rid)
            data['alphas'][i] = 1.0

        self.delete = delete_synth
        self.hide = hide_synth
        self.show = show_synth
    
    #=================================================================
    def chosen(self, studio):
        lm = self.line_mask
        ids = studio.synth_source.data['run_ids']
        i = ids.index(id(self))
        reorder_data(studio.synth_source, i, len(ids))  # move to back to display on top
        
        set_fit_form(self)
        set_rerun_form(self)

        # it is important that this happens on the main thread
        # probably because of javascript issues on worker threads
        if hasattr(self, 'llist') and not hasattr(self, 'llist_html'):
            self.llist_html = get_llist_table(self.llist)

        llist_vbox = widgets.VBox([self.llist_html, rerun_form.vbox])
        llist_vbox.add_class('llist-vbox')
        instructions = {'xmin': lm.xmin,
                        'xmax': lm.xmax,
                        'tabs': (fit_form.vbox, llist_vbox),
                        'titles': ('Fit Settings', 'Line List')}
        return instructions
    
    #=================================================================
    def convert_to_CoG(self, eid):
        off = 0
        suffix = f'_{off*100:+04.0f}'
        parent, name = (self.folder.parent, self.folder.name)
        cog_folder = unique_folder_name(parent / f'{eid}_CoG_{name}')
        logging.info(f'Creating folder {cog_folder}.')
        cog_folder.mkdir()
        new_folder = unique_folder_name(cog_folder / f'{name}{suffix}')
        logging.info(f'Renaming {self.folder} to {new_folder}.')
        self.rename_folder(new_folder)
        cog_name = cog_folder.name.split(suffix)[0]
        a0 = self.abundances[eid]
        cog = StudioCoG([self], cog_name, eid, a0, abnds=[a0])
        cog.atmos = self.atmos
        return cog
    
    #=================================================================
    def rerun(self, studio, name, settings, **kwargs):
        default = self.folder
        line_mask = self.line_mask
        llist = tumult.LineList(self.llist_html.get_dataframe())
        folder = self.folder.parent.relative_to(self.datadir)
        run = run_synth(studio, name, settings, line_mask, default, llist, folder, **kwargs)
        if run:
            if hasattr(self, 'fit_params'):
                run.fit_params = copy.deepcopy(self.fit_params)
            self.llist_html.reset_html()
            run.atmos = self.atmos
        return run
    
    #=================================================================
    def plot_llist(self, fig, source, fig_width):
        return plot_line_annotations(self, fig, source, fig_width)


class StudioCoG(tumult.CoG, StudioSynth):
    #=================================================================
    def __init__(self, runs, cog_name, eid, a0, **kwargs):
        super().__init__(runs, **kwargs)
        self.cog_name = cog_name
        self.a0 = a0
        self.eid = eid

    #=================================================================
    def add_abnd(self, anew, tol=0.005, runner=None):
        if anew in self.abnds:
            logging.info(f'Abundance {anew} already exists in CoG {self.cog_name}.')
            return self
        elif np.min(abs(self.abnds - anew)) < tol:
            logging.info(f'Abundance {anew} is too close to existing abundances in CoG {self.cog_name}.')
            return self
        
        off = anew - self.a0
        suffix = f'_{off*100:+04.0f}'
        name = self.cog_name + suffix
        settings = {'composition_params': {'abundances': f'{self.eid}={anew:.2f}'}}
        run = self.rerun(self.studio, name, settings, runner=runner)
        if run:
            self.add_run(run, abnd=anew)
        
    
    #=================================================================
    def iterate_fit(self, oxx, oyy, aguess=None, prec=0.005, runner=None, **kwargs):
        p, _ = self.fit_observed(oxx, oyy, **kwargs)
        afit = p['abnd']

        if aguess is None:
            aguess = find_nearest(self.abnds, afit)

        while (abs(afit - aguess) > prec):
            aguess = afit
            self.add_abnd(aguess, runner=runner)
            p, _ = self.fit_observed(oxx, oyy, **kwargs)
            afit = p['abnd']

        return p
    
    #=================================================================
    def fit_observed(self, oxx, oyy, bounds=None, **kwargs):
        #-----------------------------------------------------------------
        bounds = bounds or {}
        b = {'abnd':(-10.0, 12.0)}
        b.update(bounds)
        return super().fit_observed(oxx, oyy, bounds=b, **kwargs)
    
    @property
    def label(self):
        return f'{self.cog_name}'
        
        