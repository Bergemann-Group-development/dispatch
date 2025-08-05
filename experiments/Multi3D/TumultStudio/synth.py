from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .llist_annotations import get_llist_table
from tumult.auxfuncs import update_dict
import tumult
import numpy as np
import tempfile
import logging


def run_synth(folder, name, settings, line_mask=None, default=None, llist=None, progress=None):
    logging.info(f'Executing run {name}')
    folder.exists() or folder.mkdir()
    input_name = folder / name
    log_file = 'test.log'

    update_dict(settings, {'composition_params':{'abund_file':saved_settings['abundance file']},})

    if line_mask:
        lam = line_mask.get_wavelength_grid(daa=0.01)

    if not default:
        default = module_dir / 'default.nml'

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as llist_file:
        if llist is not None: llist.write_to_ts_format(llist_file.name)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as lam_file:
            if line_mask is not None:
                np.savetxt(lam_file.name, lam)
                settings['spectrum_params'] = {'lam_file':lam_file.name}

            if llist is not None:
                settings['linelist_params'] = {'line_lists':llist_file.name}

            @run_in_thread
            def execute_m3dis():
                folder, ok = tumult.launch(input_name, nml_settings=settings, default=default,
                                        log_file=log_file, threads=2, overwrite=True, verbose=0)
                return folder, ok
            
            thread = execute_m3dis()

            if progress:
                while thread.is_alive():
                    time.sleep(0.25)
                    monitor_synth(input_name, log_file, progress)
                    
            try:
                folder, ok = thread.get_result()
            except ValueError as e:
                print(f"Caught error: {e}")
                thread.print_output()  # Shows error details
        
    logging.info(f'Reading run from {folder}')
    run = tumult.read(folder, read_llist=True)
    run = StudioSynth.from_run(run)
    
    if line_mask is not None:
        run.line_mask = line_mask
    if hasattr(run, 'llist'):
        run.llist_html = get_llist_table(run)

    return run


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
    def from_run(cls, run, **kwargs):
        # Copy all attributes from run
        b = cls.__new__(cls)
        b.__dict__.update(run.__dict__)
        # Initialize any new attributes
        for key, value in kwargs.items():
            setattr(b, key, value)
        return b
    
    @property
    def label(self):
        return self.folder.name

class StudioCoG(tumult.CoG, StudioSynth):
    pass