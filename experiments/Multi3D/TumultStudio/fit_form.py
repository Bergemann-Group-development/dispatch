from .auxfuncs import *
from .ServerManager import *
from .synth import run_synth, StudioSynth, StudioCoG
import ipyformkit as ifk
import numpy as np
import copy
import tumult
from tumult.auxfuncs import atomic_weights, molecule_dict
from tumult.linelist import LineList

element_ids = ('None',) + tuple(atomic_weights.keys())
fields = {('fitting region', 'element'): (('line mask', 'selection'), element_ids),
          ('abnd', 'vshift', 'vdop'): (0.0, 0.0, 0.0),
          ('vrot', 'fac', 'tilt'): (0.0, 1.0, 0.0),
          ('fit abnd', 'fit vshift', 'fit vdop'): (False, True, False),
          ('fit vrot', 'fit fac', 'fit tilt'): (False, False, False),
          'apply': None}

disable = {'vshift': lambda d: d['fit vshift']==True,
           'vdop': lambda d: d['fit vdop']==True,
           'vrot': lambda d: d['fit vrot']==True,
           'fac': lambda d: d['fit fac']==True,
           'tilt': lambda d: d['fit tilt']==True,
           'abnd': lambda d: d['fit abnd']==True or d['element']=='None',
          }

formatters = {'abnd': lambda x: f'{x:.2f}',
              'vshift': lambda x: f'{x:+.2f}',
              'vdop': lambda x: f'{x:+.2f}',
              'vrot': lambda x: f'{x:+.2f}',
              'fac': lambda x: f'{x:.4f}',
              'tilt': lambda x: f'{x:.4f}'}

fit_form = ifk.Form(fields, disable=disable)
element_wid = fit_form.widgets_dict['element'].wid
abundance_wid = fit_form.widgets_dict['abnd'].wid
fit_button = fit_form.widgets_dict['apply'].wid
    
def init_fit_form(studio):
    @log_callback_errors
    def on_fit_click(btn):
        old_run = studio.browser.selected.object
        i = studio.synth_source.data['run_ids'].index(id(old_run))
        studio.progress.description = f'Fitting {old_run.label}:'
        new_run, xfit, yfit = fit_synth(old_run, studio.synth_folder, studio.fig_source, studio.progress)

        @update_data_source(studio.synth_source)
        def update_synth(data):
            data['xs'][i] = xfit
            data['ys'][i] = yfit
            data['run_ids'][i] = id(new_run)
            data['labels'][i] = new_run.label

            if new_run is not old_run:
                studio.browser.remove_object(old_run)
                studio.browser.add_object(new_run)
                studio.browser.selected = studio.browser.list_items[-1]  # select the new run
                widget_label = studio.browser.selected.children[0].label
                studio.set_synth_controls(new_run, widget_label)

        update_synth()
        studio.progress.value = studio.progress.max
        studio.progress.description = 'Done!'

    fit_button.on_click(on_fit_click)

def set_fit_form_values(values):
    """Set the values of the fit form."""
    val = copy.deepcopy(values)
    for k, v in val.items():
        if k in formatters and type(v) is not str:
            val[k] = formatters[k](v)
    
    fit_form.set_values(val)

def set_fit_form(run):
    if hasattr(run, 'fit_params'):
        set_fit_form_values(run.fit_params)

    llist = run.llist.sort_values(by='strength', ascending=False)
    eids = llist.ion.str.split(expand=True)[0].unique()
    output_ids = []
    for id_ in eids:
        if id_ in molecule_dict:
            meids = molecule_dict[id_]
            for eid in meids:
                if eid not in output_ids:
                    output_ids.append(eid)
        elif id_ not in output_ids:
            output_ids.append(id_)

    element_wid.options = ('None',) + tuple(output_ids)

def fit_synth(run, folder, source, progress=None):
    values = fit_form.check_and_return_values()
    lm = run.line_mask

    # LineMask data is always normalised
    if values['fitting region'] == 'line mask':
        oxx = lm.xdata
        oyy = lm.ydata
    elif values['fitting region'] == 'selection':
        islc = np.array(source.selected.indices)
        oxx = source.data['x'][islc]
        oyy = np.interp(oxx, lm.xdata, lm.ydata)

    fit = []
    bounds = {}
    if values['fit vshift']: fit.append('vshift')
    if values['fit vdop']: fit.append('vdop')
    if values['fit vrot']: fit.append('vrot')
    if values['fit fac']: fit.append('fac')
    if values['fit tilt']: fit.append('tilt')
    if values['fit abnd']: fit.append('abnd')
        
    kwargs = {'fit':fit}
    if not values['fit vshift']: kwargs['vshift']=values['vshift'] 
    if not values['fit vdop']: kwargs['vdop']=values['vdop'] 
    if not values['fit vrot']: kwargs['vrot']=values['vrot']
    if not values['fit fac']: kwargs['fac']=values['fac']
    if not values['fit tilt']: kwargs['tilt']=values['tilt']
    if not values['fit abnd'] and 'abnd' in values:
        kwargs['abnd']=values['abnd']

    LTE = run.mode=='LTE'
    eid = values['element']

    if 'abnd' in fit:
        bounds['abnd'] = (-10.0, 12.0)

    def add_run_to_cog(run, anew):
        off = anew - run.a0
        suffix = f'_{off*100:+04.0f}'
        name = run.cog_name + suffix
        settings = {'composition_params': {'abundances': f'{eid}={anew:.2f}'}}
        new_run = run_synth(folder, name, settings, run.line_mask, run.folder, run.llist, progress=progress)
        if hasattr(run, 'runs'):
            runs = list(run.runs) + [new_run]
        else:
            runs = [run, new_run]
        abnds = [r.abundances[eid] for r in runs]
        cog = tumult.CoG(runs, abnds=abnds)
        cog = StudioCoG.from_run(cog)
        cog.a0 = run.a0
        cog.cog_name = run.cog_name
        p, _ = cog.fit_observed(oxx, oyy, LTE=LTE, bounds=bounds, **kwargs)
        return cog, p
    
    def iterate_fit(cog, p, afit, anew):
        while (abs(afit - anew) > 0.005):
            anew = afit
            cog, p = add_run_to_cog(cog, anew)
            afit = p['abnd']
        return cog, p

    if isinstance(run, StudioCoG):
        cog = run
        if values['fit abnd']:
            p, _ = cog.fit_observed(oxx, oyy, LTE=LTE, bounds=bounds, **kwargs)
            afit = p['abnd']
            anew = find_nearest(cog.abnds, afit)
            cog, p = iterate_fit(cog, p, afit, anew)
        elif 'abnd' in kwargs:
            add_run = beyond_sorted(cog.abnds, kwargs['abnd'])
            add_run = add_run or len(cog.abnds) < 3
            add_run = add_run or np.min(abs(cog.abnds - kwargs['abnd'])) > 0.2
            add_run = add_run and not kwargs['abnd'] in cog.abnds
            if add_run:
                cog, p = add_run_to_cog(cog, kwargs['abnd'])
            else:
                p, _ = cog.fit_observed(oxx, oyy, LTE=LTE, bounds=bounds, **kwargs)
        else:
            kwargs['abnd'] = cog.a0
            p, _ = cog.fit_observed(oxx, oyy, LTE=LTE, bounds=bounds, **kwargs)
                
        run = cog # replace original run by cog
    else:
        if eid == 'None':
            p, _ = run.fit_observed(oxx, oyy, LTE=LTE, bounds=bounds, **kwargs)
        else:
            name = init_cog_folder(run, '_+000')
            a0 = run.abundances[eid]
            run.a0 = a0
            run.cog_name = name
            
            if values['fit abnd']:
                anew = a0+0.1
                cog, p = add_run_to_cog(run, anew)
                cog, p = iterate_fit(cog, p, p['abnd'], anew)
            else:
                anew = kwargs['abnd']
                cog, p = add_run_to_cog(run, anew)

            cog.element_id = eid
            cog.ic = run.ic
            run = cog # replace original run by cog

    del kwargs['fit']
    kwargs.update(p)
    xfit, yfit = run.crop(LTE=LTE, **kwargs)
    run.fit_params = kwargs
    set_fit_form_values(run.fit_params)
    return run, xfit, yfit


def init_cog_folder(run, suffix):
    parent, name = (run.folder.parent, run.folder.name)
    cog_folder = unique_folder_name(parent / f'CoG_{name}')
    logging.info(f'Creating folder {cog_folder}.')
    cog_folder.mkdir()
    new_folder = unique_folder_name(cog_folder / f'{name}{suffix}')
    logging.info(f'Renaming {run.folder} to {new_folder}.')
    run.rename_folder(new_folder)
    # without the suffix and possible unique suffixes
    return cog_folder.name.split(suffix)[0]