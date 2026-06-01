from .auxfuncs import *
from .ServerManager import *
from .threading import execute_in_thread
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
          'NLTE elements': ('No departures in atmos!',),
          'apply': void}

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

fit_form = ifk.Form(fields, disable=disable, formatters=formatters)
defaults = fit_form.get_values()
element_wid = fit_form.inputs_dict['element']
nlte_wid = fit_form.inputs_dict['NLTE elements']
fit_button = fit_form.inputs_dict['apply']
    
def init_fit_form(studio):
    @log_callback_errors
    def on_fit_click(btn):
        item = studio.browser.selected
        old_run = item.obj
        studio.progress.description = f'Fitting {old_run.label}:'

        @execute_in_thread(spinner=item.icon)
        def async_callback(execution):
            new_run, xfit, yfit = fit_synth(old_run, studio.synth_folder, studio.fig_source, runner=execution)

            @update_source_data(studio.synth_source)
            def update_synth(data):
                i = data['run_ids'].index(id(old_run))
                data['xs'][i] = xfit
                data['ys'][i] = yfit
                data['run_ids'][i] = id(new_run)
                data['labels'][i] = new_run.label

                if new_run is not old_run:
                    item.replace_object(new_run)
                    #studio.set_synth_controls(new_run)

            update_synth()
            studio.progress.description = 'Done!'

        async_callback()

    fit_button.on_click(on_fit_click)

    @log_callback_errors
    def on_element_change(change):
        eid = element_wid.value
        if eid == 'None':
            fit_form['abnd'] = 0.0
        else:
            run = studio.current
            if hasattr(run, 'eid') and eid == run.eid:
                abnd = run.fit_params['abnd']
            else:
                abnd = run.abundances[eid]
            fit_form['abnd'] = abnd

    element_wid.observe(on_element_change, names='value')

def set_fit_form(run):
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

    atomids = [atom.id for atom in run.atmos.atoms if atom.id in eids]

    if len(atomids) == 0:
        nlte_wid.options = ('No departures in atmos!',)
    else:
        value = ', '.join(atomids)
        nlte_wid.options = ('LTE it is!', value,)
    nlte_wid.value = nlte_wid.options[-1]

    if hasattr(run, 'fit_params'):
        fit_form.set_values(run.fit_params)

    if hasattr(run, 'eid'):
        element_wid.value = run.eid
    else:
        element_wid.value = 'None'


def fit_synth(run, folder, source, runner=None):
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
    if not values['fit abnd'] and 'abnd' in values and 'eid' != 'None':
        kwargs['abnd']=values['abnd']

    kwargs['LTE'] = 'LTE' in nlte_wid.value or 'No departures' in nlte_wid.value
    eid = values['element']

    if isinstance(run, tumult.CoG):
        if eid == 'None':
            kwargs['abnd'] = run.a0
            p, _ = run.fit_observed(oxx, oyy, **kwargs)  
        else:
            cog = run
            if values['fit abnd']:
                p = cog.iterate_fit(oxx, oyy, runner=runner, **kwargs)
            elif 'abnd' in kwargs:
                if beyond_sorted(cog.abnds, kwargs['abnd']) or len(cog.abnds) < 3:
                    cog.add_abnd(kwargs['abnd'], tol=0.2, runner=runner)
                p, _ = cog.fit_observed(oxx, oyy, **kwargs)
                
            run = cog # replace original run by cog
    else:
        if eid == 'None':
            p, _ = run.fit_observed(oxx, oyy, **kwargs)
        else:
            cog = run.convert_to_CoG(eid)
            
            if values['fit abnd']:
                anew = cog.a0+0.1
                cog.add_abnd(anew, runner=runner)
                p = cog.iterate_fit(oxx, oyy, runner=runner, **kwargs)
            else:
                anew = values['abnd']
                cog.add_abnd(anew, runner=runner)
                p, _ = cog.fit_observed(oxx, oyy, **kwargs)

            run = cog # replace original run by cog

    del kwargs['fit']
    kwargs.update(p)
    xfit, yfit = run.crop(**kwargs)
    run.fit_params = {k:v for k, v in kwargs.items()}
    del kwargs['LTE']
    fit_form.set_values(kwargs)
    return run, xfit, yfit

