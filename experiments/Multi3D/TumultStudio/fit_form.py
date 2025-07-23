from .auxfuncs import *
from .ServerManager import *
import ipyformkit as ifk
import numpy as np
import m3dis

fields = {'fitting region':('line mask', 'selection'),
          ('vshift', 'vdop', 'vrot'): (0.0,0.0,0.0),
          ('fac', 'tilt'): (1.0, 0.0),
          ('fit vshift', 'fit vdop', 'fit vrot'): (True, False, False),
          ('fit fac', 'fit tilt'): (False, False),
          'fit': None}

disable = {'vshift': lambda d: d['fit vshift']==True,
           'vdop': lambda d: d['fit vdop']==True,
           'vrot': lambda d: d['fit vrot']==True,
           'fac': lambda d: d['fit fac']==True,
           'tilt': lambda d: d['fit tilt']==True,
          }

fit_form = ifk.Form(fields, disable=disable)

def fit_synth(run, source):
    values = fit_form.check_and_return_values()
    lm = run.line_mask

    if values['fitting region'] == 'line mask':
        oxx = lm.xdata
        oyy = lm.ydata
    elif values['fitting region'] == 'selection':
        islc = np.array(source.selected.indices)
        oxx = source.data['x'][islc]
        oyy = source.data['y'][islc]

    fit = []
    if values['fit vshift']: fit.append('vshift')
    if values['fit vdop']: fit.append('vdop')
    if values['fit vrot']: fit.append('vrot')
    if values['fit fac']: fit.append('fac')
    if values['fit tilt']: fit.append('tilt')
        
    kwargs = {'fit':fit}
    if not values['fit vshift']: kwargs['vshift']=values['vshift'] 
    if not values['fit vdop']: kwargs['vdop']=values['vdop'] 
    if not values['fit vrot']: kwargs['vrot']=values['vrot']
    if not values['fit fac']: kwargs['fac']=values['fac']
    if not values['fit tilt']: kwargs['tilt']=values['tilt']
    

    LTE = run.mode == 'LTE'
    p, _ = run.fit_observed(oxx, oyy, LTE=LTE, **kwargs)
    del kwargs['fit']
    kwargs.update(p)
    xfit, yfit = run.crop(LTE=LTE, **kwargs)
    run.fit_params = kwargs
    fit_form.set_values(p)
    
    return xfit, yfit