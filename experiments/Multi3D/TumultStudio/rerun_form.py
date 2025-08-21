from .auxfuncs import *
from .ServerManager import *
import ipyformkit as ifk
from tumult.auxfuncs import update_dict
from .llist_annotations import get_llist_table

ion_ids = np.array(['I', 'II', 'III', 'IV', 'V ', 'VI', 'VII', 'VIII'])

fields = {('rerun', 'x-match levels'): (None, None)}
rerun_form = ifk.Form(fields)
rerun_form.vbox.add_class('rerun-form')
rerun_button = rerun_form.widgets_dict['rerun'].wid
xmatch_button = rerun_form.widgets_dict['x-match levels'].wid

def init_rerun_form(studio):
    def rerun_on_click(btn):
        run = studio.current
        studio.add_object('Synth', studio.plot_synth, rerun_synth, run)
        
    rerun_button.on_click(rerun_on_click)

    def xmatch_on_click(btn):
        run = studio.current
        xmatch_levels(studio, run)
    xmatch_button.on_click(xmatch_on_click)

def set_rerun_form(run):
  xmatch_button.disabled = not hasattr(run.atmos, 'atoms') or len(run.atmos.atoms) == 0

def rerun_synth(studio, name, old_run, **kwargs):
  settings = {}
  
  if hasattr(old_run, 'eid'):
    update_dict(settings, {'composition_params':
                           {'abundances': f'{old_run.eid}={old_run.fit_params["abnd"]:.2f}'}})
  
  run = old_run.rerun(studio, name, settings, **kwargs)
  return run
  
def xmatch_levels(studio, run):
  df1 = run.llist.copy().reset_index(drop=True)
  df1['ion'] = df1.ion.str.strip()
  df1['idx'] = np.arange(len(df1))

  for atom in run.atmos.atoms:
    bbt = atom.bbt
    bbt['evi'] = atom.level.ev[bbt.i].values
    bbt['ion'] = atom.id + ' ' + ion_ids[atom.level.ion[bbt.i]-1]
    df2 = bbt[['ion', 'lam0', 'evi']].copy()
    df2['i'] = (bbt['i'] + 1).astype(int)
    df2['j'] = (bbt['j'] + 1).astype(int)

    # Merge all combinations with matching id
    merged = df1.merge(df2, on='ion', suffixes=('_1', '_2'))
    merged = merged[np.abs(merged['evi_1'] - merged['evi_2']) <= 0.01]

    # Find the row in df1 with minimal |x_1 - x_2|
    idx = merged.groupby(['ion', 'evi_2', 'lam0_2'])['lam0_1'].apply(lambda g: (np.abs(g - merged.loc[g.index, 'lam0_2'])).idxmin())
    nearest = merged.loc[idx]
    df1.loc[nearest.idx, ['i', 'j']] = nearest[['i_2', 'j_2']].values

  df1.drop(columns=['idx'], inplace=True)
  row = run.llist_html.highlighted_row
  run.llist = df1
  run.llist_html = get_llist_table(df1)
  instructions = run.chosen(studio)
  studio.right_sidebar.children = instructions['tabs']
  run.llist_html.highlight_row(row)
