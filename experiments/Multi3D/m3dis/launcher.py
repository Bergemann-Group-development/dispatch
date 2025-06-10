import ipywidgets as widgets
from IPython.display import display, HTML, clear_output
import ipyformkit as ifk
import os, copy
from .launch import launch
from .m3dis import m3dis

params = {}
hides = {}
disables = {}
checks = {}

params['atmos_params'] = {
    ('atmos_format', 'atmos_file'): (('Marcs', 'Stagger', 'Stagger2', 'Co5bold', 'Text'),
                                    'input_multi3d/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod'),
    ('vmic', 'FeH'): (-1.0, 0.0),
    ('nx', 'ny', 'nz'): (1, 1, 128),
    'advanced':{
      'dims':8,
      ('use_rho', 'use_ne'): (True, True),
    }
}
hides['atmos_params'] = {
  'nx': lambda d: d['atmos_format'] in ['Marcs', 'Text'],
  'ny': lambda d: d['atmos_format'] in ['Marcs', 'Text'],
}
checks['atmos_params'] = {
  'nx': lambda d: d['nx'] > -1,
  'ny': lambda d: d['ny'] > -1,
  'nz': lambda d: d['nz'] > 0,
}

params['atom_params'] = {
    'atom_file': './input_multi3d/atoms/atom.ba06',
    ('abundance', 'relative_abundance'): (2.11, False),
    'advanced': {
      ('convlim', 'convmax'): (1e-2, 1e-3),
    }
}

fixed = ('disk_center', 'set_a2', 'set_a4', 'set_a6', 'set_a8', 'set_b4', 'set_b6', 'set_b8')

params['m3d_params'] = {
    ('maxiter', 'ng_step'): (99, -1),
    'advanced': {
      ('m1d_legacy_mode','decouple_continuum'): (False, True),
      ('rotate_atmos', 'rotate_continuum'): (True, True),
      ('short_scheme', 'short_ntheta', 'short_nphi'): (('radau', 'lobatto', 'gauss', *fixed), 2, 4),
      ('long_scheme', 'long_ntheta', 'long_nphi', 'custom_mu'): (('radau', 'lobatto', 'gauss', 'custom', *fixed), 4, 4, '1.0 0.8 0.6 0.4 0.2'),
    }
}

hides['m3d_params'] = {
  'long_ntheta': lambda d: d['long_scheme'] in ['custom', *fixed],
  'long_nphi': lambda d: d['long_scheme'] in ['custom', *fixed],
  'short_ntheta': lambda d: d['short_scheme'] in fixed,
  'short_nphi': lambda d: d['short_scheme'] in fixed,
  'custom_mu': lambda d: d['long_scheme'] != 'custom',
}

params['composition_params'] = {
    'abundance':'C=8.52, O=8.75, Ni=6.22',
    'absmet_file': './input_multi3d/absmet',
    'advanced': {
      'abund_file': './input_multi3d/abund/abund_magg',
    }
}

params['spectrum_params'] = {
    ('aa_blue', 'aa_red', 'R'): (4200, 4300, 1e5),
    'lam_file': './file_with_wavelengths.txt',
}
disables['spectrum_params'] = {
    'aa_blue': lambda d: d['lam_file'] != '',
    'aa_red': lambda d: d['lam_file'] != '',
    'R': lambda d: d['lam_file'] != '',
}

launch_params = {
    'input_file': 'ba_test1', 
    ('overwrite', 'verbose'): (False, False),
    'Launch': None,
}

#=====================================================================
class launcher():
  #=====================================================================
  def __init__(self):
    forms = []
    for key, p in params.items():
        form = ifk.Form(p, title=key,
                        hide=hides.get(key, None),
                        check=checks.get(key, None),
                        disable=disables.get(key, None),
                        )
        forms.append(form)

    launch_form = ifk.Form(launch_params)
    launch_button = launch_form.widgets_dict['Launch'].wid
    forms.append(launch_form)

    self.gui = ifk.Masonry(forms)
    self.gui.display()

    def _on_launch_click(btn):
        """Launch the M3DIS simulation."""
        # Get the values from the widgets
        nml = {}
        for form in self.gui.forms[:-1]:
          settings = form.check_and_return_values()
          if settings != {}:
            nml[key] = settings

        lp = launch_form.check_and_return_values()

        with output:
          clear_output(wait=True)
          name = lp['input_file']
          
          v = 2 if lp['verbose'] else 1
          o = lp['overwrite']
          # Launch M3DIS with the parameters
          if name:
            output_folder, ok = launch(name, nml_settings=nml, overwrite=o, verbose=v)
            if ok:
                if os.sep in name:
                  name = name.split(os.sep)[-1]

                print('reading ' + output_folder)
                self.__dict__[name] = m3dis.read(output_folder)
          else:
            print("Please enter a name for the input file.")

    launch_button.on_click(_on_launch_click)

    output = widgets.Output()
    output.add_class('ifk-widget-label')
    display(output)
