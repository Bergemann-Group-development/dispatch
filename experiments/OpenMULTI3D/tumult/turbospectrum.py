from .spectrum import Spectrum
from .auxfuncs import atomic_weights, near_equal_split, tmp_txtfile, BK
from .atmos import AtmosClass
from astropy import units as u
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import subprocess
from pathlib import Path
import multiprocessing as mp
import tempfile
from contextlib import ExitStack
from functools import partial
from io import StringIO
import pickle

TS_PATH = '/mnt/beegfs/gemini/groups/bergemann/users/hoppe/TSFitPy/turbospectrum/'

class turbospectrum(object):
  @property
  #=================================================================
  def atmos(self):
    if not hasattr(self, '_atmos'):
      self._atmos = self.read_atmos()
    return self._atmos
  
  @classmethod
  #=================================================================
  def read(cls, filename):
    self = cls()
    self.filename = filename

    with open(filename, 'r') as f:
      line = f.readline()
      self.mus = [float(s) for s in line.strip('# mu-points').split()]
      df = pd.read_csv(f, sep=r'\s+', header=None)
      self.df = df.apply(pd.to_numeric, errors='coerce')
      self.lam = self.df.iloc[:,0].values

    return self
  
  @classmethod
  #=================================================================
  def launch(cls, input, llist=None, cwd=None, settings=None, verbose=False):
    cwd = Path(cwd or TS_PATH)

    with ExitStack() as stack:
      if llist is not None:
        llist_file = tmp_txtfile(stack)
        llist_list = tmp_txtfile(stack)
        mask = (llist.lam0 >= input['LAMBDA_MIN']) & (llist.lam0 <= input['LAMBDA_MAX'])
        llist[mask].write_to_ts_format(llist_file)
        with open(llist_list, 'w') as f:
          f.write(str(llist_file))
        input['LIST_OF_LINELISTS'] = str(llist_list)

      if settings is not None:
        input.update(settings)

      babsma_in = '\n'.join([f"{pad_string(k)} {pad_string(v)}" for k,v in input.items()])
      babsma_stdout, babsma_stderr = run_process('exec/babsma_lu', babsma_in, cwd=cwd, verbose=verbose)
      bsyn_stdout, bsyn_stderr = run_process('exec/bsyn_lu', babsma_in, cwd=cwd, verbose=verbose)
      self = cls.read(cwd / input['RESULTFILE'])
      self.input = input
      self.llist = llist
      self.stdout = '\n'.join([babsma_stdout, bsyn_stdout])
      self.stderr = '\n'.join([babsma_stderr, bsyn_stderr])
      return self
  
  @classmethod
  #=================================================================
  def concatenate(cls, runs):
    self = cls()
    self.runs = runs
    self.df = pd.concat([run.df for run in runs], axis=0)
    self.lam = self.df.iloc[:,0].values
    return self
  
  @classmethod
  #=================================================================
  def from_tumult_run(cls, run, nthreads=0, **kwargs):
    nml = run.nml

    protons = {i+1:k for i,k in enumerate(atomic_weights.keys())}
    abnds = {k:run.abundances[v] for k,v in protons.items() if v in run.abundances.keys()}

    vmic = nml['atmos_params']['vmic']
    vmic = vmic if vmic > 0 else np.mean(run.atmos.vmic)

    llist = getattr(run, 'llist', None)

    lmin  = nml['spectrum_params']['aa_blue']
    lmax  = nml['spectrum_params']['aa_red']
    lstep = nml['spectrum_params']['daa']
    lpad  = nml['linelist_params']['dlam']

    input = {
        'PURE-LTE': nml['m3d_params']['ilambd'] == 0,
        'LAMBDA_MIN': lmin - lpad,
        'LAMBDA_MAX': lmax + lpad,
        'LAMBDA_STEP': lstep,
        'LOGFILE': 'tmp.log',
        'MODELINPUT': nml['atmos_params']['atmos_file'],
        'MARCS-FILE': True,
        'MODELOPAC': 'tmp.opac',
        'RESULTFILE': 'tmp.spec',
        'INTENSITY/FLUX': 'INTENSITY',
        'SPHERICAL':'False\n30\n300.00\n15\n1.30',      
        #'LIST_OF_LINELISTS': 'linelists_used.txt',
        'XIFIX':f'True\n{vmic:.1f}'
    }
    input['INDIVIDUAL ABUNDANCES'] = len(abnds)
    input.update(abnds)

    # bind kwargs to cls.launch method
    launch = partial(cls.launch, **kwargs)

    if nthreads == 0:
      ts = launch(input, llist)
      lam = ts.df.iloc[:,0].values
      mask = (lam >= lmin) & (lam <= lmax)
      ts.df = ts.df[mask].reset_index(drop=True)
      ts.lam = ts.df.iloc[:,0].values
      return ts
    
    else:
      splits = np.concatenate([[0], np.cumsum(near_equal_split(run.lam.size, nthreads))])
      splits[-1] = run.lam.size-1
      splits = np.round(run.lam[splits], 2)

      args = []
      for i, (start, stop) in enumerate(zip(splits[:-1], splits[1:])):
        input_chunk = input.copy()
        input_chunk['LAMBDA_MIN'] = start - lpad
        input_chunk['LAMBDA_MAX'] = stop + lpad
        input_chunk['LOGFILE'] = f'tmp_{i:03d}.log'
        input_chunk['MODELOPAC'] = f'tmp_{i:03d}.opac'
        input_chunk['RESULTFILE'] = f'tmp_{i:03d}.spec'
        args.append(input_chunk)

      with mp.Pool(nthreads) as pool:
          results = pool.starmap(launch, [(arg, llist) for arg in args])

      for result in results:
          lmin = result.lam[0] + lpad
          lmax = result.lam[-1] - lpad
          lam = result.df.iloc[:,0].values
          mask = (lam >= lmin) & (lam <= lmax)
          result.df = result.df[mask].reset_index(drop=True)
          result.lam = result.df.iloc[:,0].values

      return cls.concatenate(results)
    
  def get_spectrum(self, imu=None, ang=None, norm=True):
    imu = imu or ang
    idx = 1 if norm else 2
    if imu is not None:
      idx = idx + 2*(imu+1)

    xx = self.df.iloc[:,0].values
    yy = self.df.iloc[:,idx].values
    yunit = u.dimensionless_unscaled if norm else u.erg/(u.s*u.cm**2*u.AA)
    return Spectrum(xx, yy, xunit=u.AA, yunit=yunit)
    
  #=================================================================
  def crop(self, imu=None, ang=None, cont=False, norm=True, **kwargs):
    if cont:
      nspec = self.get_spectrum(imu=imu, ang=ang, norm=True)
      uspec = self.get_spectrum(imu=imu, ang=ang, norm=False)
      cspec = uspec / nspec
      spec = cspec / cspec if norm else cspec
    else:
      spec = self.get_spectrum(imu=imu, ang=ang, norm=norm)

    return spec.crop(**kwargs)
  
  #=================================================================
  def read_atmos(self):
    key = 'tau,T,Pe,Pg,ro,P6'
    lines = self.stdout.split('\n')
    data = [l[len(key):] for l in lines if l.startswith(key)]
    string = '\n'.join(data)
    df = pd.read_csv(StringIO(string), sep=r'\s+', header=None, names=['tau','temp','pe','pg','rho','p6'])
    
    atmos = AtmosClass()
    atmos.dim = '1D'
    for col in df.columns:
      setattr(atmos, col, df[col].values)
    
    atmos.ltau = np.log10(atmos.tau)
    atmos.ne = atmos.pe / (BK*atmos.temp)
    return atmos
  
  #=================================================================
  def save(self, filename):
    with open(filename, 'wb') as f:
      pickle.dump(self, f)

  @classmethod
  #=================================================================
  def load(cls, filename):
    with open(filename, 'rb') as f:
      self = pickle.load(f)
    return self



def run_process(executable, input, cwd='./', verbose=True):
    cwd = str(Path(cwd).resolve()) + '/'
    pr = subprocess.Popen(
            cwd + executable,
            cwd=cwd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

    stdout, stderr = pr.communicate(input=input)
    if verbose:
      print(stdout)
      print(stderr)

    return stdout, stderr

def pad_string(s):
  if type(s) == str:
    if any(c in s for c in [' ','/']):
      s = f"'{s}'"
  return s