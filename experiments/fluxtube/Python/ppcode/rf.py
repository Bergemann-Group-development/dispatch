__all__ = ['rf']

def rf (filename):
  """
  Reads fields from PPcode data files.

  Usage: (f, h) = ppcode.rf ("Data/fields-000000.dat")
  """

  from numpy import float32, empty, all as np_all
  from ppcode import open_data
  from bunch import Bunch

  (f, h) = open_data (filename)

  do_em = h.do_lorentz or h.do_maxwell

  if h.fmts.params > 4:
    tmp = f.readInts ()
    if h.fmts.io >= 8:
      assert tmp[0] == do_em
      assert np_all (tmp[1:4] == h.lb)
      assert np_all (tmp[4:7] == h.ub)
    else:
      assert np_all (tmp[0:3] == h.lb)
      assert np_all (tmp[3:6] == h.ub)
    del tmp

  ngb = h.gn

  if do_em:
    Bx = f.readReals ('f').reshape (ngb).squeeze ()
    By = f.readReals ('f').reshape (ngb).squeeze ()
    Bz = f.readReals ('f').reshape (ngb).squeeze ()

    Ex = f.readReals ('f').reshape (ngb).squeeze ()
    Ey = f.readReals ('f').reshape (ngb).squeeze ()
    Ez = f.readReals ('f').reshape (ngb).squeeze ()

  nsp = h.nspecies

  rho = empty ([nsp] + list (ngb))
  rux = empty ([nsp] + list (ngb))
  ruy = empty ([nsp] + list (ngb))
  ruz = empty ([nsp] + list (ngb))

  for isp in range (nsp):
    rho[isp, ...] = f.readReals ('f').reshape (ngb)
    rux[isp, ...] = f.readReals ('f').reshape (ngb)
    ruy[isp, ...] = f.readReals ('f').reshape (ngb)
    ruz[isp, ...] = f.readReals ('f').reshape (ngb)

  # Make sure end of file is reached
  assert not f.read ()

  f.close ()

  rho = rho.squeeze ()
  rux = rux.squeeze ()
  ruy = ruy.squeeze ()
  ruz = ruz.squeeze ()

  varnames = ["Bx","By","Bz","Ex","Ey","Ez","rho","rux","ruy","ruz"]

  return (Bunch ([(name, eval (name)) for name in varnames]), h)
