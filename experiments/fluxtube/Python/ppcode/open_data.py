__all__ = ['open_data']

def open_data (filename):
  """
  This routine opens the data file, reads the header, and returns
  the open file and the header.

  Example: (f, header) = ppcode.open_data ("Data/fields-000000.dat")
  """

  from fortranfile import FortranFile
  from bunch import Bunch
  from os.path import basename, splitext
  from numpy import byte

  f = FortranFile (filename)

  ids = {}
  fmts = {}

  while True:

    fmt, mid = f.readInts ()
    if fmt < 0: break
    sid = f.readString ()
    path = sid.split (' ')[0]
    module = splitext (basename (path))[0]
    n, = f.readInts ()

    ids[module] = sid
    fmts[module] = fmt

    if module == 'io':
      ptime = read_data (f)
      indices = read_data (f)

    elif module == 'params':
      mpi = read_data (f)
      dims = read_data (f)
      ipar = read_data (f)
      rpar = read_data (f)
      lpar1 = read_data (f)
      lpar2 = read_data (f)
      lpar3 = read_data (f)
      run_id = dims[6] if fmt > 1 else -1
      do_vth = lpar1[21] if fmt > 6 else 0

    elif module == 'species':
      rspec, ispec = read_data (f)
      pname = []
      for i in xrange (ispec.size/18):
        pname += [ispec[14+i*18:18+i*18].tostring ().strip ()]

    elif module == 'grid':
      rgrid, igrid = read_data (f)

    elif module == 'units':
      runits, iunits = read_data (f)

    elif module == 'stat':
      rstat = read_data (f)

    elif module == 'slice':
      rslice, islice = read_data (f)

    elif module == 'synthetic_spectra':

      ifname = read_data (f)
      isynth_spectra = read_data (f)
      rsynth_spectra = read_data (f)
      lsynth_spectra = read_data (f)
      obs_synth_spectra = read_data (f)

    else:
      print 'unknown module: %s' % module
      print sid
      print fmt
      for i in range (n): read_data (f)

  dictionary = {}

  dictionary.update (**{
        "run_id":run_id,
           "ids":Bunch (**ids),
          "fmts":Bunch (**fmts),
         "nodes":mpi[0],
          "rank":mpi[1],
      "periodic":mpi[5:8],
          "time":rpar[0],
            "dt":rpar[1],
           "Cdt":rpar[5],
         "ptime":ptime[0],
       "indices":indices[0],
      "nspecies":ipar[2],
          "mdim":dims[0],
        "mcoord":dims[1],
    "do_lorentz":lpar1[0],
    "do_maxwell":lpar1[8],
        "do_vth":do_vth,
             "n":igrid[0:3],
            "lb":igrid[3:6],
            "ub":igrid[6:9],
            "gs":rgrid[9:12],
             "s":rgrid[12:15],
            "ds":rgrid[15:18],
            "gn":igrid[18:21],
        "offset":igrid[21],
           "rlb":rgrid[22:25],
           "rub":rgrid[25:28],
          "grlb":rgrid[28:31],
          "grub":rgrid[31:34],
            "gV":rgrid[34],
             "V":rgrid[35],
            "dV":rgrid[36],
         #"units":units,
     #"constants":constants,
           #"elm":elm
  })
  if len (pname) > 0:
    dictionary.update (**{
              "tnp":ispec[0::18],
             "imax":ispec[2::18],
             "mass":rspec[4::18],
           "charge":rspec[5::18],
         "lifetime":rspec[6::18],
      "temperature":rspec[7::18],
       "min_weight":ispec[8::18],
               "mp":ispec[9::18],
               "np":ispec[10::18],
             "nmax":ispec[11::18],
             "nmin":ispec[12::18],
          "ntarget":ispec[13::18],
            "pname":pname
  })

  header = Bunch (**dictionary)

  return (f, header)

def read_data (f):
  """
  Utility routine to parse header part of data files
  """
  from numpy import float32

  length, typ = f.readInts ()
  if typ == 1:
    ints = f.readInts ()
    assert ints.size == length
    return ints
  if typ == 2:
    reals = f.readReals (prec = 'f')
    assert reals.size == length
    return reals
  if typ == 3:
    reals = f.readReals (prec = 'f')
    assert reals.size == length
    ints = f.readInts ()
    assert ints.size == length
    return (reals, ints)
