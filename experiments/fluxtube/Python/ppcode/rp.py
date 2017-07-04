# -*- coding: utf-8 -*-
__all__ = ['rp']

def rp (filename):
  """
  Reads particles from PPcode data files.

  Usage: (p, h) = ppcode.rp ("Data/particles-000000.dat")
  """

  from numpy import int16, int64, float32, recarray, squeeze
  from ppcode import open_data
  (f, h) = open_data (filename)

  # Number of particles written out
  np = f.readInts ('l')

  # Marker size is hard coded to 8
  f.HEADER_PREC = 'l'

  # Define structured array
  particle = recarray ([h.nspecies, np[0]], dtype = [
    ('x', float32),
    ('y', float32),
    ('z', float32),
    ('r', float32, (h.mdim,)),
    ('q', int16,   (h.mdim,)),
    ('v', float32, (h.mdim,)),
    ('e', float32),
    ('w', float32),
    ('i', int64)
    ])

  for isp in xrange (h.nspecies):

    for idim in xrange (h.mdim): particle.r[isp,:,idim] = f.readReals ('f')
    for idim in xrange (h.mdim): particle.q[isp,:,idim] = f.readInts ('h')
    for idim in xrange (h.mdim): particle.v[isp,:,idim] = f.readReals ('f')

    particle.e[isp, :] = f.readReals ('f')
    particle.w[isp, :] = f.readReals ('f')

    print "VALUES", h.ds[0], particle.q[...,0], particle.r[...,0], h.rlb[0]
    particle.x = h.ds[0]*(particle.q[...,0] + particle.r[...,0]) + h.rlb[0]
    particle.y = h.ds[1]*(particle.q[...,1] + particle.r[...,1]) + h.rlb[1]
    particle.z = h.ds[2]*(particle.q[...,2] + particle.r[...,2]) + h.rlb[2]

    if h.indices: particle.i[isp, :] = f.readInts ('l')

  # Make sure end of file is reached
  assert not f.read ()

  f.close ()

  return (squeeze (particle), h)
