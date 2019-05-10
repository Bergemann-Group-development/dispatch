"""
   Produce plots showing average linear momentum components
   from runs with &turbulence_params .. verbose=1 .. /
"""
import os
import sys
import matplotlib.pyplot as pl

for f in sys.argv[1:]:
  cmd='grep update_gl '+f
  #print(cmd)

  time=[]
  vx=[]
  vy=[]
  vz=[]
  for line in os.popen(cmd,'r'):
    words = line.strip().split()
    time.append(float(words[4]))
    vx.append(float(words[5]))
    vy.append(float(words[6]))
    vz.append(float(words[7]))

  pl.clf()
  pl.plot(time,vx)
  pl.plot(time,vy)
  pl.plot(time,vz)
  fp=f+'.png'
  print(fp)
  pl.savefig(fp)
