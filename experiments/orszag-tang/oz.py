import matplotlib.pyplot as pl
import read_file
reload(read_file)

pl.figure(figsize=(10,8))
pl.ion()
pl.set_cmap('hsv')
#pl.set_cmap('gray')
title='nu = 0.001, 0.0, 1.0, 0.2, 0.1  t ={:5.2f}'

name='data/orszag-tang_{:03d}_{:04d}'
nt=51
for it in range(51):
    file=name.format(1,it)
    f=read_file.read_file(file)
    pl.clf()
    pl.imshow(f.d[2:130,2:130,0].squeeze().transpose(),origin='lower')
    pl.title(title.format(0.1*it))
    pl.colorbar()
    pl.draw()
    pl.pause(0.1)