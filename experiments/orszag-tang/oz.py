import matplotlib.pyplot as pl
import read_file
reload(read_file)

pl.figure(figsize=(10,8))
pl.ion()
pl.set_cmap('hsv')
#pl.set_cmap('gray')
#pl.set_cmap('bone')

case = 'case3'
name='data/'+case+'_{:03d}_{:04d}'
nt=51
for it in range(51):
    file=name.format(1,it)
    f=read_file.read_file(file)
    pl.clf()
    pl.imshow(f.d[2:130,2:130,0].squeeze().transpose(),origin='lower')
    pl.title(case+'  nu ='+str(f.nu))
    pl.colorbar()
    pl.draw()
    pl.pause(0.1)