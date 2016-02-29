import matplotlib.pyplot as pl
import read_file
reload(read_file)

#pl.set_cmap('gray')
#pl.set_cmap('bone')

i0=10
for i in range(8):
    case='case{:02d}'.format(i+i0)
    name='data/'+case+'_{:03d}_{:04d}'
    file=name.format(1,50)
    f=read_file.read_file(file)
    pl.figure(i+1,figsize=(10,8))
    pl.ion()
    pl.set_cmap('hsv')
    pl.clf()
    pl.imshow(f.d[2:130,2:130,0].squeeze().transpose(),origin='lower')
    pl.title(case+'  nu ='+str(f.nu))
    pl.colorbar()
    pl.draw()
    pl.pause(0.1)