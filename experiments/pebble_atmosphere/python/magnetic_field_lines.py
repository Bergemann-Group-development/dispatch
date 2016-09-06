""" Based on a Mayavi example
"""
import numpy as np
from mayavi import mlab
from scipy import special
from numpy import fromfile, float32, prod, size

# Read the data (you may want to have a for loop here)
n=64
shape=(n,n,n)
fd=open('snapshot.dat','rb')
for i in range(5):
    dd=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    ds=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    Uz=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    Uy=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    Ux=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    Bz=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    By=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    Bx=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)

# Find the max
Bmax=np.sqrt(np.max(Bx**2+By**2+Bz**2))

# Open a figure
fig = mlab.figure(1, size=(800,800), bgcolor=(1,1,1), fgcolor=(0,0,0))

# Specify the components
field = mlab.pipeline.vector_field(Bx,By,Bz)
magnitude = mlab.pipeline.extract_vector_norm(field)

# Create contours
contours = mlab.pipeline.iso_surface(magnitude,
                                    contours=[0.7*Bmax,0.95*Bmax, ],
                                    transparent=True,
                                    opacity=0.4,
                                    colormap='YlGnBu',
                                    vmin=0, vmax=2)
                                    
# Create field lines, with interactive widget
field_lines = mlab.pipeline.streamline(magnitude, seedtype='line',
                                       integration_direction='both',
                                       vmin=0, vmax=1)
# Tweak the streamline.
field_lines.stream_tracer.maximum_propagation = 100.
field_lines.seed.widget.point1 = [n/2, n/2-2, n/2]
field_lines.seed.widget.point2 = [n/2, n/2+2, n/2]
field_lines.seed.widget.resolution = 20
field_lines.seed.widget.enabled = True

mlab.view()
mlab.show()
