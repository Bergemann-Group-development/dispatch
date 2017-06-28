import numpy as np
import image
reload(image)

class task():
    dt=0.0
    def __init__(self,id=0,llc=[0,0],size=[1,1]):
        self.id=id
        self.im=image.image()
        self.llc=llc
        self.size=size
        self.extent=[llc[0],llc[0]+size[0],llc[1],llc[1]+size[1]]
        self.dt=1.0*(np.random.uniform()-1.0)-0.4
        self.dt=10.**self.dt
        self.t=0.0
    def update(self):
        self.t=self.t+self.dt
    def show(self,color=0.7):
        #print self.extent
        self.im.rectangle(self.extent,color)
    def info(self):
        print 'id:',self.id,' llc:',self.llc,' t:',self.t,self.dt
