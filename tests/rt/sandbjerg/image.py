import matplotlib.pyplot as pl
import numpy as np

class image():
    def __init__(self,img=[[1,1],[1,1]],size=[20,10],\
                 cmap='Spectral',color=0.5):
        ''' Create an image object '''
        self.extent=[0,size[0],0,size[1]]
        self.size=size
        self.cmap=cmap
        self.img=np.array(img)
        self.color=color
        self.vmin=0.0
        self.vmax=1.0
        pl.xlim(0,self.size[0])
        pl.ylim(0,self.size[1])
        self.show()
    def show(self):
        ''' Show a full size background image '''
        pl.imshow(self.img*self.color,extent=self.extent,\
            vmin=self.vmin,vmax=self.vmax,cmap=self.cmap)        
    def rectangle(self,extent=[0,1,0,1],color=0.5):
        ''' Show a rectangle '''
        pl.imshow(self.img*color,extent=extent,\
            vmin=self.vmin,vmax=self.vmax,cmap=self.cmap)
        
        