import copy as copy
import numpy as np
import pylab

class point_on_sphere():
    def __init__(self,dtheta,r=1.0,verbose=0):
        self.verbose=0
        self.dtheta=dtheta
        self.set_cartesian(1.0,0.0,0.0)
        self.id=1
        self.nr=[0.0,0.0,1.0]
    def cartesian(self,theta,phi,r):
        x=r*np.cos(theta)*np.cos(phi)
        y=r*np.cos(theta)*np.sin(phi)
        z=r*np.sin(theta)
        return x,y,z
    def to_cartesian(self):
        self.x,self.y,self.z=self.cartesian(self.theta,self.phi,self.r)
    def spherical(self,x,y,z):
        r=np.sqrt(x**2+y**2+z**2)
        theta=np.arcsin(z/r)
        phi=np.arctan2(y,x)
        return theta,phi,r
    def to_spherical(self):
        self.theta,self.phi,self.r=self.spherical(self.x,self.y,self.z)
    def set_cartesian(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        self.to_spherical()
    def set_spherical(self,theta,phi,r):
        self.theta=theta
        self.phi=phi
        self.r=r
        self.to_cartesian()
        
    def go_east(self,dtheta):
        #self.size=self.r*dtheta
        #self.theta-=dtheta/2
        self.phi+=dtheta/np.cos(self.theta)
        #self.theta+=dtheta/2
        self.to_cartesian()

    def go_west(self, dtheta):
        #self.size=self.r*dtheta
        #self.theta-=dtheta/2
        self.phi -= dtheta/np.cos(self.theta)
        #self.theta+=dtheta/2
        self.to_cartesian()
    
    def go_north(self,dtheta):
        self.theta+=dtheta
        self.to_cartesian()

    def go_south(self, dtheta):
        self.theta -= dtheta
        self.to_cartesian()
        
    def east_west(self):
        self.phi=-self.phi
        self.to_cartesian()
    def north_south(self):
        self.theta=-self.theta
        self.to_cartesian()

    def get_corner_points(self):
        corners=[]
        self.to_cartesian()
        pos=self.xyz()
        self.unit_vectors()
        c=np.zeros(2)
        c[0]=-0.5*self.size
        c[1]=+0.5*self.size
        for i in (0,1):
            for j in (0,1):
                for k in (0,1):
                    new = copy.deepcopy(self)
                    new.x,new.y,new.z=np.matmul(self.e,[c[i],c[j],c[k]])+pos
                    new.to_spherical()
                    corners.append(new)


        return corners
        
    def permute(self):
        x = self.x
        self.x=self.y
        self.y=self.z
        self.z=x
        n1=self.nr[0]
        self.nr[0]=self.nr[1]
        self.nr[1]=self.nr[2]
        self.nr[2]=n1
        self.to_spherical()
    def reflect(self):
        self.x=-self.x
        self.y=-self.y
        self.z=-self.z
        self.to_spherical()
        
    def info(self,label,verbose=0):
        if self.verbose>0 or verbose>0:
            print(label)
            print(' x:{:5.3f}  theta:{:5.3f}'.format(self.x,self.theta))
            print(' y:{:5.3f}    phi:{:5.3f}'.format(self.y,self.phi))
            print(' z:{:5.3f}      r:{:5.3f}'.format(self.z,self.r))
    def is_inside(self,theta0,dtheta):
        """ For a patch to be inside a volleyball face its central point 
            must have a latitude less that the limit, and its latitude
            with respect to the east face must be larger than the limit.
            A more careful constraint is that it is outside if *all* of
            its corners are outside the limit. The last corner to leave
            is has a position displaced by size/2 to the west, nord, and
            down.
        """
        inside=self.theta<theta0
        # move point to south-west lower corner before converting to east VB coords
        new=copy.deepcopy(self)
        corners=new.get_corner_points()
        inside=np.zeros(8)
        inside2=np.zeros(8)
        i=0
        for c in corners:
            inside[i]=c.theta<=theta0
            #print (i,c.theta, theta0, c.xyz(), self.xyz())
            i+=1
        i = 0
        for c in corners:
            c.permute()
            inside2[i] = c.theta >= theta0
            #print('another',i, c.theta, theta0, c.xyz(), self.xyz())
            i += 1
        #new.theta = new.theta - dtheta/2
        #new.phi   = new.phi   - dtheta/(2*np.cos(new.theta))
        #new.r     = new.r*(1-dtheta/2)
        #new.permute()
        return self.theta<=theta0 and any(inside) and any(inside2)
    #
    def distance_to(self,other):
        return other.xyz()-self.xyz()
    #
    def unit_vectors(self):
        """ Define the unit vectors of a coordinate system where
            x1 is East, x2 is North, and x3 is Up
        """
        pos=self.xyz()
        self.e=np.zeros((3,3))
        i3=np.argmax(abs(np.array(pos)))          # local up
        i1=(i3+1)%3                         # local east
        i2=(i1+1)%3                         # local north
        i1=0;i2=1;i3=2
        p=[pos[i1],pos[i2],pos[i3]]
        self.e[:,i3]=p/np.linalg.norm(p)                # unit vector up
        self.e[:,i2]=self.nr#[0.0,1.0,0.0]          # north proxy
        e2=np.cross(self.e[:,i2],self.e[:,i3])
        self.e[:,i1]=e2/np.linalg.norm(e2) # unit vector east
        e3=np.cross(self.e[:,i3],self.e[:,i1])
        self.e[:,i2]=e3/np.linalg.norm(e3) # unit vector north
        #
        #theta=self.theta
        #phi  =self.phi
        #self.e[:,0]=[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
        #self.e[:,1]=[np.cos(theta)*np.cos(phi),np.cos(theta)*np.cos(phi),-np.sin(theta)]
        #self.e[:,2]=[-np.sin(phi),np.cos(phi),0]

    def xyz(self):
        return self.x,self.y,self.z
