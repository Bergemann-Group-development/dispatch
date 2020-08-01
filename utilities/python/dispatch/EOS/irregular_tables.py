import numpy as np
import scipy
import scipy.io
import scipy.interpolate

class linear_scale_c(object):
    """ Linear scale for independent variables """
    def __init__(self,n,min,step,name=None):
        self.n=n
        self.min=min
        self.max=min+(n-1)*step
        self.step=step
        self.values=min+step*np.arange(n)
        self.offset=0
        if type(name) is str:
            self.name=name
    def accumulate(self,n,mn,step):
        """ Accumulate index ranges """
        self.step=step
        self.min=min(self.min,mn)
        self.max=max(self.max,mn+(n-1)*step)
        self.n=int((self.max-self.min)/self.step+1.5)
        self.values=self.min+self.step*np.arange(self.n)
    def index(self,x):
        """ Return interpolation index and weight """
        p=(x-self.min)/self.step
        if type(p) is np.ndarray:
            i=p.astype('int')
        else:
            i=int(p)
        p=p-i
        return i,p
    def nint(self,x):
        """ Return nearest integer of value x """
        return int((x-self.min)/self.step+0.5)
    def print(self):
        print('{}: min={:10.3} step:{:10.3}    o:{}  n:{} max:{:10.3}'.format(
            self.name,self.min,self.step,self.offset,self.n,self.max))
    def write(self,fd):
        fd.write_record(np.int(self.n))
        fd.write_record(np.array([self.min,self.step,self.max],dtype=np.float32))

class table_c(object):
    """ Table class, the name and values of the dependent variable, 
        a linear scale for the independent variable, and an interpolation
        function to compute values from the table.
    """
    kind='cubic'
    fill_value='extrapolate'
    def __init__(self,scale,values,name):
        self.scale=scale
        self.values=values
        self.name=name
        self.interpolate=scipy.interpolate.interp1d(
            scale.values,values,self.kind,fill_value=self.fill_value)
    def write(self,fd):
        fd.write_record(len(self.name))
        fd.write_record(self.name)
        fd.write_record(len(self.values))
        fd.write_record(self.values)

class irregular_table(object):
    """ Template EOS object, with methods to
        1) set up a linear scale of the 1st independent variable
        2) set up 2nd linear independent variable for each 1st variable value
        3) dictionaries with tables for each density
        4) a dictionary with non-mandatory info
    """      
    def set_independent0(self,nd,min,step,name='0'):
        """ Define values of the first independent variable (e.g. log density) """
        self.scale0=linear_scale_c(nd,min,step,name)
        self.scale1=linear_scale_c(0,1e9,0.0)
        self.scale1.max=-1e9
        self.index=self.scale0.index
        self.table={}
        self.scales={}
        self.variables={}

    def set_independent1(self,id,n,min,step,name='x1'):
        """ Define values of the 2nd independent variable """
        self.scales[id]=linear_scale_c(n,min,step,name)
        self.scale1.name=self.scales[id].name
        self.scale1.accumulate(n,min,step)
        self.table[id]={}
        
    def set_dependent(self,id,iv,name,values):
        """ Define values of dependent variables """
        table=table_c(self.scales[id],values,name)
        self.table[id][iv]=table
        self.table[id][name]=table
        self.variables[iv]=name
        self.nv=len(self.variables)
        table.scale.offset=self.scale1.nint(table.scale.min)
    
    def set_info(self,info):
        """ Set non-mandatory information """
        self.info=info
    
    def values(self,ir,iv):
        """ Return x,y for y(x) at log density and variable indices ir,iv """
        table=self.table[ir][iv]
        return table.scale.values,table.values

    def lookup(self,v,x,y,verbose=0):
        """ Lookup a table(x,y) value for variable v
            Vectorized in the computational part (a few mus/pt)
        """
        if verbose:
            self.scale0.print()
            self.scale1.print()
        ix,px=self.scale0.index(x)
        iy,py=self.scale1.index(y)
        if type(ix) is np.ndarray:
            value=np.zeros_like(x)
            t=np.zeros((4,len(x)))
            for j,yj in enumerate(y):
                """ Get the four corner values in t[0:4], using y-indices
                    that are the difference from the global y-index and
                    the indices of the local table minima.
                """
                ixj =ix[j]
                table0=self.table[ixj+0][v]
                table1=self.table[ixj+1][v]
                iyj0=iy[j]-table0.scale.offset
                iyj1=iy[j]-table1.scale.offset
                t[0:2,j]=table0.values[iyj0:iyj0+2]
                t[2:4,j]=table1.values[iyj1:iyj1+2]
            qy=1.0-py
            return (1.-px)*(qy*t[0]+py*t[1])+px*(qy*t[2]+py*t[3])
        else:
            v1=self.table[ix  ][v].interpolate(y)
            v2=self.table[ix+1][v].interpolate(y)
            return (1-px)*v1+px*v2
        
    def write(self,file):
        with scipy.io.FortranFile(file,'w') as fd:
            fd.write_record(self.nv)
            self.scale0.write(fd)
            self.scale1.write(fd)
            for id in range(self.scale0.n):
                self.scales[id].write(fd)
                for iv in range(self.nv):
                    self.table[id][iv].write(fd)
                    
    def all_table_derivatives(self,verbose=0):
        import dispatch.EOS
        dispatch.EOS.all_table_derivatives(self,verbose=verbose)
                    
    def table_derivatives(self,v,verbose=0):
        import dispatch.EOS
        dispatch.EOS.table_derivatives(self,v,verbose=verbose)
