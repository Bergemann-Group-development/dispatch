import numpy as np
import scipy
import scipy.io
import scipy.interpolate
from copy import deepcopy

def write_str(fd,s):
    """ Write a string to an unformatted Fortran file """
    fd.write_record(np.int(len(s)))
    fd.write_record(bytes(s,'utf-8'))

def read_str(fd):
    """ Read a string from an unformatted Fortran file """
    return(''.join(map(chr,fd.read_record('b'))))

class independent_c(object):
    """ Linear scale for independent variables """
    def __init__(self,n=0,min=0.0,step=0.0,name=None,fd=None,verbose=0):
        if not fd is None:
            self.read(fd,verbose=verbose)
        else:
            self.setup(n,min,step,name,verbose=verbose)

    def setup(self,n,min,step,name='',verbose=0):
        self.name=name
        self.n=n
        self.min=min
        self.max=min+(n-1)*step
        self.step=step
        self.values=min+step*np.arange(n)
        self.offset=0
        if verbose>2:
            self.print()

    def clone(self,min,max):
        """ Create a new, shorter scale, with integer offset """
        new=deepcopy(self)
        i0=int((min-self.min)/self.step+1)
        i1=int((max-self.min)/self.step)
        new.offset=i0
        new.min=self.min+i0*new.step
        new.max=self.min+i1*new.step
        new.n=i1-i0+1
        new.values=new.min+new.step*np.arange(new.n)
        return new

    def accumulate(self,scale,verbose=0):
        """ Accumulate index ranges """
        step=scale.step
        if verbose>2:
            self.print()
            scale.print()
        mn=min(self.min,scale.min)
        mx=max(self.max,scale.max)
        n=int((mx-mn)/step+1.5)
        self.setup(n,mn,step,scale.name,verbose=verbose)

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

    def print(self,label=None):
        s=self.name if label is None else label
        print('{}   min={:10.3} step:{:10.3}    o:{}  n:{} max:{:10.3}'.format(
            s,self.min,self.step,self.offset,self.n,self.max))

    def write(self,fd):
        write_str(fd,self.name)
        fd.write_record(np.int(self.n))
        fd.write_record(np.array([self.min,self.step,self.max],dtype=np.float64))

    def read(self,fd,verbose=0):
        name=read_str(fd)
        n=fd.read_record('i4')[0]
        min,step,max=fd.read_record('f8,f8,f8')[0]
        if verbose>2:
            print(name,min,step,max)
        self.setup(name,n,min,step)

class dependent_c(object):
    """ Table class, the name and values of the dependent variable,
        a linear scale for the independent variable, and an interpolation
        function to compute values from the table.
    """
    kind='linear'
    fill_value='extrapolate'
    def __init__(self,scale=None,values=None,name=None,fd=None,verbose=0):
        if not scale is None:
            self.setup(scale,values,name)
        if not fd is None:
            self.read(fd,verbose=verbose)

    def setup(self,scale,values,name):
        self.scale=scale
        self.values=values
        self.name=name
        self.interpolate=scipy.interpolate.interp1d(
            scale.values,values,self.kind,fill_value=self.fill_value)

    def write(self,fd):
        write_str(fd,self.name)
        fd.write_record(len(self.values))
        fd.write_record(self.values.astype(np.float64))

    def read(self,fd,verbose=0):
        self.name=read_str(fd)
        self.n=fd.read_record('i4')[0]
        self.values=fd.read_record('f8')

class itable_c(object):
    """ Template EOS object, with methods to
        1) set up a linear scale of the 1st independent variable
        2) set up 2nd linear independent variable for each 1st variable value
        3) dictionaries with tables for each density
        4) a dictionary with non-mandatory info
    """
    def __init__(self,file=None,verbose=0):
        if not file is None:
            self.read(file,verbose)

    def setup(self,nv):
        """ Set the initial number of variables """
        self.nv=nv
        self.table={}
        self.scales={}
        self.variables={}

    def set_independent0(self,nd,min,step,name='0',verbose=0):
        """ Define values of the first independent variable (e.g. log density) """
        self.scale0=independent_c(nd,min,step,name,verbose=verbose)
        self.scale1=independent_c(0,1e9,0.0,verbose=verbose)
        self.scale1.max=-1e9
        self.index=self.scale0.index
        self.tables={}
        self.scales={}
        self.variables={}

    def set_independent1(self,id,n,min,step,name='x1',verbose=0):
        """ Define values of the 2nd independent variable """
        self.scales[id]=independent_c(n,min,step,name,verbose=verbose)
        self.scale1.accumulate(self.scales[id],verbose=verbose)
        self.tables[id]={}

    def set_dependent(self,id,iv,values,name=None,verbose=0):
        """ Define values of dependent variables """
        table=dependent_c(self.scales[id],values,name,verbose=verbose)
        if verbose==3:
            print('dependent:',id,iv,name)
        self.tables[id][iv]=table
        self.tables[id][name]=table
        if type(name) is str:
            self.variables[iv]=name
            self.nv=len(self.variables)
        table.scale.offset=self.scale1.nint(table.scale.min)
        if verbose>3:
            print(id,iv,self.scales[id].n,table.values.shape[0])

    def first_independent1(self,n,min,max,name='y',verbose=0):
        """ Define values of the 2nd independent variable """
        step=(max-min)/(n-1)
        min=int(min/step+1)*step
        self.scale1=independent_c(n,min,step,name,verbose=verbose)
        self.scales={}

    def next_independent1(self,id,min,max,new_tables=True,verbose=0):
        """ Clone a local scale off the main one """
        self.scales[id]=self.scale1.clone(min,max)
        if new_tables:
            self.tables[id]={}

    def set_info(self,info):
        """ Set non-mandatory information """
        self.info=info

    def values(self,ir,iv,verbose=0):
        """ Return x,y for y(x) at log density and variable indices ir,iv """
        if verbose>2:
            print(ir,iv,len(self.tables.keys()),self.tables[ir].keys())
        table=self.tables[ir][iv]
        return table.scale.values,table.values

    def lookup(self,v,x,y,verbose=0):
        """ Lookup a table(x,y) value for variable v
            Vectorized in the computational part (a few mus/pt)
        """
        if verbose>1:
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
                table0=self.tables[ixj+0][v]
                table1=self.tables[ixj+1][v]
                iyj0=iy[j]-table0.scale.offset
                iyj1=iy[j]-table1.scale.offset
                t[0:2,j]=table0.values[iyj0:iyj0+2]
                t[2:4,j]=table1.values[iyj1:iyj1+2]
            qy=1.0-py
            return (1.-px)*(qy*t[0]+py*t[1])+px*(qy*t[2]+py*t[3])
        else:
            v1=self.tables[ix  ][v].interpolate(y)
            v2=self.tables[ix+1][v].interpolate(y)
            return (1-px)*v1+px*v2

    def write(self,file,verbose=0):
        """ Write out an EOS table as unformatted Fortran records """
        with scipy.io.FortranFile(file,'w') as fd:
            fd.write_record(self.nv)
            if verbose>2:
                print('nv:',self.nv)
                print('scale0')
            self.scale0.write(fd)
            if verbose>2:
                print('scale1')
            self.scale1.write(fd)
            for id in range(self.scale0.n):
                self.scales[id].write(fd)
                nv=0
                for iv in self.tables[id].keys():
                    if type(iv) is int:
                        nv+=1
                if verbose>1:
                    print(id,nv)
                for iv,table in self.tables[id].items():
                    if type(iv) is int:
                        self.tables[id][iv].write(fd)
                        if verbose>2:
                            print('write:',id,iv,self.tables[id][iv].values.shape[0])

    def read(self,file,verbose=0):
        """ Read an EOS file produce by write() """
        if verbose>1:
            print('reading',file)
        self.file=file
        with scipy.io.FortranFile(file,'r') as fd:
            nv=fd.read_record('i4')[0]
            self.setup(nv)
            self.scale0=independent_c()
            self.scale1=independent_c()
            self.scale0.read(fd=fd,verbose=verbose)
            self.scale1.read(fd=fd,verbose=verbose)
            self.scale1.max=-1e9
            self.index=self.scale0.index
            if verbose>2:
                print('nv:',self.nv,'n:',self.scale0.n)
            for id in range(self.scale0.n):
                """ Setup local scale for dependent variables """
                s1=independent_c(fd=fd,verbose=verbose)
                s1.offset=self.scale1.nint(s1.min)
                if verbose>3:
                    print(id,s1.file,s1.n,s1.values[0])
                self.scales[id]=s1
                self.scale1.accumulate(s1)
                self.tables[id]={}
                """ Read dependent variables """
                nv=fd.read_ints()[0]
                for iv in range(nv):
                    table=dependent_c(fd=fd)
                    table.scale=s1
                    if verbose>3:
                        print('iv:',iv,table.name)
                    self.tables[id][iv]=table
                    self.tables[id][table.name]=table
                    self.variables[iv]=table.name
                    if verbose>3:
                        print('read:',id,iv,nv,table.values.shape[0])
        self.nv=len(self.variables)

    def limits(self,iv):
        """ Return the independent variable and limits of the variable iv """
        indep=[]
        min=[]
        max=[]
        for i in range(self.scale0.n):
            if iv in self.tables[i]:
                v=self.tables[i][iv].values
                min.append(v.min())
                max.append(v.max())
                indep.append(self.scale0.values[i])
        return np.array(indep),np.array(min),np.array(max)

    def square(self,iv,n=0,verbose=0):
        """ Return the full square extent and values of a variable """
        if n==0:
            n=self.scale0.n
        n1=self.scale1.n+1
        f=np.zeros((n,n1))
        if verbose>2:
            print('shape of f:',f.shape)
        n1max=0
        n0max=0
        for i in range(n):
            if iv in self.tables[i]:
                o=self.scales[i].offset
                n=self.scales[i].n
                if verbose>3:
                    sz=self.tables[i][iv].values.shape[0]
                    print(i,n,sz,o+n-1,f.shape[1])
                f[i,o:o+n]=self.tables[i][iv].values
                n0max=max(n0max,i)
                n1max=max(n1max,o+n)
        s1=self.scale1
        x=self.scale0.values[0:n0max]
        y=self.scale1.values[0:n1max]
        f=f[:n0max,:n1max]
        e=[x.min(),x.max(),y.min(),y.max()]
        return e,f

    def extent(self,iv,n=0,verbose=0):
        """ Return the extent of the table """
        e,f=self.square(iv,n,verbose)
        return e

    def imshow(self,iv,n=0,extent=None,cmap=None,verbose=0):
        """ Show an colormap of a variable """
        if iv in self.variables.keys() or iv in self.variables.values():
            import matplotlib.pyplot as pl
            e,f=self.square(iv,n,verbose=verbose)
            if extent is None:
                extent=e
            if type(iv) is not str:
                iv=self.variables[iv]
            if iv.find('ln')>=0:
                f=f/np.log(10)
                iv=iv.replace('ln','log10')
            name0=self.scale0.name
            if name0.find('ln')>=0:
                extent[0]=extent[0]/np.log(10)
                extent[1]=extent[1]/np.log(10)
                name0=name0.replace('ln','log10')
            pl.xlabel(name0)
            name1=self.scale1.name
            if name1.find('ln')>=0:
                extent[2]=extent[2]/np.log(10)
                extent[3]=extent[3]/np.log(10)
                name1=name1.replace('ln','log10')
            pl.ylabel(name1)
            pl.imshow(f.T,origin='lower',cmap=cmap,extent=extent,aspect='auto')
            pl.title(iv)
            pl.colorbar()
            return np.array(extent)
        else:
            print(iv,'is not among variables',self.variables)

    def all_table_derivatives(self,verbose=0):
        """ Add all table derivatives """
        import dispatch.EOS
        dispatch.EOS.all_table_derivatives(self,verbose=verbose)

    def table_derivatives(self,v,verbose=0):
        """ Add derivates for one table variable; the derivatives
            in the two axis directions (expressed in grid units)
            are added as attributes deriv0 and deriv1 of the table.
        """
        import dispatch.EOS
        dispatch.EOS.table_derivatives(self,v,verbose=verbose)

    def show(self,verbose=0,plot=False):
        """ Show text and optionally graphic information about an EOS """
        if verbose>1:
            print('type:',self.type)
            print('file:',self.file)
            self.scale0.print(label='axis 0: '+self.scale0.name)
            self.scale1.print(label='axis 1: '+self.scale1.name)
        if verbose>0:
            print('variables:',self.variables)
        s0=self.scale0
        s1=self.scale1
        if verbose>2:
            print('independent variables, axis 0:')
            s0.print()
            print('independent variables, axis 1:')
            s1.print()
        if plot:
            import matplotlib.pyplot as pl
            if s0.name.find('ln')>=0:
                c0=1./np.log(19)
                name0=s0.name.replace('ln','log10')
            else:
                c0=1.
                name0=s0.name
            if s1.name.find('ln')>=0:
                c1=1./np.log(19)
                name1=s1.name.replace('ln','log10')
            else:
                c1=1.
                name1=s1.name
            s0=self.scale0
            for i in range(s0.n):
                si=self.scales[i]
                pl.plot(s0.values[i]*c0,si.min*c0,'+',color='g')
                pl.plot(s0.values[i]*c1,si.max*c1,'+',color='r')
            pl.xlabel(name0)
            pl.ylabel(name1)
            pl.title('type:'+self.type+'     file:'+self.file)
            pl.tight_layout()
