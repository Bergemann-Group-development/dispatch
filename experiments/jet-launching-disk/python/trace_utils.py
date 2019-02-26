import fortranfile as ftf
import numpy as np
import os
import time

def timeit(f):
    def timed(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print '{} {} took: {:2.4f} sec'.format(f.__name__, args, te-ts)
        return result
    return timed

def trace_utils():
    '''
    f=trace_files(run):  Returns a list of *.trc trace data files in data/run
    tr=trace(file):      Returns a structure with info for one particle
    trs=traces(dir):     Returns a list of structures, one for each particle
    cache_traces(trs):   Dress a list of trace objects with actual paths
    
    # Example traces() use, plotting trace x- and y-components:
    for tr in traces(dir):
        if tr.rmin < 100:
            s=tr.trace()
            plot(s.x,s.y)
        ...            
    # The construct is analogous to one that loops over files:
    for file in trace_files(dir):
        tr=trace(file)
        if tr.rmin < 100:
            plot(tr.x,tr.y)
        ...      
    # To select a subset based on e.g. rmin, use syntax such as:  
    subset=[tr for tr in trace(run) if tr.rmin < 200]
    '''
    print "Type trace_utils?<RETURN> for help text, or trace_utils(\
 for pop-up help in Canopy"

class trace:
    def __init__(self, file, cache=False, verbose=0):
        '''
        Return a structure with trace data from a trace file, with fields
        {n,p,v,vg,c,d,t}, where p[0:3,:], v[0:3,:], vg[0:3,:] are the positions,
        velocities, and gas velocities, and c[:], d[:], t[:] are the drag 
        coefficients, densities, and times.  Fields x, y, and z are
        shorthand for the arrays p[0], p[1], and p[2], respectively.

        Example use:

            s=traces()
            path=s[0].trace()
            plot(path.t,path.x)
            mayavi.mlab.plot3d(path.x,path.y,path.z)
        '''
        self.read (file, cache=cache, verbose=verbose)

    def read (self, file, cache=False, verbose=0):
        self.file=file;
        with ftf.FortranFile(file,'rb') as f:
            p=[]
            v=[]
            vg=[]
            c=[]
            d=[]
            t=[]
            nt=0
            try:
                j=f.readInts()
                self.status=j[1]
                sz=f.readReals()
                self.size=sz[0]
                id,rec,n=f.readInts()
                self.id=id
                if verbose>2:
                    print id
                elif verbose>1:
                    nt+=1
                    if mod(nt,100)==0:
                        print id
                # flatten records, which contain several trace points, into list
                t_prev=0.0
                while id:
                    j=f.readInts()
                    a=f.readReals()
                    if n>0:
                        a=a.reshape(n,np.size(a)/n)
                        for i in range(n):
                            p.append(a[i,0:3])
                            t_now=a[i,11]
                            if t_now < t_prev:
                                print '{:s} time reversal: {:.4f} > {:.4f}'.format(file,t_prev,t_now)
                            t_prev=t_now
                            t.append(t_now)
                            if cache:
                                v.append(a[i,3:6])
                                vg.append(a[i,6:9])
                                c.append(a[i,9])
                                d.append(a[i,10])
                    else:
                        print 'skipping zero length trace',file
                    ii=f.readInts()
                    id,rec,n=ii[0:3]
            except:
                print file,"read error after",np.size(list)
        # extract structure fields from list
        if list==[]:
            return None
        else:
            self.n =np.size(t)
            p=np.array(p)
            if cache:
                self.p=p.transpose()
                self.t=np.array(t)
                self.v=np.array(v).transpose()
                self.vg=np.array(vg).transpose()
            r=np.sqrt(np.sum(p**2,1))
            self.rmin=r.min()
            self.rmax=r.max()
            self.first=p[0,:]
            self.last=p[-1,:]
            if verbose>0: 
                print file, self.n, 'points'

    def cache (self, verbose=0):
        '''Dress a trace object with actual path'''
        self.path=self.read (self.file, cache=True, verbose=verbose)

def aspect(select):
    Wx=0.0; Wy=0.0; Wz=0.0
    for tr in select:
        Wx=max(Wx,abs(tr.first[0]))
        Wy=max(Wy,abs(tr.first[1]))
        Wz=max(Wz,abs(tr.first[2]))
    return Wx,Wy,Wz

class traces:
    def __init__(self, run='', data='../data', cache=True, verbose=0):
        '''
        Return a list of structures {file,trace,rmin,rmax,first,last},
        where self.trace() is a function returning a trace(self.file)
        structure, rmin,rmax are the min and max radii for the trace, 
        and first,last are the first and last 3D trace points

        To select a subset based on e.g. rmin and cache the paths,
        use syntax such as:  

            subset=[tr for tr in trace(run) if tr.rmin < 200]
            for tr in subset: tr.path=tr.trace()

        See trace( for more about the contents of the path structure
        '''
        rundir=data+'/'+run
        filenames = os.listdir(rundir)
        files=[rundir+'/'+f for f in filenames if f.endswith(".trc")]
        self.files=np.sort(files)
        self.run=run
        if (verbose>=0):
            print np.size(files),'trace file names in self.files'
        if cache:
            self.info(cache=True)
        elif verbose>=0:
            print 'use self.info() to read basic information into self.traces'
            print 'use self.cache() to cache full path information into all self.traces[].path,'
            print 'or use self.traces[i].cache() to cache individual path information'

    def info (self, cache=False, verbose=0):
        self.traces=[]
        self.tmin=1e10
        self.tmax=0.0
        for file in self.files:
            tr=trace(file, cache=cache, verbose=verbose)
            if any(np.isnan(tr.first)):
                print 'skipping trace with nan', tr.id
            else:
                self.traces.append (tr)
                self.tmin=min(self.tmin,tr.t[0])
                self.tmax=max(self.tmax,tr.t[-1])
            if verbose>1:
                print '{}: {:10.3f} {:10.3f} (s)'.format(file,tr.rmin,tr.rmax)

    #@timeit
    def cache (self, verbose=0):
        '''Dress a list of objects with path data'''
        if not hasattr(self,'traces'):
            if verbose:
                print 'getting trace data'
            self.info (verbose=verbose)
        for tr in self.traces:                                
            tr.cache(verbose=verbose)
