import fortranfile as ftf
import numpy as np
import os

def trace_utils():
    '''
    f=trace_files(dir):  Returns a list of *.trc trace data files in dir
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


def trace_files (run='', data='../data', verbose=1):
    '''
    Get a list of .trc trace files
    '''
    rundir=data+'/'+run
    filenames = os.listdir(rundir)
    files=[rundir+'/'+f for f in filenames if f.endswith(".trc")]
    if verbose==1: print np.size(files), 'files found'
    return np.sort(files)

def trace(file,verbose=1):
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
    class p:
        pass
    with ftf.FortranFile(file,'rb') as f:
        j=f.readInts()
        p.status=j[1]
        size=f.readReals()
        p.size=size[0]
        id,rec,n=f.readInts()
        list=[]
        # flatten records, which contain several trace points, into list
        while id:
            j=f.readInts()
            a=f.readReals()
            a=a.reshape(n,np.size(a)/n)
            for i in range(n):
                q=p()
                q.ip=j[i]
                q.a=a[i,:]
                list.append(q)
            ii=f.readInts()
            id,rec,n=ii[0:3]
    # extract structure fields from list
    p.n=np.size(list)
    p.ip=np.array([i.ip     for i in list])
    p.p =np.array([i.a[0:3] for i in list]).transpose()
    p.v =np.array([i.a[3:6] for i in list]).transpose()
    p.vg=np.array([i.a[6:9] for i in list]).transpose()
    p.c =np.array([i.a[9]   for i in list])
    p.d =np.array([i.a[10]  for i in list])
    p.t =np.array([i.a[11]  for i in list])
    p.x = p.p[0]
    p.y = p.p[1]
    p.z = p.p[2]
    r=np.sqrt(np.sum(p.p**2,0))
    p.rmin=r.min()
    p.rmax=r.max()
    p.first=p.p[:,0]
    p.last=p.p[:,-1]
    if verbose==1: print file, p.n, 'points'
    return p

def traces(run='', dir='../data', verbose=1):
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
    files=trace_files(dir=dir,run=run,verbose=verbose)
    s=[]
    class tr:
        file=''
        def trace(self):
            return trace(self.file,verbose=0)
    for f in files:
        p=trace(f,verbose=0)
        a=tr()
        a.rmin=p.rmin
        a.rmax=p.rmax
        a.first=p.first
        a.last=p.last
        a.ip=p.ip[0]
        a.file=f
        s.append(a)
        if verbose>1: print '{}: {:10.3f} {:10.3f}'.format(a.file,a.rmin,a.rmax)
    return s

def cache_path(tr):
    '''Dress a trace object with actual path'''
    tr.path=tr.trace()
    return tr

def cache_paths(trs):
    '''Dress a list of objects with path data'''
    for tr in trs:
        cache_path(tr)