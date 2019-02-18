import fortranfile as ftf
import numpy as np

class OnePebble(object):

    def __init__ (self, id, s, w, m, t, p, v, o, d_to_g=0.01, ppc=0.05):
        self.id = id
        self.s = s                            # radius
        self.w = w                            # statistical weight
        self.m = m                            # mass
        self.t = []                           # list of times (scalar)
        self.p = []                           # list of positions (3-tuplet)
        self.v = []                           # list of velocities (3-tuplet)
        self.d_to_g = d_to_g                  # dust to gas ratio
        self.ppc = ppc                        # initial no. of particles per cell
        self.m_swarm = self.w * d_to_g / ppc  # mass of the swarm particle
        self.owner = []
        
        # add first data point
        self.update(t, p, v, o)

    def update (self, t, p ,v, o):
        self.t.append(t)
        self.p.append(p)
        self.v.append(v)
        self.owner.append(o)

    def pretty_print (self):
        print "ID: ",self.id
        print "Radius:", self.s
        print "Stat. weight:", self.w
        print "Mass: ", self.m
        print "Time(s): ", self.t
        print "Position(s): ", self.p
        print "Velocit(y/ies): ", self.v
        print "Swarm mass: ", self.m_swarm
        print "Owner patch: ", self.o


class PebbleDict(dict):

    def __init__ (self, *args, **kw):
        super(PebbleDict,self).__init__(*args, **kw)
        self.ntot = 0


def read_pebbles_dict(pebbles=None, ip=1, it=0, verbose=False, run='', data='data', origin=(0.0,0.0,0.0), rlimit=None, ppc=None):
    """
    Read pebble data and store them in a dictionary of `OnePebble` objects.
    If `pebbles` is an input keyword and it exists, update it.

    If `rlimit` is not None, only return pebbles that are within `rlimit` radii
    of `origin`.

    """

    name = data + '/' + run + '/{:05d}'.format(it) + '/{:05d}.peb'.format(ip)
    aorigin = np.array(origin)

    # open the file
    with ftf.FortranFile(name,'rb') as f:

        # read the format key
        fmt = f.readInts()

        # read the dimension info
        dim=f.readInts()
        n_repeat=dim[2] # what's this?
        dim=dim[0:2]

        # verbose info
        if verbose:
            print "file: ", name
            print " fmt: ", fmt
            print "dims: ", dim
    
        # read the pebble data
        a = f.readReals()
        a = a.reshape(dim[1],dim[0]).transpose()
    
        # read the index data
        idv = f.readInts()
        idx = idv[0:dim[0]]
        ids = idv[dim[0]:]

        # sort on increasing indices
        s = idx.argsort()
        idx = idx[s]
        a = a[s,:]

        # create dict if not supplied
        if pebbles is None:
            pebbles = PebbleDict()

        # update or intialise dict entries
        n      = dim[0]
        pebbles.ntot += n
        pos    = a[:,0:3]
        vel    = a[:,3:6]
        time   = a[:,6]
        size   = a[:,7]
        weight = a[:,8]
        o = ip
        for i,p,v,t,s,w in zip(idx,pos,vel,time,size,weight):
            if w == 0.0: continue # skip particles with zero weight
            if rlimit is not None:
                if ((p-aorigin)**2).sum() > rlimit: continue # skip particles beyond a certain radius `rlimit`
            if pebbles.has_key(i):
                pebbles[i].update(t,p,v,o)
            else:
                m=s**3*3.0
                if (ppc is not None):
                    pebbles[i] = OnePebble(i, s, w, m, t, p, v, o, ppc=ppc)
                else:
                    pebbles[i] = OnePebble(i, s, w, m, t, p, v, o)

        if verbose:
            print "   index range: [", idx.min(), ",", idx.max(), "]"
            print "position range: [", pos.min(), ",", pos.max(), "]"
            print "velocity range: [", vel.min(), ",", vel.max(), "]"
            print "    time range: [", time.min(), ",", time.max(), "]"
            print "    size range: [", size.min(), ",", size.max(), "]"

        return pebbles


if __name__ == "__main__":
    """
    This script reads pebble files from snapshots and then puts the data into dictionaries
    before writing said dictionaries to file as pickles.
    
    """

    import dispatch_utils as du
    import dispatch_data as dd
    import numpy as np
    import cPickle
    from yt_utils import planet
    
    #run = 'm095r1t10d10p-conv'            # experiment name
    #run = 'm095r1t10d10p2'
    run = 'm095r1t10d10p-conv-prod'
    #run = 'm05r1t00p'
    #run = 'm01r1t00p'
    datadir = '/groups/astro/popovas/lustre/dispatch2/experiments/pebble_atmosphere/data' # location of the data
    planet_mass = 0.95
    radius_multiplier = 1.0
    orbit_radius_au = 1.0  
    pl = planet(mass = planet_mass, rad_mult = radius_multiplier, distance = orbit_radius_au)
    ppc = 0.2
    rlimit = 4.0 * pl.hill_rad

    # change this line!
    snapshots = range(2000,2500,5)

    # Since we're using static refinement, the patch IDs are constant in time and we only need to load
    # them once (to get the IDs).
    pp=du.patches(run=run,data=datadir,iout=snapshots[0])
    pp=[p for p in pp if p.kind[0:9]!='rt_solver']

    import os
    if not os.path.exists(run):
        os.makedirs(run)

    for snap in snapshots:
        pebbles = None
        for p in pp:
            ec = np.sqrt((p.pos**2).sum() + ((p.size*0.5)**2).sum())
            if ec > rlimit: continue # skip patches beyond a certain radius `rlimit`
            pebbles = read_pebbles_dict(pebbles=pebbles,ip=p.id,it=snap,run=run,data=datadir,ppc=ppc)
        # save pebbles
        with open('./{0}/pebbles_{1:05d}.pickle'.format(run,snap),'wb') as handle:
            cPickle.dump(pebbles, handle, protocol=cPickle.HIGHEST_PROTOCOL)
        print snap, len(pebbles)

