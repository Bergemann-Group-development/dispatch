import numpy as np

#=====================================================================
class MemmapList:
    #-----------------------------------------------------------------
    def __init__(self, memmaps, axis=-1):
        self.mlist = [m for m in memmaps if m.size > 0]
        self.shape = np.array(self.mlist[0].shape)
        self.shape[axis] = np.sum([m.shape[axis] for m in self.mlist])
        self.size = np.prod(self.shape)
        self.ndim = self.shape.size
        self.caxis = axis
        if axis == -1: self.caxis = self.shape.size - 1
        self.cstart = np.cumsum([0,] + [m.shape[axis] for m in self.mlist[:-1]])

    def __getitem__(self, index):
    #-----------------------------------------------------------------
        if self.ndim == 1:
            return np.array(self)[index]
        
        if np.issubdtype(type(index), np.integer):
            index = (index,)

        if type(index) == slice:
            index = (index,)

        if Ellipsis in index:
            Edim = self.ndim - len(index) + 1
            iE = index.index(Ellipsis)
            index = index[:iE] + (slice(None),)*Edim + index[iE+1:]

        dim = len([i for i in index if i != None])
        if dim > self.ndim:
            raise IndexError("too many indices")
        elif dim-1 >= self.caxis:
            cutc = index[self.caxis] != slice(None)
        else:
            cutc = False

        idim = [i for i, v in enumerate(index) if isinstance(v, int)]
        ncut = len([i for i in idim if i < self.caxis])
            
        if cutc:
            slc = index[self.caxis]
            idx = np.array(index)
            if np.issubdtype(type(slc), np.integer):
                if slc < 0: slc += self.shape[self.caxis]
                i = np.searchsorted(self.cstart, slc, side='right') - 1
                m = self.mlist[i]
                idx[self.caxis] = slc - self.cstart[i]
                return self.mlist[i][tuple(idx)]
            else:
                mlist = []
                for i, m in enumerate(self.mlist):
                    #start = slc.start
                    #stop = slc.stop
                    #if start is not None: start = max(0, start - self.cstart[i])
                    #if stop is not None: stop = max(0, stop - self.cstart[i])

                    i1 = self.cstart[i]
                    o1 = i1 + m.shape[self.caxis]
                    i2 = slc.start
                    o2 = slc.stop
                    s2 = slc.step
                    n = m.shape[self.caxis]

                    if s2 is None: s2 = 1
                    if s2 > 0:
                        start = np.ceil(i1/s2).astype(int)*s2 - i1
                        stop  = np.ceil(o1/s2).astype(int)*s2 - i1
                        if i2 is not None: start = max(start, i2 - i1)
                        if o2 is not None: stop  = min(stop,  o2 - i1)
                        if stop < 0: stop = 0
                    else:
                        start = n - np.floor((n-o1)/s2).astype(int)*s2 - i1 - 1
                        stop  = n - np.floor((n-i1)/s2).astype(int)*s2 - i1 - 1
                        if i2 is not None: start = min(start, i2 - i1)
                        if o2 is not None: stop  = max(stop,  o2 - i1)
                        if stop < 0: stop = None

                    idx[self.caxis] = slice(start, stop, slc.step)
                    mlist.append(m[tuple(idx)])

                if slc.step is not None:
                    if slc.step < 0: mlist = mlist[::-1]

        else:
            mlist = [m[index] for m in self.mlist]

        if len(mlist)>1:
            return MemmapList(mlist, axis=self.caxis-ncut)
        else:
            return mlist[0]
        
    def __array__(self, dtype=None):
        return np.concatenate(self.mlist, axis=self.caxis)
    
    def squeeze(self, axis=None):
        slc = np.array([slice(None)]*self.ndim)
        slc[self.shape==1] = 0
        slc[self.caxis] = slice(None)
        return self[tuple(slc)]
    
    def flatten(self):
        return np.array(self).flatten()
    
    def ravel(self):
        return self.flatten()

    def __add__(self, other):
        return np.array(self) + other

    def __sub__(self, other):
        return np.array(self) - other

    def __mul__(self, other):
        return np.array(self) * other

    def __truediv__(self, other):
        return np.array(self) / other
    
    def __transpose__(self):
        return np.array(self).transpose()