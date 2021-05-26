import mmap
import numpy as np

class memmap2(np.memmap):
    """ Attempt at minimal np.memmap subclass """
    def __new__(subtype, *args, **kwargs):
        obj = super(memmap2, subtype).__new__(subtype,  *args, **kwargs)
        obj._saved = {'len':obj._mmap.__len__(),'off':obj.offset,
                      'file':obj.filename}
        #fo=open(obj.filename,'rb+')
        obj._mmap.close()
        return obj
    def reopen(self):
        s=self._saved
        if self._mmap.closed:
            print('was closed, reopen')
            fo=open(s['file'],'rb+')
            self._mmap=mmap.mmap(fo.fileno(),offset=s['off'],length=s['len'])
        else:
            print('was open')
    def __array_prepare__ (self, *args, **kwargs):
        s=self._saved
        if self._mmap.closed:
            print('was closed, reopen')
            fo=open(s['file'],'rb+')
            self._fileobj=fo
            self._mmap=mmap.mmap(fo.fileno(),offset=s['off'],length=s['len'])
        else:
            print('was open')
        obj = super(memmap2, self).__array_prepare__(self, *args, **kwargs)
        return obj

class memmap3(np.memmap):
    """ Working np.memmap subclass, with option to close after access """
    def __new__(subtype, file, dtype=np.float32, mode='r+', offset=0, 
                shape=(10,), order='F'):
        self = super(memmap3, subtype).__new__(subtype, file, dtype=dtype, 
                    mode=mode, offset=offset, shape=shape, order=order)
        #print('-> memmap3.__new__')
        length = self._mmap.__len__()
        offset = self.offset
        saved = (file,length,offset)
        print('-> memmap3.__new__: saved=',saved)
        self._saved = saved
        self._mmap.close()
        return self
    def reopen(self):
        print("-> memmap3.reopen")
        print(self._saved)
        (file,length,offset)=self._saved
        if self._mmap.closed:
            #print('--> was closed, reopen')
            fo=open(file,'rb+')
            block=mmap.ALLOCATIONGRANULARITY
            o=(offset//block)*block
            self._mmap=mmap.mmap(fo.fileno(),length,offset=o)
    def __array_wrap__(self, arr, context=None):
        #print('-> memmap3.__array_wrap__')
        arr = super(memmap3, self).__array_wrap__(arr)
        return arr
    def __array_prepare__ (self, obj, context=None):
        #print('-> memmap3.__array_prepare__')
        self.reopen()
        obj = super(memmap3, self).__array_prepare__(obj)
        return obj
    def __array_finalize__ (self, obj):
        #print('-> memmap3.__array_finalize__')
        super(memmap3,self).__array_finalize__(obj)
        if obj is None:
            return
        self._saved = obj._saved
        self._mmap = obj._mmap
    def __repr__(self):
        #print('-> memmap3.__repr__')
        self.reopen()
        output = super(memmap3,self).__repr__()
        return output
    def __str__(self):
        #print('-> memmap3.__str__')
        self.reopen()
        output = super(memmap3,self).__str__()
        return output

class memmap4(np.memmap):
    """ Working np.memmap subclass, with verbose option """
    def __new__(subtype, file, dtype=np.float32,mode='w+',offset=0, 
                shape=(10,), order='F',verbose=False):
        self = super(memmap4, subtype).__new__(subtype,file,dtype=dtype, 
                    mode=mode,offset=offset,shape=shape,order=order)
        self._verbose=verbose
        self._saved=(file,self._mmap.__len__(),self.offset)
        self._mmap.close()
        if self._verbose:
            print('-> memmap4.__new__: saved =',self._saved,
                  type(self.offset),type(offset))
        return self
    def verbose(self):
        return hasattr(self,'_verbose') and self._verbose
    def reopen(self):
        (file,length,offset)=self._saved
        if self.verbose():
            print("-> memmap4.reopen: length, offset =",length,offset)
        if self._mmap.closed:
            if self.verbose():
                print('--> was closed, reopen')
            fo=open(file,'rb+')
            self._fileobj=fo
            block=mmap.ALLOCATIONGRANULARITY
            o=(offset//block)*block
            self._mmap=mmap.mmap(fo.fileno(),length,offset=o)
        else:
            if self.verbose():
                print('--> was already open')
    def __array_wrap__(self, arr, context=None):
        if self.verbose():
            print('-> memmap4.__array_wrap__')
        arr = super(memmap4, self).__array_wrap__(arr, context)
        return arr
    def __array_prepare__ (self, obj, context=None):
        if self.verbose():
            print('-> memmap4.__array_prepare__')
        self.reopen()
        obj = super(memmap4, self).__array_prepare__(obj, context)
        return obj
    def __array_finalize__ (self, obj):
        super(memmap4,self).__array_finalize__(obj)
        if self.verbose():
            print('-> memmap4.__array_finalize__')
        if obj is None:
            return
        self._saved = obj._saved
        self._mmap = obj._mmap
    def __repr__(self):
        self.reopen()
        if self.verbose():
            print('-> memmap4.__repr__')
        output = super(memmap4,self).__repr__()
        return output
    def __str__(self):
        if self.verbose():
            print('-> memmap4.__str__')
        self.reopen()
        output = super(memmap4,self).__str__()
        return output
