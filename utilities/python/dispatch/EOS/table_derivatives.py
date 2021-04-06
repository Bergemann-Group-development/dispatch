import numpy as np
from dispatch.EOS import deriv

def table_derivatives(eos,v,verbose=0):
    if verbose>0:
        print(v,eos.variables[v])
    """ Add derivates for one table variable; the derivatives
        in the two axis directions (expressed in grid units)
        are added as attributes deriv0 and deriv1 of the tables.
    """

    def left1(eos,i,o,o10,o00,verbose=0):
        """ 1st order, left axis=0 derivative """
        if verbose>1:
            print('left1',i,o)
        d=eos.tables[i][v].deriv0
        d[o-o10]=(eos.tables[i  ][v].values[o-o10]-
                  eos.tables[i-1][v].values[o-o00])

    def right1(eos,i,o,o10,o20,verbose=0):
        """ 1st order, right axis=0 derivative """
        if verbose>1:
            print('right1',i,o)
        d=eos.tables[i][v].deriv0
        d[o-o10]=(eos.tables[i+1][v].values[o-o20]-
                  eos.tables[i  ][v].values[o-o10])

    def left(eos,i):
        tables=eos.tables[i][v]
        tables.deriv0=np.zeros(tables.scale.n)
        tables.deriv1=deriv(tables.scale.values,tables.values)
        o00=eos.tables[i-1][v].scale.offset
        o10=eos.tables[i  ][v].scale.offset
        o01=eos.tables[i-1][v].scale.n+o00
        o11=eos.tables[i  ][v].scale.n+o10
        o0=max(o00,o10)
        o1=min(o01,o11)
        for o in range(o0,o1):
            left1(eos,i,o,o10,o00)

    def right(eos,i):
        tables=eos.tables[i][v]
        tables.deriv0=np.zeros(tables.scale.n)
        tables.deriv1=deriv(tables.scale.values,tables.values)
        o10=eos.tables[i  ][v].scale.offset
        o20=eos.tables[i+1][v].scale.offset
        o11=eos.tables[i  ][v].scale.n+o10
        o21=eos.tables[i+1][v].scale.n+o20
        o0=max(o20,o10)
        o1=min(o21,o11)
        for o in range(o0,o1):
            right1(eos,i,o,o10,o20)

    def central(eos,i,verbose=0):
        tables=eos.tables[i][v]
        tables.deriv0=np.zeros(tables.scale.n)
        tables.deriv1=deriv(tables.values)
        o00=eos.tables[i-1][v].scale.offset
        o10=eos.tables[i  ][v].scale.offset
        o20=eos.tables[i+1][v].scale.offset
        o01=eos.tables[i-1][v].scale.n+o00
        o11=eos.tables[i  ][v].scale.n+o10
        o21=eos.tables[i+1][v].scale.n+o20
        o0=max(o00,o10,o20)
        o1=min(o01,o11,o21)
        f=np.zeros(3)
        d=tables.deriv0
        for o in range(o0,o1):
            f[0]=eos.tables[i-1][v].values[o-o00]
            f[1]=eos.tables[i  ][v].values[o-o10]
            f[2]=eos.tables[i+1][v].values[o-o20]
            d[o-o10]=deriv(f)[1]
        for o in range(o10,o0):
            if o-o00 > 0:
                left1(eos,i,o,o10,o00,verbose=verbose)
            if o-o20 > 0:
                right1(eos,i,o,o10,o00,verbose=verbose)
        for o in range(o1,o11):
            if o-o00 < eos.tables[i-1][v].scale.n:
                left1(eos,i,o,o10,o00,verbose=verbose)
            if o-o20 < eos.tables[i+1][v].scale.n:
                right1(eos,i,o,o10,o20,verbose=verbose)

    right(eos,0)
    left(eos,eos.scale0.n-1)
    for i in range(1,eos.scale0.n-1):
        central(eos,i,verbose=verbose)

def all_table_derivatives(eos,verbose=0):
    """ Add EOS derivatives (in index units) in the axis=0 and axis=1 directions
    """
    for v in range(0,eos.nv):
        table_derivatives(eos,v,verbose)