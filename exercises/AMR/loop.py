'''
   Skeleton procedure, which one can use to visualize data in slices, for example
   the "d" structure contains enough information to place the patch data in a larger
   context
'''
from read_patch import *

f = open('data/sedov.log','r')

for i in range(0,5):
    line = f.readline()
    print line
    if line is None:
        print 'at EOF'
        exit
    column = line.split()
    fileroot = column[0]
    d = read_patch (f, fileroot)
