from read_patch import *

f = open('data/sedov.log','r')

for i in range(0,5):
    line = f.readline()
    if line is None:
        print 'at EOF'
        exit
    column = line.split()
    fileroot = column[0]
    d = read_patch (fileroot)
