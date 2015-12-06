'''
   Read info in the data/*.txt file
'''
import numpy as np

f = open('data/sedov_001_0001.txt','r')
# line 1: time pos(1:3) size(1:3)
line = f.readline()
column = line.split()
time = float(column[0])
pos = [float(column[1]),float(column[2]),float(column[3])]
size = [float(column[4]),float(column[5]),float(column[6])]
print time, pos, size
# line 2: m(1:3) n(1:3) li(1:3) ui(1:3)
line = f.readline()
column = line.split()
m = [int(column[0]),int(column[1]),int(column[2])]
n = [int(column[3]),int(column[4]),int(column[5])]
li = [int(column[6]),int(column[7]),int(column[8])]
ui = [int(column[9]),int(column[10]),int(column[11])]
print m, n, li, ui