import os
import numpy as np

class JCD_S(object):
    def __init__(self,file='JCD_model_S.txt',top='~/codes/dispatch2',verbose=0):
        os.chdir(top+'/data/solar_models')
        with open(file,'r') as fd:
            a=[]
            for line in fd:
                line=line.rstrip()
                if line.startswith('#'):
                    if verbose>0:
                        print(line)
                else:
                    a.append([float(word) for word in line.split()])
        a=np.array(a).T
        self.r     =a[0]
        self.cs    =a[1]
        self.rho   =a[2]
        self.P     =a[3]
        self.gamma1=a[4]
        self.T     =a[5]
        if verbose>0:
            print(vars(self).keys())
    