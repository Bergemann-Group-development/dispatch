# -*- coding: utf-8 -*-
"""
    Save and restore an object or dict

    save(run+'.save',dict)
    dict=restore(run+'.save')

    If the file does not exist, an empty dict is returned

Created on Wed Jan 9 03:25:50 2017

@author: Åke Nordlund
"""

import pickle

def save(file,var):
    fd=open(file,'wb')
    pickle.dump(var,fd)
    fd.close()

def restore(file):
    try:
        fd=open(file,'rb')
        var=pickle.load(fd)
        fd.close()
        return var
    except:
        return {}
