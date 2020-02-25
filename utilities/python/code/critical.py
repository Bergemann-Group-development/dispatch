# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:14:06 2019

@author: Aake
"""

import os
import sys

#%%
def critical(file,verbose=0):
    """ Find and check legality of double if """
    with open(file) as fo:
        n0=0        # line number
        watch=False
        ok=True
        if verbose>1:
            print('')
            print(file)
        for line in fo:
            n0+=1
            line=line.lstrip()

            # detect critical-block start
            if line.startswith('!$omp critical'):
                if (line.find('test_cr') < 0 and 
                    line.find('list_cr') < 0 and
                    line.find('cleanup_cr')) < 0:
                    watch=True

            # detect 'first_time ='
            elif watch:
                if line.startswith('call'):
                    if line.find('%init') > 0:
                        ok=False
                        if verbose>0:
                            print('')
                            print(file)
                            print(n0,line.rstrip())

                # detect critical-block end
                if line.startswith('!$omp end'):
                    watch=False
        if ok:
            if verbose>1:
                print(' OK:',file)
            return 0
        else:
            if verbose>1:
                print('NOK:',file)
            return 1

#%%
def find():
    n=0
    for root, dirs, files in os.walk("."):
        for file in files:
             if root != './.git' and root != './patch' and not root.endswith('temporary_files'):
                 if file.endswith('.f90'):
                    f=file if root=='.' else root.lstrip('.').lstrip(os.sep)+os.sep+file
                    f=f.strip()
                    try:
                        n+=critical(f,verbose=1)
                    except:
                        print('ERROR: '+f)
    if n==0:
        print('no problems found')

if __name__ == '__main__':
    find()
