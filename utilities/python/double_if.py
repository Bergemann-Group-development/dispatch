# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:14:06 2019

@author: Aake
"""

import os
import sys

#%%
def find_double_if(file,verbose=0):
    """ Find and check legality of double if """
    with open(file) as fo:
        ifs=0       # if-level
        n0=0        # line number
        n1=0        # lines in if-level 1
        n2=0        # line with 'first_time ='
        ok=True
        printed=False
        if verbose:
            print('')
            print(file)
        for line in fo:
            n0+=1
            words=line.split()
            if len(words)>1:

                # detect if-block start
                if line.lstrip().startswith('if (first_time)'):
                    ifs+=1
                    if ifs>1:
                        #print(file,'has double-if')
                        n1=0
                        if verbose>0:
                            print('starting count at line',n0,ifs)
                            sys.stdout.write(line)
                elif ifs>1 and line.lstrip().startswith('if'):
                        ifs+=1

                # detect 'first_time ='
                if ifs>1:
                    n1+=1
                    if line.lstrip().startswith('first_time ='):
                        n2=n1
                        if verbose>0:
                            print('fist_time set at line',n0,ifs)
                            sys.stdout.write(line)

                # detect if-block end
                elif ifs>0 and line.lstrip().startswith('end if') \
                          or line.lstrip().startswith('endif'):
                    ifs-=1
                    if verbose>0:
                        print('detected end if at line',n0,ifs)
                        sys.stdout.write(line)
                    if ifs==1:
                        # detect cases where 'fist_tim =' is not last line
                        if n1>n2+1:
                            if not printed:
                                print(file)
                                printed=True
                            ok=False
                            print('    problem at line',n0,n1,n2)
            if ok:
                if verbose:
                    print(' OK:',file)
                return 0
            else:
                if verbose:
                    print('NOK:',file)
                return 1

#%%
def find():
    n=0
    for root, dirs, files in os.walk("."):
        for file in files:
             if file.endswith('.f90'):
                f=file if root=='.' else root.lstrip('.\\')+os.sep+file
                n+=find_double_if(f)
    if n==0:
        print('no problems found')

if __name__ == '__main__':
    find()