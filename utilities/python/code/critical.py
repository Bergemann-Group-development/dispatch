# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:14:06 2019

@author: Aake
"""

import os

#%%
def critical(file,log1,log2,commit=False,verbose=0):
    """ Find and change omp critical (input_cr) """
    region='(cleanup_cr)'
    with open(file) as fo:
        lines=[]
        n0=0        # line number
        watch=False
        use_io=False
        critical=False
        filename=False
        if verbose>1:
            print('')
            print(file)
        for full_line in fo:
            n0+=1
            line=full_line.lstrip()
            
            if line.startswith('USE io_mod'):
                use_io=True

            # detect critical-block start
            if line.startswith('!$omp critical '+region):
                #print('')
                if not filename:
                    str='EDITED: '+file
                    print(str)
                    log1.write(str+'\n')
                    filename=True
                full_line=full_line.replace('!$omp critical '+region,
                                            'call io%lock%set')
                watch=True
                critical=True

            # detect 'end critical'
            elif watch:
                if line.startswith('!$omp end critical '+region):
                    full_line=full_line.replace('!$omp end critical '+region,
                                                'call io%lock%unset')
                    watch=False
            lines.append(full_line)

            # detect other crtical regions
            if line.startswith('!$omp critical') or \
               line.startswith('!$omp end critical'):
                if line.find('input_cr') < 0:
                    if not filename:
                        str='OTHER: '+file
                        print(str)
                        log2.write(str+'\n')
                        filename=True
                    str='{:4d} {}'.format(n0,full_line)
                    print(str.rstrip())
                    log2.write(str)

        if critical and commit:
            print("writing",file)
            with open(file,'w') as fo:
                for line in lines:
                    if not use_io:
                        if line.lstrip().startswith('USE'):
                            i=line.index('_mod')
                            fo.write('  USE io'+line[i:])
                            print("wrote USE io_mod")
                            use_io=True
                    fo.write(line)
        return 0

#%%
def find(log1,log2):
    for root, dirs, files in os.walk("."):
        for file in files:
             if root != './.git' and root != './patch' and not root.endswith('temporary_files'):
                 if file.endswith('.f90'):
                    f=file if root=='.' else root.lstrip('.').lstrip(os.sep)+os.sep+file
                    f=f.strip()
                    critical(f,log1,log2,verbose=1,commit=True)

if __name__ == '__main__':
    log1=open('critical_input.log','w')
    log2=open('critical_other.log','w')
    find(log1,log2)
    log1.close()
    log2.close()
