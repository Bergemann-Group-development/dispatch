import os
import time
import numpy as np
from scipy.io import FortranFile

def read_patch_metadata(io=0,run='',data='../data',n_rank=1,verbose=0):
    """ Read patch metadata from Fortran unformatted file
    """
    if verbose>0:
        start=time.time()
    
    rank=0
    patches={}
    while(rank < n_rank):
        file0='{}/{}/{:05d}/rank_{:05d}_patches.dat'.format(data,run,io,rank)
        file1='{}/{}/{:05d}/rank_{:05d}_patches.nml'.format(data,run,io,rank)
        if os.path.exists(file0):
            read_rank_metadata(file0,patches,rank)
        elif os.path.exists(file1):
            read_rank_namelist(file1,patches,rank)
        else:
            break
        if verbose>1:
            print(metadata,len(patches.keys()),'patches')
        rank+=1
    
    np=len(patches.keys())
    if verbose>0:
        sec=time.time()-start
        print('{} patches, {:.2f} sec'.format(n,sec))

    if np==0:
        print('no files found')
        return None
    else:
        return patches
            
def _parse_multiple (items):
    """ Figure out namelist records that may have N*xxx syntax
    """
    pos=items[1:]
    if len(pos)==1:
        pos[0]=pos[0].replace('3*','')
        pos=[pos[0],pos[0],pos[0]]
    elif len(pos)==2:
        if pos[0].find('2*')>=0:
            pos[0]=pos[0].replace('2*','')
            pos=[pos[0],pos[0],pos[1]]
        if pos[1].find('2*')>=0:
            pos[1]=pos[1].replace('2*','')
            pos=[pos[0],pos[1],pos[1]]
    # choose right data type for the resulting array
    try:
        pos = [int(p) for p in pos]
    except:
        pos = [float(p) for p in pos]
    return np.array(pos)

def read_nbor_list(snap,d):
    if hasattr(snap,'params_list'):
        params_list=snap.params_list
        if 'io_params' in params_list:
            io_params=params_list['io_params']
            if 'nml_version' in io_params:
                nml_version=io_params['nml_version']

    # Inform about nbor_nml, if present
    if nml_version>1:
        patches_nml=f90nml.read(file)
        if 'nbor_nml' in patches_nml:
            start=time()
            nbor_nml=patches_nml['nbor_nml']
            nbor_info={}
            for nb in nbor_nml:
                id=nb['id']
                parent_id=nb['parent_id']
                nbor_ids=nb['nbor_ids']
                nbor_info[id]={'parent_id':parent_id,'nbor_ids':nbor_ids}
            if verbose:
                print('snapshot nbor_nml info ({:.3f}s):'.format(time()-start))
                print(nbor_nml[0].keys())

        if nbor_info:
            d['parent_id']=nbor_info[id]['parent_id']
            d['nbor_ids']=nbor_info[id]['nbor_ids']

def read_rank_namelist(file,patches,rank,verbose=0):
    """
        Optimized parsing of patch namelist entries.  Updates dictionary with
        patch ID as key and a dictionary with multi propereties as the value.
    """
    n=0
    if verbose>1:
        start=time.time()
    with open(file,'r') as fo:
        watch_block = False
        for line in fo:
            # strip commas and equal sign from line and split
            line=line.replace('=',' ').replace(',',' ').replace('"',' ')
            items=line.split()
            if len(items)==0:
                pass
            elif items[0] == "&PATCH_NML":
                d={}
                watch_block = True
            elif items[0]=='ID':
                id=int(items[1])
            elif items[0]=='POSITION':
                d['position'] = _parse_multiple(items)
            elif items[0]=='SIZE':
                d['size'] = _parse_multiple(items)
            elif items[0]=='LEVEL':
                d['level']=int(items[1])
            elif items[0]=='DTIME':
                d['dtime']=float(items[1])
            elif items[0]=='TIME':
                d['time']=float(items[1])
            elif items[0]=='ISTEP':
                d['istep']=int(items[1])
            elif items[0]=='DS':
                d['ds'] = _parse_multiple(items)
            elif items[0]=='MESH_TYPE':
                d['mesh_type'] = int(items[1])
            elif items[0]=='NCELL':
                d['ncell'] = _parse_multiple(items)
            elif items[0]=='KIND':
                d['kind']=items[1]
            elif items[0]=='/' and watch_block: # the final entry of a namelist is always "/"
                d['rank']=rank
                patches[id]=d
                watch_block = False
                n+=1
                if verbose>2:
                    print('id',id)
    if verbose>1:
        used=time.time()-start
        print('{}: {} patches, {:.2f} sec'.format(file,n,used))

def read_rank_metadata(file,patches,rank,verbose=0):
    """ Read patch metadata from unformatted file """
    n=0
    if verbose>1:
        start=time.time()
    with FortranFile(file) as fd:
        key=fd.read_ints()[0]
        while(1):
            try:
                id,t,p,s,l,dt,tt,dtt,o=fd.read_record('i4,f8,(3)f8,(3)f8,i4,f8,(2)f8,(2)f8,i4')[0]
                patches[id]=dict(time=t,position=p,size=s,level=l,dtime=dt,offset=o,rank=rank)
                n+=1
                if verbose>2:
                    print('id',id)
            except:
                break
    if verbose>1:
        used=time.time()-start
        print('{}: {} patches, {:.2f} sec'.format(file,n,used))
