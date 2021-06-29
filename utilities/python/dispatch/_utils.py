import os
import time
import numpy as np
import f90nml
import time
import save_restore as sr
from scipy.io import FortranFile

def read_patch_metadata(io=0,run='',data='../data',n_rank=0,verbose=0):
    """ Read patch metadata from Fortran unformatted file
    """
    if verbose>0:
        start=time.time()

    if n_rank==0:
        dir='{}/{}/{:05d}'.format(data,run,io)
        files=[file for file in os.listdir(dir) if file.endswith('patches.nml')]
        n_rank=len(files)

    patches={}

    file0='{}/{}/{:05d}/patches.sr'.format(data,run,io)
    file1='{}/{}/{:05d}/rank_00000_patches.nml'.format(data,run,io)
    if os.path.exists(file0) and os.path.getmtime(file0) > os.path.getmtime(file1):
        patches=sr.restore(file0)
        if verbose>0:
            print(file0,len(patches.keys()),'patches')
    else:
        if verbose > 0:
            print('cashing metadata to',file0)
        rank=0
        while(rank < n_rank):
            file1='{}/{}/{:05d}/rank_{:05d}_patches.nml'.format(data,run,io,rank)
            read_rank_namelist(file1,patches,rank)
            if verbose>1:
                print(file1,len(patches.keys()),'patches')
            rank+=1
        try:
            sr.save(file0,patches)
        except:
            pass

    np=len(patches.keys())
    if verbose>0:
        sec=time.time()-start
        print('{} ranks, {} patches, {:.2f} sec'.format(n_rank,np,sec))

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

def read_rank_namelist(file,patches,rank,verbose=0):
    """
        Optimized parsing of patch namelist entries.  Updates dictionary with
        patch ID as key and a dictionary with multi propereties as the value.
    """
    n=0
    if verbose>1:
        start=time.time()
    with open(file,'r') as fo:
        nbor_block=False
        watch_block=False
        for line in fo:
            # strip commas and equal sign from line and split
            line=line.replace('=',' ').replace(',',' ').replace('"',' ')
            items=line.split()
            if len(items)==0:
                continue

            # Detect blocks to parse
            if items[0] == "&PATCH_NML":
                d={}
                watch_block=True
            elif items[0] == "&NBOR_NML":
                nbor_block=True
                if verbose>0:
                    print('starting nbor_block')

            # Parse a block of the namelist with an nbor list
            elif nbor_block:
                def parse_nbors(items,nbor_ids):
                    for item in items:
                        if item.find('*') >=0:
                            print('warning, multiple:',id,item)
                            rep=item.split('*')
                            for i in range(int(rep[0])):
                                nbor_ids.append(int(rep[1]))
                        else:
                            nbor_ids.append(int(item))
                if verbose>1:
                    print(items)
                if items[0] == 'ID':
                    id = int(items[1])
                    if verbose>0:
                        print('nbor_block id:',id)
                elif items[0] == 'PARENT_ID':
                    parent_id = int(items[1])
                    if verbose>0:
                        print('nbor_block parent_id:',parent_id)
                elif items[0] == "NBOR_IDS":
                    nbor_ids=[]
                    parse_nbors(items[1:],nbor_ids)
                    if verbose>0:
                        print('nbor_block nbor_ids:',nbor_ids)
                elif items[0] == "/":
                    nbor_ids=list(np.unique(np.sort(nbor_ids)))
                    patches[id]['nbor_ids']=nbor_ids
                    nbor_block=False
                    if verbose>0:
                        print('ending nbor_block:')
                else:
                    parse_nbors(items,nbor_ids)
                    if verbose>0:
                        print('appeded nbor_ids:',items)

            # Parse a block of the namelist with patch info
            elif watch_block:
                if len(items)==0:
                    pass
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
                elif items[0]=='RECORD':
                    d['record']=int(items[1])
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
