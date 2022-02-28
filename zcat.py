import glob
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from   astropy.io import fits as fits
from   astropy.table import Table, vstack, join
from   astropy import constants as const
from   astropy import units as u
from   astropy.table import QTable
from   desimodel.footprint import tiles2pix
from   get_fba_fromnewmtl import tileid2fbas
from   desitarget.sv3.sv3_targetmask import desi_mask

import os
import healpy as hp
import fitsio
import math

import sys

# sys.path.append(os.environ['HOME'] + '/LSS/py')

import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.geomask import get_imaging_maskbits 
from   desitarget.io import read_targets_in_tiles

zcatdatamodel =[('RA', '>f8'),\
                ('DEC', '>f8'),\
                ('TARGETID', '>i8'),
                ('NUMOBS', '>i4'),\
                ('Z', '>f8'),\
	        ('ZWARN', '>i8'),\
                ('ZTILEID', '>i4'),\
               ]

def coadd(zz):
    # HACK
    zz['TARGETID'][:10] = 1

    uids, cnts = np.unique(zz['TARGETID'], return_counts=True)

    repeats    = uids[cnts > 1]
    rcnts      = cnts[cnts > 1]

    # print(repeats)
    
    for tid, cnt in zip(repeats, rcnts):
        row    = zz[zz['TARGETID'] == tid][:1]
        rest   = zz[zz['TARGETID'] != tid]
        
        row['NUMOBS']      = cnt

        row['ZWARN']       = 0
        row['ZTILEID']     = rest['ZTILEID'].max()

        # TODO:  All good redshifts for now.
        # TODO:  All good redshifts for now; MTL does not track dChi2. 
        # row['DELTACHI2'] *= np.sqrt(cnt)

        zz                 = np.array(rest.tolist() + row.tolist(), dtype=zcatdatamodel)

    zz = zz[np.argsort(zz['TARGETID'])]
        
    return  zz 

def multitile2zcat(fs, zdir):
    print(zdir)

    for ff in fs:
        print(ff)

    zbests = [fba2zcat(fpath, zdir) for fpath in fs]
    zbests = [zbest.tolist() for zbest in zbests if zbest is not None]
    zbests = [item for sublist in zbests for item in sublist] 
    zbests = np.array(zbests, dtype=zcatdatamodel)

    zbests = coadd(zbests)

    zbests = Table(zbests)
    
    return zbests

def fba2zcat(fpath, zdir, nside=32):    
    # set values for mock zcat.
    # RA, DEC, ZTILEID, NUMOBS, DELTACHI2, ZWARN;
    fba = Table.read(fpath)
    tar = Table.read(fpath, 'FTARGETS')
    # tar.pprint()
    
    # TODO:  Why no SV3_DESI_TARGET etc?
    # bgs = (fba['SV3_DESI_TARGET'].data & desi_mask['BGS_ANY']) != 0     
    # bgs = fba[bgs]
    bgs = fba
    # bgs.pprint()
    
    tid            = fpath.split('-')[-1].split('.')[0]

    fbadir, _      = tileid2fbas(tid)
    ts             = str(tid).zfill(6)
    tpath          = fbadir+ts+'-tiles.fits'
    
    tiles          = Table.read(tpath)
    # tiles.pprint()

    pix            = tiles2pix(nside, tiles=tiles, radius=None, per_tile=False, fact=2**7)
        
    # Force zcat data model.
    zz                  = np.zeros(len(bgs), dtype=zcatdatamodel)
    zz['RA']            = bgs['TARGET_RA']
    zz['DEC']           = bgs['TARGET_DEC']
    zz['TARGETID']      = bgs['TARGETID']
    zz['NUMOBS']        = 1
    zz['ZWARN']         = 0
    zz['ZTILEID']       = tid
    
    zz                  = zz[np.argsort(zz['TARGETID'])]

    # Assume anything assigned is redshifted.                                                                                                                                                        
    pix_paths           = ['{}/sv3zs-bright-hp-{}.ecsv'.format(zdir, x) for x in pix]
    pix_paths           = [pix_path for pix_path in pix_paths if os.path.isfile(pix_path)]

    # print(pix_paths)
    
    mock_zs             = vstack([Table.read(pix_path) for pix_path in pix_paths])
    mock_zs             = mock_zs[np.isin(mock_zs['TARGETID'].data, zz['TARGETID'])]
    
    mock_zs.sort('TARGETID')
    
    # mock_zs.pprint()
    
    zz                  = zz[np.isin(zz['TARGETID'], mock_zs['TARGETID'])]    
    zz['Z']             = mock_zs['Z']
    
    # TODO:
    # https://github.com/desihub/desitarget/blob/1c7edc091ba7b8c628914826abcd5ee9c7a8bf24/py/desitarget/mtl.py#L2394
    
    '''
    https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L410

    - Sources in the zcat with `ZWARN` of `BAD_SPECQA|BAD_PETALQA|NODATA`
      are always ignored.
    - As sources with `BAD_SPECQA|BAD_PETALQA|NODATA` are ignored they
      will have ridiculous z-entries if `trimtozcat`=``False`` is passed!
    - The input `zcat` MAY BE MODIFIED. If a desideratum is that `zcat`
      remains unaltered, make sure to copy `zcat` before passing it.
    '''

    if len(zz) == 0:
        return None

    else:
        return  zz


if __name__ == '__main__':
    fpath   = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial//Univ000//fa/SV3/20210406/fba-000201.fits'
    zdir    = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/zs/'
    
    zz      = fba2zcat(fpath, zdir)

    # print(zz)

    
    fpaths  = glob.glob('/global/cscratch1/sd/mjwilson/altmtls/ledger/initial//Univ000//fa/SV3/20210406/fba*.fits')

    zz      = multitile2zcat(fpaths, zdir)
    
    print(zz)
    
