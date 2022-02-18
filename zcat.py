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

def coadd(zz):
    uids, cnts = np.unique(zz['TARGETID'].data, return_counts=True)

    repeats    = uids[cnts > 1]
    rcnts      = cnts[cnts > 1]
    
    for tid, cnt in zip(repeats, rcnts):
        row    = zz[zz['TARGETID'] == tid][0] 
        rest   = zz[zz['TARGETID'] != tid]
        
        row['NUMOBS']     = cnt

        # TODO:  All good redshifts for now. 
        row['ZWARN']      = 0
        row['ZTILEID']    = rest['ZTILEID'].max()
        row['DELTACHI2'] *= np.sqrt(cnt)
        
        zz                = rest.add_row(row)

    return zz 

def multitile2zcat(fs, zdir):
    zz = vstack([fba2zcat(fpath, zdir) for fpath in fs])
    zz = coadd(zz)

    return zz
    
def fba2zcat(fpath, zdir, nside=32):    
    # set values for mock zcat.
    # RA, DEC, ZTILEID, NUMOBS, DELTACHI2, ZWARN;
    zcatdatamodel =[('RA', '>f8'),\
                    ('DEC', '>f8'),\
                    ('TARGETID', '>i8'),
                    ('NUMOBS', '>i4'),\
                    ('Z', '>f8'),\
                    ('ZWARN', '>i8'),\
                    ('ZTILEID', '>i4'),\
                   ]
    
    fba = Table.read(fpath)
    tar = Table.read(fpath, 'FTARGETS')
    # tar.pprint()

    zz  = np.zeros(len(fba), dtype=zcatdatamodel)
    
    # TODO:  Why no SV3_DESI_TARGET etc?
    # bgs = (fba['SV3_DESI_TARGET'].data & desi_mask['BGS_ANY']) != 0     
    # bgs = fba[bgs]
    bgs = fba
    # bgs.pprint()
    
    zz['RA']       = bgs['TARGET_RA']
    zz['DEC']      = bgs['TARGET_DEC']
    zz['TARGETID'] = bgs['TARGETID']

    zz             = zz[~np.isnan(zz['RA'].data)]
    tid            = fpath.split('-')[-1].split('.')[0]

    fbadir         = tileid2fbas(tid)
    ts             = str(tid).zfill(6)
    tpath          = fbadir+ts+'-tiles.fits'
    
    tiles          = Table.read(tpath)
    # tiles.pprint()

    pix            = tiles2pix(nside, tiles=tiles, radius=None, per_tile=False, fact=2**7)
    
    # HACK
    pix            = [6398, 6399, 6570, 6740, 6741, 6743, 6912, 6914]
    
    # Assume anything assigned is redshifted.
    mock_zs             = vstack([Table.read('{}/sv3zs-bright-hp-{}.ecsv'.format(zdir, x)) for x in pix])

    # TODO:
    all_zs              = mock_zs 

    all_zs['NUMOBS']    = np.ones(len(all_zs), dtype='>i4')
    all_zs['ZWARN']     = np.zeros(len(all_zs), dtype='>i8')
    all_zs['ZTILEID']   = int(tid) * np.ones(len(all_zs), dtype='>i4')
    all_zs['DELTACHI2'] = 100.

    # all_zs.pprint()

    # Force zcat data model.
    zz                  = Table(zz)
    # zz.pprint()

    del  zz['Z']
    del  zz['NUMOBS']
    del  zz['ZWARN']
    del  zz['ZTILEID']

    '''
    https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L410

    - Sources in the zcat with `ZWARN` of `BAD_SPECQA|BAD_PETALQA|NODATA`
      are always ignored.
    - As sources with `BAD_SPECQA|BAD_PETALQA|NODATA` are ignored they
      will have ridiculous z-entries if `trimtozcat`=``False`` is passed!
    - The input `zcat` MAY BE MODIFIED. If a desideratum is that `zcat`
      remains unaltered, make sure to copy `zcat` before passing it.
    '''

    
    zz                = join(zz, all_zs, keys='TARGETID', join_type='left')
    zz['Z']           = zz['Z'].data.astype('>f8')
        
    return  zz


if __name__ == '__main__':
    fpath  = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial//Univ000//fa/SV3/20210406/fba-000201.fits'
    zdir   = '/global/cscratch1/sd/mjwilson/altmtls/ledger/zs/'
    
    zz     = fba2zcat(fpath, zdir)
    zz     = zz[zz['Z'].data > 0.1]

    fpaths = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial//Univ000//fa/SV3/20210406/fba*.fits'

    zz     = multitile2zcat(fpaths, zdir)
    
    
