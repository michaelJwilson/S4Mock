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

import os
import healpy as hp
import fitsio
import math

import sys

sys.path.append(os.environ['HOME'] + '/LSS/py')


import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.io import read_targets_in_tiles


def fba2zcat(fpath='/global/cscratch1/sd/mjwilson/altmtls/fba-000039.fits', nside=32):    
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
    bgs.pprint()
    
    zz['RA']       = bgs['TARGET_RA']
    zz['DEC']      = bgs['TARGET_DEC']
    zz['TARGETID'] = bgs['TARGETID']

    zz             = zz[~np.isnan(zz['RA'].data)]
    tid            = fpath.split('-')[-1].split('.')[0]
    
    tiles          = Table.read(os.path.dirname(fpath) + '/{}-tiles.fits'.format(tid))
    # tiles.pprint()

    pix            = tiles2pix(nside, tiles=tiles, radius=None, per_tile=False, fact=2**7)

    # Assume anything assigned is redshifted.
    all_zs              = vstack([Table.read('/global/cscratch1/sd/mjwilson/altmtls/ledger/zs/sv3zs-bright-hp-{}.ecsv'.format(x)) for x in pix])
    all_zs['NUMOBS']    = np.ones(len(all_zs), dtype='>i4')
    all_zs['ZWARN']     = np.zeros(len(all_zs), dtype='>i8')
    all_zs['ZTILEID']   = np.int(tid) * np.ones(len(all_zs), dtype='>i4')
    all_zs['DELTACHI2'] = 100.

    # all_zs.pprint()    
    zz                  = Table(zz)
    # zz.pprint()

    del zz['Z']
    del zz['NUMOBS']
    del zz['ZWARN']
    del zz['ZTILEID']

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
    zz = fba2zcat()
    zz = zz[zz['Z'].data > 0.1]
    
    zz.pprint()
