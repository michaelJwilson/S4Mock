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

from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.geomask import get_imaging_maskbits 
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

    # Assume anything assigned is redshifted.
    zz['NUMOBS']   = 1
    zz['ZWARN']    = 0

    tid            = fpath.split('-')[-1].split('.')[0]
    
    tiles          = Table.read(os.path.dirname(fpath) + '/{}-tiles.fits'.format(tid))
    pix            = tiles2pix(nside, tiles=tiles, radius=None, per_tile=False, fact=2**7)
    
    all_zs         = vstack([Table.read('/global/cscratch1/sd/mjwilson/altmtls/iledger/sv3zs-bright-hp-{}.ecsv'.format(x)) for x in pix])
    all_zs.pprint()
    
    zz                = Table(zz)

    del  zz ['Z']
    
    zz                = join(zz, all_zs, keys='TARGETID', join_type='left')
    zz['Z']           = zz['Z'].data.astype('>f8')
    zz['ZTILEID']     = tid
    zz['DELTACHI2']   = 100.
    
    zz.meta['AUTHOR'] = 'Leah Bigwood' 
    zz.meta['MOCK']   = 1
    
    return  zz


if __name__ == '__main__':
    zz = fba2zcat()
    zz.pprint()
    
