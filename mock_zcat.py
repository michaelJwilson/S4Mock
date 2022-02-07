import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from   astropy.io import fits as fits
from   astropy.table import Table, vstack
from   astropy import constants as const
from   astropy import units as u
from   astropy.table import QTable

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


def fba2zcat(fpaths=['/global/cscratch1/sd/mjwilson/altmtls/fba-000037.fits']):    
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
    
    fba = vstack([Table.read(x) for x in fpaths])
    tar = vstack([Table.read(x, 'FTARGETS') for x in fpaths])
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
    
    # TODO: update redshifts, ztileid.
    zz['Z']        = 0
    zz['ZTILEID']  = 0

    zz             = Table(zz)
    
    zz.meta['AUTHOR'] = 'Leah Bigwood' 
    zz.meta['MOCK']   = 1
    
    return  zz

if __name__ == '__main__':
    zz = fba2zcat()
    zz.pprint()
    
