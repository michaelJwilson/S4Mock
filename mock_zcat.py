import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from   astropy.io import fits as fits
from   astropy.table import Table
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


def fba2zcat(fpath='/global/cscratch1/sd/mjwilson/altmtls/fba-000037.fits'):    
    # set values for mock zcat.
    # RA, DEC, ZTILEID, NUMOBS, DELTACHI2, ZWARN;
    zcatdatamodel = np.array([], dtype=[
                                       ('RA', '>f8'),\
                                       ('DEC', '>f8'),\
                                       ('TARGETID', '>i8'),
                                       ('NUMOBS', '>i4'),\
                                       ('Z', '>f8'),\
                                       ('ZWARN', '>i8'),\
                                       ('ZTILEID', '>i4')
                                       ])
    zz  = Table(zcatdatamodel) 

    fba = Table.read(fpath) 
    tar = Table.read(fpath, 'FTARGETS')
    # tar.pprint()
    
    # bgs = (fba['SV3_DESI_TARGET'].data & desi_mask['BGS_ANY']) != 0     
    # bgs = fba[bgs]
    bgs = fba
    bgs.pprint()

    # FIBER      TARGETID      LOCATION FIBERSTATUS LAMBDA_REF PETAL_LOC DEVICE_LOC DEVICE_TYPE     TARGET_RA           TARGET_DEC           FA_TARGET      FA_TYPE FIBERASSIGN_X FIBERASSIGN_Y
    
    exit(0)
    '''
    for i, row in enumerate(single_pixel_mxxl):
        t.add_row((row['RA'],\
                   row['DEC'],\
                   prev_maxtid,\
                   0,\
                   row['Z']\
                   0,\
                   ztileid))
    '''
    t.meta['AUTHOR']  = 'Leah Bigwood' 
    t.meta['Mock']    = True 
    
    return t 

if __name__ == '__main__':
    fba2zcat()
