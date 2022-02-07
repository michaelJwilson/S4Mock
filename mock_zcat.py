import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits as fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from astropy.table import QTable

import os
import healpy as hp
import math

import sys

sys.path.append(os.environ['HOME'] + '/LSS/py')

import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.geomask import get_imaging_maskbits 


def create_mock_zcat_hp(fpath='/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/bright_v0.9.fits', healpix=2286, nside=32):    
    npix = hp.nside2npix(nside)
    pixel_area = hp.nside2pixarea(nside, degrees=True)
    
    f = fits.open(fpath)
    mxxl=f[1].data

    theta = np.pi / 2. - np.radians(mxxl['DEC'].data)
    phi = np.radians(mxxl['RA'].data)

    #indices of pixels with non-zero density, unorganised list.
    all_pixel_indices = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)

    #indice of filled pixels and corrosponding targets in pixel
    filled_pixel_index, filled_targets_per_pixel = np.unique(all_pixel_indices, return_counts=True) 

    #no. targets per pixel, initially 0 
    targets_per_pixel = np.zeros(hp.nside2npix(nside))

    #update no. targets per pixel 
    targets_per_pixel[filled_pixel_index] = filled_targets_per_pixel/pixel_area
    targets_per_pixel[targets_per_pixel == 0] = np.NaN 
    
    #########################
    # cut to a single pixel
    
    single_mask = (all_pixel_indices==healpix)
    single_pixel_mxxl = mxxl[single_mask]
    single_pixel_mxxl = Table(single_pixel_mxxl)
    
    #########################

    # set values for mock zcat
    zcatdatamodel = np.array([], dtype=[
                                       ('RA', '>f8'),\
                                       ('DEC', '>f8'),\
                                       ('TARGETID', '>i8'),
                                       ('NUMOBS', '>i4'),\
                                       ('Z', '>f8'),\
                                       ('ZWARN', '>i8'),\
                                       ('ZTILEID', '>i4')
                                       ])
    t = Table(mtldatamodel) 

    prev_maxtid=0 
    ztileid = 24 #int

    for i, row in enumerate(single_pixel_mxxl):
        t.add_row((row['RA'],\
                   row['DEC'],\
                   prev_maxtid,\
                   0,\
                   row['Z']\
                   0,\
                   ztileid))

    t.meta['AUTHOR']  = 'Leah Bigwood' 
    t.meta['Mock']    = True 
    
    return t 
