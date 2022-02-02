import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits as fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from astropy.table import QTable

import healpy as hp
import math

import sys
sys.path.append('/global/homes/l/lbigwood/LSS/py')
import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from desitarget.geomask import get_imaging_maskbits 

nside = 32
orig_density_per_deg = 2500

npix = hp.nside2npix(nside)
pixel_area = hp.nside2pixarea(nside,degrees=True)

def create_mock_ledger_hp(file = '/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/bright_v0.9.fits',healpix=2286):
    
    #map of all healpix in MXXL file
    
    f = fits.open(file)
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
    
    #cut to a single pixel
    
    single_mask = (all_pixel_indices==healpix)
    single_pixel_mxxl = mxxl[single_mask]
    single_pixel_mxxl = Table(single_pixel_mxxl)
    
    #########################
    
    #set values for mock ledger 
    
    #true/false array for bright/faint objects
    single_pixel_mxxl['BGS_BRIGHT'] = single_pixel_mxxl['RMAG_DRED'] <= 19.5
    
    #set subpriorities for all 
    single_pixel_mxxl['SUBPRIORITY'] = np.random.uniform(0, 1, len(single_pixel_mxxl))
    
    #set some other headings for all
    for x in ['PARALLAX', 'PMRA', 'PMDEC', 'REF_EPOCH']:
        single_pixel_mxxl[x] = 0.0
        
    #mask for brights
    is_bright =  single_pixel_mxxl['BGS_BRIGHT'] == True
    #mask for faints
    is_faint =  single_pixel_mxxl['BGS_BRIGHT'] == False
    
    for x in ['PRIORITY', 'PRIORITY_INIT','BGS_TARGET','DESI_TARGET']:
        single_pixel_mxxl[x] = -99

    #bright columns using modal values
    single_pixel_mxxl['PRIORITY_INIT'][is_bright] = 102100
    single_pixel_mxxl['PRIORITY'][is_bright] = 102100
    single_pixel_mxxl['BGS_TARGET'][is_bright] = 514
    single_pixel_mxxl['DESI_TARGET'][is_bright] = 1152921504606846976 

    #faint columns using modal values 
    single_pixel_mxxl['PRIORITY_INIT'][is_faint] = 102000
    single_pixel_mxxl['PRIORITY'][is_faint] = 102000
    single_pixel_mxxl['BGS_TARGET'][is_faint] = 257
    single_pixel_mxxl['DESI_TARGET'][is_faint] = 1152921504606846976
    
    #promote the faint higher priority ones i.e 20\% of faints
    draws    = np.random.uniform(0, 1, len(single_pixel_mxxl))
    is_hip   = (draws > 0.8) & is_faint

    single_pixel_mxxl['PRIORITY_INIT'][is_hip] = 102100
    single_pixel_mxxl['PRIORITY'][is_hip] = 102100
    single_pixel_mxxl['DESI_TARGET'][is_hip] = 1152921504606846976 
    single_pixel_mxxl['BGS_TARGET'][is_hip] = 265 
    
    #########################
    
    #create ledger 
    
    mtldatamodel = np.array([], dtype=[ 

        ('RA', '>f8'), ('DEC', '>f8'), ('PARALLAX', '>f4'), 

        ('PMRA', '>f4'), ('PMDEC', '>f4'), ('REF_EPOCH', '>f4'), 

        ('DESI_TARGET', '>i8'), ('BGS_TARGET', '>i8'), ('MWS_TARGET', '>i8'), 

        ('SCND_TARGET', '>i8'), ('TARGETID', '>i8'), 

        ('SUBPRIORITY', '>f8'), ('OBSCONDITIONS', 'i4'), 

        ('PRIORITY_INIT', '>i8'), ('NUMOBS_INIT', '>i8'), ('PRIORITY', '>i8'), 

        ('NUMOBS', '>i8'), ('NUMOBS_MORE', '>i8'), ('Z', '>f8'), ('ZWARN', '>i8'), 

        ('TIMESTAMP', 'U25'), ('VERSION', 'U14'), ('TARGET_STATE', 'U30'), 

        ('ZTILEID', '>i4') 

        ]) 


    t = Table(mtldatamodel) 

    # Entries correspond to the datamodel above.  
    # RA and DEC are first two entries, need replaced by the mock value.  
    # TARGETID needs to start at 0 and increment by 1 with every add row.  
    # SUBPRIORITY is a column with values equivalent to np.uniform(0, 1, len(mxxl_healpixel)) 
    # PRIORITY_INIT = 102100 for BGS BRIGHT, 102000 FOR BGS FAINT.  
    # PRIORITY = PRIORITY_INIT  

    #be careful with target ids overlapping for faint and bright targets???
    prev_maxtid=0 

    for i, row in enumerate(single_pixel_mxxl):
        t.add_row((row['RA'],\
                   row['DEC'],\
                   row['PARALLAX'],\
                   row['PMRA'],\
                   row['PMDEC'],\
                   row['REF_EPOCH'],\
                   row['DESI_TARGET'],\
                   row['BGS_TARGET'],\
                   0,\
                   0,\
                   prev_maxtid,\
                   row['SUBPRIORITY'],\
                   516,\
                   row['PRIORITY_INIT'],\
                   9,\
                   row['PRIORITY'],\
                   0,\
                   9,\
                   -1.0,\
                   -1,\
                   '2021-04-04T23:05:09',\
                   '0.57.0',\
                   'BGS|UNOBS',\
                   -1))

    t.meta['AUTHOR']  = 'Leah Bigwood' 
    t.meta['Mock']    = True 
    
    return t 
