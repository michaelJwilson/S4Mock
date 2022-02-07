# Get rid of the libraries that aren't required, e.g. pyplot.
import os
import sys
import math
import argparse
import healpy as hp
import pandas as pd
import numpy  as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from   astropy.io import fits as fits
from   astropy.table import Table
from   astropy import constants as const
from   astropy import units as u
from   astropy.table import QTable

sys.path.append(os.environ['HOME'] + '/LSS/py')

import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.geomask import get_imaging_maskbits 


def load_mxxl(nside=32):
    fpath = '/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/bright_v0.9.fits'
    mxxl  = Table.read(fpath)

    theta = np.pi / 2. - np.radians(mxxl['DEC'].data)
    phi   = np.radians(mxxl['RA'].data)

    mxxl['HPX'] = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)
    
    return  mxxl
    
def create_mock_ledger_hp(outdir, healpix=2286, nside=32, mxxl=None, overwrite=False):    
    # TODO: Check nside matches desitarget file split NSIDE.     
    if mxxl == None:
        mxxl = load_mxxl()

    opath = outdir + '/testledger-{:06d}.ecsv'.format(healpix)

    if os.path.isfile(opath) & ~overwrite:
        print(f'Warning: {opath} exists; skipping.')

        return 0
        
    npix       = hp.nside2npix(nside)
    pixel_area = hp.nside2pixarea(nside,degrees=True)

    print('npix: {}; pixel_area: {} for nside: {}'.format(npix, pixel_area, nside))
    
    single_mask = (mxxl['HPX'].data == healpix)
    single_pixel_mxxl = mxxl[single_mask]
    
    #true/false array for bright/faint objects
    single_pixel_mxxl['BGS_BRIGHT'] = single_pixel_mxxl['RMAG_DRED'] <= 19.5
    single_pixel_mxxl['BGS_FAINT']  = (single_pixel_mxxl['RMAG_DRED'] > 19.5) & (single_pixel_mxxl['RMAG_DRED'] <= 20.175)
    
    print('Selected {:.3f} as BGS Bright'.format(np.mean(single_pixel_mxxl['BGS_BRIGHT'])))
    
    #TODO: what is the resulting target density. 
    
    #set subpriorities for all 
    single_pixel_mxxl['SUBPRIORITY'] = np.random.uniform(0, 1, len(single_pixel_mxxl))
    
    #set some other headings for all
    for x in ['PARALLAX', 'PMRA', 'PMDEC', 'REF_EPOCH']:
        single_pixel_mxxl[x] = 0.0
        
    #mask for brights
    is_bright =  single_pixel_mxxl['BGS_BRIGHT'] == True

    #mask for faints
    is_faint =  single_pixel_mxxl['BGS_FAINT']   == False
    
    for x in ['PRIORITY', 'PRIORITY_INIT','BGS_TARGET','DESI_TARGET']:
        single_pixel_mxxl[x] = -99

    # TODO:  BGS_TARGET, DESI_TARGET -> SV3_BGS_TARGET; SV3_DESI_TARGET.
        
    #bright columns using modal values, for initial ledger. 
    single_pixel_mxxl['PRIORITY_INIT'][is_bright] = 102100
    single_pixel_mxxl['PRIORITY'][is_bright]      = 102100
    single_pixel_mxxl['BGS_TARGET'][is_bright]    = 514
    single_pixel_mxxl['DESI_TARGET'][is_bright]   = 1152921504606846976 

    #faint columns using modal values 
    single_pixel_mxxl['PRIORITY_INIT'][is_faint]  = 102000
    single_pixel_mxxl['PRIORITY'][is_faint]       = 102000
    single_pixel_mxxl['BGS_TARGET'][is_faint]     = 257
    single_pixel_mxxl['DESI_TARGET'][is_faint]    = 1152921504606846976
    
    #promote the faint higher priority ones i.e 20\% of faints
    draws    = np.random.uniform(0, 1, len(single_pixel_mxxl))
    is_hip   = (draws > 0.8) & is_faint

    single_pixel_mxxl['PRIORITY_INIT'][is_hip]    = 102100
    single_pixel_mxxl['PRIORITY'][is_hip]         = 102100
    single_pixel_mxxl['DESI_TARGET'][is_hip]      = 1152921504606846976 
    single_pixel_mxxl['BGS_TARGET'][is_hip]       = 265 

    print('Check: {:.3f}'.format(np.mean(is_hip)))
    
    # TODO: Check with BGS_TARGET & bgs_mask ... what flags the model values satisfy.
    
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

    t.meta['AUTHOR']   = 'L. Bigwood' 
    t.meta['ISMOCK']   = 1 
    t.meta['SURVEY']   = 'SV3'
    t.meta['OBSCON']   = 'BRIGHT'
    # t.meta['OVERRIDE'] = 'False'
    
    print(f'Writing {opath}')

    # E.g. /global/cfs/cdirs/desi/survey/catalogs/SV3/LSS//altmtl/debug_jl/alt_mtls_run128/Univ000/sv3/bright/sv3mtl-bright-hp-2286.ecsv
    t.write(opath, format='ascii.ecsv', overwrite=overwrite)
    
    return  0


if __name__ == '__main__':
    # python mock_ledger.py --healpixel 1 --nside 32
    parser    = argparse.ArgumentParser(description='Create mock ledger for a given healpixel.')
    parser.add_argument('--healpixel',  type=int, default=2286, help='Healpixel.')
    parser.add_argument('--nside',      type=int, default=32,   help='nside.')
    parser.add_argument('--overwrite',  help='Overwrite existing files', action='store_true')
    parser.add_argument('--outdir',     type=str, help='Output directory.', required=True)
    
    args      = parser.parse_args()
    hpixel    = args.healpixel
    nside     = args.nside
    overwrite = args.overwrite
    outdir    = args.outdir 

    mxxl      = load_mxxl()

    hps       = [6399, 6570, 6741, 6743, 6912, 6914]
    hps      += [6398, 6399, 6570, 6740, 6741, 6743, 6912, 6914]
    
    for ii in hps:
        create_mock_ledger_hp(outdir, healpix=ii, nside=nside, mxxl=mxxl, overwrite=overwrite)

    print('\n\nDone.\n\n')
