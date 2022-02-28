# Get rid of the libraries that aren't required, e.g. pyplot.
import os
import sys
import math
import argparse
import healpy as hp
import pandas as pd
import numpy  as np
import matplotlib.gridspec as gridspec

from   astropy.io import fits as fits
from   astropy.table import Table
from   astropy import constants as const
from   astropy import units as u
from   astropy.table import QTable
from   geometry import radec2pix
from   mock_io import read_sv3tiles 

sys.path.append(os.environ['HOME'] + '/LSS/py')

import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools

from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask
from   desitarget.geomask import get_imaging_maskbits 
from   LSS.SV3.cattools import tile2rosette


def load_mxxl(nside=32, subsample=1):
    fpath = '/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/bright_v0.9.fits'
    mxxl  = Table.read(fpath)

    mxxl  = mxxl[::subsample]
    
    theta = np.pi / 2. - np.radians(mxxl['DEC'].data)
    phi   = np.radians(mxxl['RA'].data)

    mxxl['HPX'] = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)

    # single_pixel_mxxl['BGS_BRIGHT'] = single_pixel_mxxl['RMAG_DRED'] <= 19.5
    # single_pixel_mxxl['BGS_FAINT']  = (single_pixel_mxxl['RMAG_DRED'] > 19.5) & (single_pixel_mxxl['RMAG_DRED'] <= 20.175)

    return  mxxl
    
def create_mock_ledger_hp(outdir, healpix=2286, nside=32, mxxl=None, overwrite=False):    
    # TODO: Check nside matches desitarget file split NSIDE.     
    if mxxl == None:
        mxxl = load_mxxl()

    # E.g. /global/cfs/cdirs/desi/survey/catalogs/SV3/LSS//altmtl/debug_jl/alt_mtls_run128/Univ000/sv3/bright/sv3mtl-bright-hp-2286.ecsv
    opath = outdir + '/sv3mtl-bright-hp-{:d}.ecsv'.format(healpix)

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
    is_faint =  single_pixel_mxxl['BGS_FAINT']   == True
    
    for x in ['PRIORITY', 'PRIORITY_INIT','BGS_TARGET','DESI_TARGET']:
        single_pixel_mxxl[x] = -99

    # TODO:  BGS_TARGET, DESI_TARGET -> SV3_BGS_TARGET; SV3_DESI_TARGET.

    # HACK
    # single_pixel_mxxl.rename_column('BGS_TARGET',  'SV3_BGS_TARGET')
    # single_pixel_mxxl.rename_column('DESI_TARGET', 'SV3_DESI_TARGET')
    # single_pixel_mxxl.rename_column('MWS_TARGET',  'SV3_MWS_TARGET')

    # HACK
    del single_pixel_mxxl['BGS_TARGET']
    del	single_pixel_mxxl['DESI_TARGET']

    single_pixel_mxxl['SV3_MWS_TARGET']  = 0
    single_pixel_mxxl['SV3_BGS_TARGET']  = 0
    single_pixel_mxxl['SV3_DESI_TARGET'] = 0
    
    #bright columns using modal values, for initial ledger. 
    single_pixel_mxxl['PRIORITY_INIT'][is_bright]   = 102100
    single_pixel_mxxl['PRIORITY'][is_bright]        = 102100
    single_pixel_mxxl['SV3_BGS_TARGET'][is_bright]  = 514
    single_pixel_mxxl['SV3_DESI_TARGET'][is_bright] = 1152921504606846976 

    #faint columns using modal values 
    single_pixel_mxxl['PRIORITY_INIT'][is_faint]   = 102000
    single_pixel_mxxl['PRIORITY'][is_faint]        = 102000
    single_pixel_mxxl['SV3_BGS_TARGET'][is_faint]  = 257
    single_pixel_mxxl['SV3_DESI_TARGET'][is_faint] = 1152921504606846976
    
    #promote the faint higher priority ones i.e 20\% of faints
    draws    = np.random.uniform(0, 1, len(single_pixel_mxxl))
    is_hip   = (draws > 0.8) & is_faint

    single_pixel_mxxl['PRIORITY_INIT'][is_hip]     = 102100
    single_pixel_mxxl['PRIORITY'][is_hip]          = 102100
    single_pixel_mxxl['SV3_DESI_TARGET'][is_hip]   = 1152921504606846976 
    single_pixel_mxxl['SV3_BGS_TARGET'][is_hip]    = 265 

    print('Check: {:.3f}'.format(np.mean(is_hip)))
    
    # TODO: Check with BGS_TARGET & bgs_mask ... what flags the model values satisfy.
    
    #########################
    
    #create ledger 
    
    mtldatamodel = np.array([], dtype=[('RA', '>f8'),\
                                       ('DEC', '>f8'),\
                                       ('PARALLAX', '>f4'),\
                                       ('PMRA', '>f4'),\
                                       ('PMDEC', '>f4'),\
                                       ('REF_EPOCH', '>f4'),\
                                       ('SV3_DESI_TARGET', '>i8'),\
                                       ('SV3_BGS_TARGET', '>i8'),\
                                       ('SV3_MWS_TARGET', '>i8'),\
                                       ('TARGETID', '>i8'),\
                                       ('SUBPRIORITY', '>f8'),\
                                       ('OBSCONDITIONS', 'i4'),\
                                       ('PRIORITY_INIT', '>i8'),\
                                       ('NUMOBS_INIT', '>i8'),\
                                       ('PRIORITY', '>i8'),\
                                       ('NUMOBS', '>i8'),\
                                       ('NUMOBS_MORE', '>i8'),\
                                       ('Z', '>f8'),\
                                       ('ZWARN', '>i8'),\
                                       ('TIMESTAMP', 'U25'),\
                                       ('VERSION', 'U14'),\
                                       ('TARGET_STATE', 'U30'),\
                                       ('ZTILEID', '>i4'),\
                                       ('SV3_SCND_TARGET', '>i8')]) 

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
                   row['SV3_DESI_TARGET'],\
                   row['SV3_BGS_TARGET'],\
                   0,\
                   prev_maxtid,\
                   row['SUBPRIORITY'],\
                   516,\
                   row['PRIORITY_INIT'],\
                   3,\
                   row['PRIORITY'],\
                   0,\
                   3,\
                   row['Z'],\
                   -1,\
                   '2021-04-04T23:05:09',\
                   '0.57.0',\
                   'BGS|UNOBS',\
                   -1,\
                   0))

    t.meta['ISMOCK']     = 1 
    t.meta['SURVEY']     = 'sv3'
    t.meta['OBSCON']     = 'BRIGHT'
    t.meta['FILENSID']   = 32
    
    # t.meta['OVERRIDE'] = 'False'

    # HACK
    t['TARGETID']        = 1000 * healpix + np.arange(len(t))
    
    print(f'Writing {opath}')

    o = Table(t, copy=True)
    o['Z'] = -1.
    
    o.write(opath, format='ascii.ecsv', overwrite=overwrite)

    z = t['TARGETID', 'Z']

    opath = os.path.dirname(opath) + '/' + os.path.basename(opath).replace('mtl', 'zs')

    print(f'Writing {opath}')
    
    t['TARGETID', 'Z'].write(opath, format='ascii.ecsv', overwrite=overwrite)
    
    return  0


if __name__ == '__main__':
    # python mock_ledger.py --healpixel 1 --nside 32
    parser    = argparse.ArgumentParser(description='Create mock ledger for a given healpixel.')
    parser.add_argument('--healpixel',  type=int, default=None, help='Healpixel.')
    parser.add_argument('--nside',      type=int, default=32,   help='nside.')
    parser.add_argument('--overwrite',  help='Overwrite existing files', action='store_true')
    parser.add_argument('--outdir',     type=str, help='Output directory.', default='/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/')
    
    args      = parser.parse_args()
    hpixel    = args.healpixel
    nside     = args.nside
    overwrite = args.overwrite
    outdir    = args.outdir 

    tiles     = read_sv3tiles()
    hps       = radec2pix(tiles['RA'].data, tiles['DEC'].data, unique=True)

    tiles.pprint()
    
    if hpixel is not None:
        hps  += hpixel

    print(hps)

    mxxl      = load_mxxl(subsample=1)
    
    for ii in hps:
        try:
            create_mock_ledger_hp(outdir, healpix=ii, nside=nside, mxxl=mxxl, overwrite=overwrite)

        except Exception as E:
            print('ERROR on HPX {}'.format(ii))
            print(E)
            
    print('\n\nDone.\n\n')
