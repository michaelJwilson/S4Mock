#functions to help with running fiberassign, using SV3 parameters/targets

import fitsio
import numpy as np
from astropy.table import Table,join
# system
import os
import subprocess
import sys
import tempfile
import shutil
import re

# time
from time import time
from datetime import datetime, timedelta
from altcreate_mtl import altcreate_mtl

#import some functions from fiberassign
#from fiberassign.assign import minimal_target_columns
#from fiberassign.fba_launch_io import (
#    mv_temp2final,
#    force_finite_pm,
#    force_nonzero_refepoch,
#    gaia_ref_epochs,
#    mv_write_targets_out
#)

#from desitarget
import desitarget
from desitarget import io 

#hardcode target directories; these are fixed

skydir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/skies'
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'

# AR default REF_EPOCH for PMRA=PMDEC=REF_EPOCH=0 objects
gaia_ref_epochs = {"dr2": 2015.5}


minimal_target_columns= ['RELEASE','BRICKNAME','BRICKID','BRICK_OBJID','MORPHTYPE','RA',\
'DEC','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2','FLUX_IVAR_G','FLUX_IVAR_R',\
'FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z',\
'FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','REF_EPOCH','MASKBITS','SERSIC',\
'SHAPE_R','SHAPE_E1','SHAPE_E2','REF_ID','REF_CAT','GAIA_PHOT_G_MEAN_MAG',\
'GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG','PARALLAX','PMRA','PMDEC','PHOTSYS',\
'TARGETID','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT','NUMOBS_INIT','SV3_DESI_TARGET',\
'SV3_BGS_TARGET','SV3_MWS_TARGET','SV3_SCND_TARGET']

def tileid2fbas(tileid):
    ts = str(tileid).zfill(6)
    #get info from origin fiberassign file                                                                                                                                                            
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    indir = fht['OUTDIR']
    if fht['DESIROOT'] == '/data/datasystems':
        indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/' +fht['PMTIME'][:10].translate({ord('-'): None})  +'/'
        try:
            f = fitsio.read(indir+ts+'-targ.fits')
        except:

            date = int(fht['PMTIME'][:10].translate({ord('-'): None}))-1
            indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/'+str(date)+'/'

    return indir, fht 


def get_fba_fromnewmtl(tileid,mtldir=None,getosubp=False,outdir=None,faver=None):
    ts = str(tileid).zfill(6)

    indir, fht = tileid2fbas(tileid)

    print(indir)        
    tilef = indir+ts+'-tiles.fits'
    try:
        fitsio.read(tilef)
    except:
        return('Error! tile file does not appear to exist for tile '+ts+' '+tilef)    
    skyf = indir+ts+'-sky.fits'
    try:
        fitsio.read(skyf)
    except:
        print('Error! sky file does not appear to exist')    
    scndf = indir+ts+'-scnd.fits'
    scnd = True 
    try:
        fitsio.read(scndf)
    except:
        print(' secondary file does not appear to exist')
        scnd = False 
    gfaf = indir+ts+'-gfa.fits'
    try:
        fitsio.read(gfaf)
    except:
        print('Error! gfa file does not appear to exist')    
    toof = indir+ts+'-too.fits'
    too = os.path.isfile(toof)
    if too:
        print('will be using too file '+toof)
    if outdir is None:
        outdir = '/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/SV3rerun/'
    if getosubp == True or mtldir == None:
        outdir += 'orig/'
    if mtldir == None:
        tarfn = indir+ts+'-targ.fits' 
    else:
        tarfn = outdir+ts+'-targ.fits'    
    prog = fht['FAPRGRM'].lower()
    gaiadr = None
    if np.isin('gaiadr2',fht['FAARGS'].split()):
        gaiadr = 'dr2'
    if np.isin('gaiaedr3',fht['FAARGS'].split()):
        gaiadr = 'edr3'

    print('Ready to create alt. mtl.')
                
    if mtldir is not None:
        altcreate_mtl(tilef,
        mtldir+prog,        
        gaiadr,
        fht['PMCORR'],
        tarfn,
        tdir+prog)
        
    if getosubp:
        otar = Table.read(indir+ts+'-targ.fits')
        otar.keep_columns(['TARGETID','SUBPRIORITY'])
        ntar = Table.read(tarfn)
        ntar.remove_columns(['SUBPRIORITY'])
        ntar = join(ntar,otar,keys=['TARGETID'])
        ntar.write(tarfn,format='fits', overwrite=True)

    print('Ready to write: {}'.format(outdir+'fa-'+ts+'.sh'))
    
    fo = open(outdir+'fa-'+ts+'.sh','w')
    fo.write('#!/bin/bash\n\n')
    fo.write('source /global/project/projectdirs/desi/software/desi_environment.sh master\n')

    if faver == None:
        faver = float(fht['FA_VER'][:3])
        if faver == 2.4:
            fo.write('export SKYBRICKS_DIR=${DESI_ROOT}/target/skybricks/v2\n')

        if faver < 2.4:
            if int(indir[-7:-1]) > 210413:
                fo.write("module swap fiberassign/2.3.0\n") #inspection of results revealed tiles that used 2.2.dev* after 20210413 are reproduced using 2.3.0 and those before using 2.2.0
            else:
                fo.write("module swap fiberassign/"+fht['FA_VER'][:3]+'.0'+"\n")
        else:
            fo.write("module swap fiberassign/"+fht['FA_VER']+"\n")
    else:
        fo.write("module swap fiberassign/"+str(faver)+"\n")
        faver = float(faver[:3])
    fo.write("fba_run")
    fo.write(" --targets "+tarfn)
    if scnd:
        fo.write(" "+scndf)
    if too:
        fo.write(" "+toof)
    fo.write(" --sky "+skyf)
    fo.write(" --footprint "+tilef)
    rundate= fht['RUNDATE']
    if rundate == '2021-04-10T21:28:37':
        rundate = '2021-04-10T20:00:00'
    fo.write(" --rundate "+rundate)
    fo.write(" --fieldrot "+str(fht['FIELDROT']))
    fo.write(" --dir "+outdir)
    fo.write(" --sky_per_petal 40 --standards_per_petal 10")
    #fo.write(" --by_tile true")
    if faver >= 2.4:
        fo.write(" --sky_per_slitblock 1")
    if faver >= 3:
        fo.write(" --ha "+str(fht['FA_HA']))
        fo.write(" --margin-gfa 0.4 --margin-petal 0.4 --margin-pos 0.05")
    fo.close()    

    return  'fba-' + ts + '.fits'
    
#     if float(fht['FA_VER'][:3]) < 2.4:
#         fo.write("module swap fiberassign/2.3.0\n")
#     else:
#         fo.write("module swap fiberassign/"+fht['FA_VER']+"\n")
#     fo.write("fba_run")
#     fo.write(" --targets "+tarfn+" "+scndf)
#     if too:
#         fo.write(" "+toof)
#     fo.write(" --sky "+skyf)
#     fo.write(" --footprint "+tilef)
#     fo.write(" --rundate "+fht['RUNDATE'])
#     fo.write(" --fieldrot "+str(fht['FIELDROT']))
#     fo.write(" --dir "+outdir)
#     #fo.write(" --by_tile true")
#     if float(fht['FA_VER'][:3]) >= 3:
#         fo.write(" --ha "+str(fht['FA_HA']))
#     fo.close()    
