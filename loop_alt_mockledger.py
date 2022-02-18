import os
import sys
import glob
import numpy as np
import desitarget
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import numpy.lib.recfunctions as rfn
import fitsio
import subprocess
import LSS.SV3.fatools as fatools

from   numpy import random as rand
from   desitarget import io 
from   desitarget import mtl
from   desitarget import mtl
from   desitarget.targets import initial_priority_numobs
from   desitarget.targetmask import obsconditions, obsmask
from   desitarget.targetmask import desi_mask
from   desitarget.mtl import get_mtl_dir, get_mtl_tile_file_name,get_mtl_ledger_format
from   astropy.table import Table,join,unique,vstack
from   desiutil.log import get_logger
from   desitarget.mtl import get_zcat_dir, get_ztile_file_name, tiles_to_be_processed
from   desitarget.mtl import make_zcat,survey_data_model
from   desitarget.targets import decode_targetid
from   LSS.SV3.altmtltools import *
from   get_fba_fromnewmtl import get_fba_fromnewmtl
from   zcat import multitile2zcat
from   update_ledger import update_ledger


log = get_logger()

print(desitarget.__file__)
print(mtl.__file__)
print(io.__file__)
print(pf.__file__)
print(fatools.__file__)

os.environ['DESIMODEL'] = '/global/common/software/desi/cori/desiconda/current/code/desimodel/master'

zcatdatamodel = np.array([], dtype=[
    ('RA', '>f8'), ('DEC', '>f8'), ('TARGETID', '>i8'),
    ('NUMOBS', '>i4'), ('Z', '>f8'), ('ZWARN', '>i8'), ('ZTILEID', '>i4')
    ])

mtltilefiledm = np.array([], dtype=[
    ('TILEID', '>i4'), ('TIMESTAMP', 'U25'),
    ('VERSION', 'U14'), ('PROGRAM', 'U6'), ('ZDATE', 'U8')
    ])

def loop_alt_mockledger(obscon,\
                        survey='main',\
                        zcatdir=None,\
                        mtldir=None,
                        altmtlbasedir=None,\
                        ndirs = 1,\
                        numobs_from_ledger=True,\
                        secondary=False,\
                        singletile = None,\
                        singleDate = None,\
                        debugOrig = False,\
                        getosubp = False,\
                        quickRestart = False,\
                        redoFA = False,\
                        multiproc = False,\
                        nproc = None):
    """
    Execute full MTL loop, including reading files, updating ledgers.

    Parameters
    ----------
    obscon : :class:`str`
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set when merging targets.
    survey : :class:`str`, optional, defaults to "main"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    zcatdir : :class:`str`, optional, defaults to ``None``
        Full path to the "daily" directory that hosts redshift catalogs.
        If this is ``None``, look up the redshift catalog directory from
        the $ZCAT_DIR environment variable.
    mtldir : :class:`str`, optional, defaults to ``None``
        Full path to the directory that hosts the MTL ledgers and the MTL
        tile file. If ``None``, then look up the MTL directory from the
        $MTL_DIR environment variable.
    altmtlbasedir : :class:`str`, optional, defaults to ``None``
        Formattable path to a directory that hosts alternate MTL ledgers  
        If ``None``, then look up the MTL directory from the
        $ALT_MTL_DIR environment variable. This will fail since that variable
        is not currently set in the desi code setup.
    ndirs : :class:`int`, optional, defaults to ``3``
        Number of alternate MTLs to process within altmtlbasedir
    numobs_from_ledger : :class:`bool`, optional, defaults to ``True``
        If ``True`` then inherit the number of observations so far from
        the ledger rather than expecting it to have a reasonable value
        in the `zcat.`
    secondary : :class:`bool`, optional, defaults to ``False``
        If ``True`` then process secondary targets instead of primaries
        for passed `survey` and `obscon`.
    quickRestart : :class:`bool`, optional, defaults to ``False``
        If ``True`` then copy original alternate MTLs from 
        altmtlbasedir/Univ*/survey/obscon/orig and  
    redoFA : :class:`bool`, optional, defaults to ``False``
        If ``True`` then automatically redo fiberassignment regardless of
        existence of fiberassign file in alternate fiberassign directory
    multiproc : :class:`bool`, optional, defaults to ``False``
        If ``True`` then run a single MTL update in a directory specified by
        nproc.
    nproc : :class:`int`, optional, defaults to None
        If multiproc is ``True`` this must be specified. Integer determines 
        directory of alternate MTLs to update.

    Returns
    -------
    :class:`str`
        The directory containing the ledger that was updated.
    :class:`str`
        The name of the MTL tile file that was updated.
    :class:`str`
        The name of the ZTILE file that was used to link TILEIDs to
        observing conditions and to determine if tiles were "done".
    :class:`~numpy.array`
        Information for the tiles that were processed.

    Notes
    -----
    - Assumes all of the relevant ledgers have already been made by,
      e.g., :func:`~LSS.SV3.altmtltools.initializeAlternateMTLs()`.
    """

    import multiprocessing as mp
    import logging

    logger=mp.log_to_stderr(logging.DEBUG)

    # ADM first grab all of the relevant files.
    # ADM grab the MTL directory (in case we're relying on $MTL_DIR).
    mtldir = get_mtl_dir(mtldir)

    print(f'MTL dir.: {mtldir}')
    
    # ADM construct the full path to the mtl tile file.
    # E.g. /global/cfs/cdirs/desi/target/mtl/mtl-done-tiles.ecsv 
    mtltilefn = os.path.join(mtldir, get_mtl_tile_file_name(secondary=secondary))

    print(f'MTL tile filename: {mtltilefn}')
    
    # ADM construct the relevant sub-directory for this survey and
    # ADM set of observing conditions..
    form    = get_mtl_ledger_format()
    resolve = True
    msg     = "running on {} ledger with obscon={} and survey={}"

    if secondary:
        log.info(msg.format("SECONDARY", obscon, survey))
        resolve = None

    else:
        log.info(msg.format("PRIMARY", obscon, survey))
        
    # ADM grab the zcat directory (in case we're relying on $ZCAT_DIR).
    zcatdir = get_zcat_dir(zcatdir)

    print(f'ZCAT dir.: {zcatdir}')
    
    # ADM And contruct the associated ZTILE filename.
    # E.g. /global/cfs/cdirs/desi/spectro/redux/daily/tiles-specstatus.ecsv  
    ztilefn = os.path.join(zcatdir, get_ztile_file_name())
    
    print(f'ZCAT tiles status: {ztilefn}')
        
    if altmtlbasedir is None:
        print('This will automatically find the alt mtl. dir in the future but fails now. Bye.')
        assert(0)
        
    iterloop = range(ndirs)

    for n in iterloop:
        print('')
        print('')
        print('')
        print(f'*******  ALT MTl. REALIZATION (DIRECTORY)  NUMBER  {n}  *******')
        print('')
        print('')
        print('')
        print('')

        altmtldir = altmtlbasedir + '/Univ{0:03d}/'.format(n)

        print(f'ALT MTL dir: {altmtldir}')
        
        # fn: filename    
        altmtltilefn = os.path.join(altmtldir, get_mtl_tile_file_name(secondary=secondary))

        print(f'ALT MTL tile status filename: {altmtltilefn}')
        
        althpdirname = io.find_target_files(altmtldir, flavor="mtl", resolve=resolve,
                                            survey=survey, obscon=obscon, ender=form)
        
        print(f'ALT MTL dir/survey/program/: {althpdirname}')
        
        # ADM grab an array of tiles that are yet to be processed.
        print(f'Directory for redshift information: {zcatdir}')
        print(f'Alt. MTL tile status path: {altmtltilefn}')

        print(survey, obscon)

        #  SETUP STEP:
        #  cp /global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv /global/cscratch1/sd/mjwilson/altmtls/ledger/initial//Univ000/
        
        # Find tiles that are "done" but aren't yet in the MTL tile record.
        # https://github.com/desihub/desitarget/blob/1c7edc091ba7b8c628914826abcd5ee9c7a8bf24/py/desitarget/mtl.py#L2076
        # using e.g. ZDONE == True in /global/cfs/cdirs/desi/spectro/redux/daily/tiles-specstatus.ecsv
        # Main survey must have been archived: tilelookup["ARCHIVEDATE"] > 0
        tiles = tiles_to_be_processed(zcatdir, altmtltilefn, obscon, survey)

        #     ADM look up the time.
        #     newtiles["TIMESTAMP"] = get_utc_date(survey=survey)
        
        #     ADM the final processed date for the redshifts.
        #     newtiles["ZDATE"] = tiles["LASTNIGHT"]
        # 
        #     ADM the date the tile was archived.
        #     newtiles["ARCHIVEDATE"] = tiles["ARCHIVEDATE"]
        
        print('Checkpoint A: found atl. mtl ledger & {} tiles to be processed:'.format(len(tiles)))
            
        # ADM stop if there are no tiles to process.
        if len(tiles) == 0:
            if n != ndirs - 1:
                # More realizations to be done. 
                continue
            else:
                # All done. 
                return althpdirname, mtltilefn, ztilefn, tiles
            
        if not (singletile is None):
            # Process a single tile, assumnig its in the tiles to be done.
            assert singletile in tiles['TILEID'], f'TILEID {singletile} is invalid.'
            
            tiles = tiles[tiles['TILEID'] == singletile]
            
        # sorttiles = np.sort(tiles, order = 'ZDATE')

        print('Checkpoint B: Applied single tile cut, if necessary')
        
        # Dates to be reduced.                                                                                                                                                                      
        dates = np.sort(np.unique(tiles['ZDATE']))

        print(f'Checkpoint A1:  Available dates to reduce: {dates}')

        if not singleDate == None:
            # Process a single date, assumnig its in the dates to be done. 
            assert singleDate in tiles['ZDATE'], f'{singleDate} is invalid.'
            
            tiles = tiles[tiles['ZDATE'] == singleDate]

        else:
            assert(0)

        print('Checkpoint C: Applied single zdate cut, if necessary')

        # Dates to be reduced. 
        dates = np.sort(np.unique(tiles['ZDATE']))

        print(f'Checkpoint C1:  Ready to reduce dates: {dates}')

        fba_paths = []
        
        for date in dates:
            print(f'\n\n  -------  Reducing {date}  -------')

            dateTiles = tiles[tiles['ZDATE'] == date]

            print('Tiles whose cumulative redshifts lastnight were ZDATE {} on SV3 MTL update History: {}'.format(date, dateTiles['TILEID']))
            
            OrigFAs = []
            AltFAs  = []
            
            for t in dateTiles:
                print('Checkpoint D: reducing TILEID {} from lastnight {}'.format(t['TILEID'], date))

                ts         = str(t['TILEID']).zfill(6)
                
                FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
                fhtOrig    = fitsio.read_header(FAOrigName)
                fadate     = fhtOrig['RUNDATE']

                # reformat rundate 
                fadate     = ''.join(fadate.split('T')[0].split('-'))
                
                # prepare alt. mtl. directory to write  to. 
                fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/'
                
                print(f'Checkpoint E:  Retrieved SV3 fiberassign file (and run date) & ready to create directory {fbadirbase}.')
                
                if getosubp:
                    fbadir    = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/orig/'
                    FAAltName = fbadir + '/fba-' + ts+ '.fits'

                    print(f'Retrieving original subpriorities from {FAAltName}')
                    
                else:
                    FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/fba-' + ts+ '.fits'
                    fbadir    = fbadirbase

                    print(f'Looking for fiberassign file {FAAltName}')

                fba_paths.append(FAAltName)
                    
                # Note: forced to redo fiberassign. 
                if redoFA or (not os.path.exists(FAAltName)):
                    print(f'Checkpoint F: Rerunning fiberassign for {FAAltName}')
                    '''
                    if getosubp:
                        redo_fba_fromorig(ts,outdir=fbadirbase)
                    else:
                        print
                    '''

                    # get_fba_fromnewmtl(ts, mtldir=altmtldir + survey.lower() + '/', outdir=fbadirbase, getosubp = getosubp)
                    # https://github.com/desihub/LSS/blob/62446e54b10a5b531ef3b495b8f2f775d476ac25/py/LSS/SV3/fatools.py#L244
                    # 
                    # if mtldir == None:
                    #   tarfn = indir+ts+'-targ.fits' 
                    # else:
                    #   tarfn = outdir+ts+'-targ.fits'

                    # https://github.com/desihub/LSS/blob/62446e54b10a5b531ef3b495b8f2f775d476ac25/py/LSS/SV3/fatools.py#L375
                    #
                    # if mtldir is not None:
                    #   altcreate_mtl(tilef,         # path to a tiles fits file (string)
                    #                 mtldir+prog,   # mtldir: folder with ledger files        
                    #                 gaiadr,        # Gaia dr ("dr2" or "edr3")
                    #                 fht['PMCORR'], # apply proper-motion correction? ("y" or "n")
                    #                 tarfn,         # fits file name to be written 
                    #                 tdir + prog)   # targ. dir.

                    # mtltime: MTL isodate (string formatted as yyyy-mm-ddThh:mm:ss+00:00); this needs be considered carefully for alt mtls
                    
                    ledg        = altmtldir + survey.lower() + '/'

                    tarfn       = fbadirbase + ts + '-targ.fits'
                    fba_bash    = fbadirbase + 'fa-' + ts + '.sh'
                    
                    print(f'Checkpoint pre-G: Ready to write {tarfn} & {fba_bash} based on\n{ledg}.')
                    
                    fba_path    = get_fba_fromnewmtl(ts, mtldir=ledg, outdir=fbadirbase, getosubp = getosubp)
                    
                    print(f'Checkpoint G: Written {tarfn} & {fba_bash} based on {ledg}')
                    
                    command_run = (['bash', fbadir + 'fa-' + ts + '.sh'])

                    result      = subprocess.run(command_run, capture_output=True)

                    # https://docs.python.org/3/library/subprocess.html
                    stde        = result.stde
                    stdo        = result.stdout
                    retcode     = result.returncode
                    
                    print(f'Checkpoitn H:  Ran fiberassign with result {retcode}')

                    if retcode:
                        print(stde)
                        print(stdo)
                        
                    OrigFAs.append(pf.open(FAOrigName)[1].data)
                    AltFAs.append(pf.open(FAAltName)[1].data)

            print(f'Checkpoint I:  Finished reducing date {date}.  Time to create a redshift catalog.')

            # Where to find the mock redshifts.
            # HACK 
            mock_zdir  = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/zs/' 
            zcat       = multitile2zcat(fba_paths, mock_zdir)
            
            opath      = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/zbest-{}.txt'.format(fadate)

            # print(opath)
            
            np.savetxt(opath, zcat)

            # ADM create the catalog of updated redshifts.
            # HACK
            # zbest_dir = os.path.dirname(opath) 
            # zcat      = make_zcat(zbest_dir, dateTiles, obscon, survey)
            
            # ADM insist that for an MTL loop with real observations, the zcat
            # ADM must conform to the data model. In particular, it must include
            # ADM ZTILEID, and other columns addes for the Main Survey. These
            # ADM columns may not be needed for non-ledger simulations.
            # ADM Note that the data model differs with survey type.

            print(f'Checkpoint I2:  Created redshift catalog, checking datamodel')

            zcatdm  = survey_data_model(zcatdatamodel, survey=survey)
            
            if zcat.dtype.descr != zcatdm.dtype.descr:
                msg = "zcat data model must be {} not {}!".format(
                    zcatdm.dtype.descr, zcat.dtype.descr)
                log.critical(msg)
                raise ValueError(msg)

            print(f'Checkpoint I3:  All good.')
            
            # ADM useful to know how many targets were updated.
            _, _, _, _, sky, _ = decode_targetid(zcat["TARGETID"])
            ntargs, nsky       = np.sum(sky == 0), np.sum(sky)
            
            msg    = "Update state for {} targets".format(ntargs)
            msg   += " (the zcats also contain {} skies with +ve TARGETIDs)".format(nsky)

            log.info(msg)
            
            A2RMap = {}
            R2AMap = {}

            print('Checkpoint J:  Applying fiber remapping.')

            for ofa, afa in zip(OrigFAs, AltFAs):
                # Mapping of Real2Alt and Alt2Real Targetid based on fiber.
                A2RMapTemp, R2AMapTemp = createFAmap(ofa, afa)
                
                A2RMap.update(A2RMapTemp)
                R2AMap.update(R2AMapTemp)

                print('Checkpoint K:  Created fiber remapping, making alternate Z cat.')

                # print(type(zcat))
                # print(zcat.dtype)

                # Take the data zcatalog, clone it and rewrite the TARGETIDs using real2alt.
                altZCat = makeAlternateZCat(zcat, R2AMap, A2RMap)

                print('Checkpoint K2:  Created alternate Z cat.')

            # HACK
            if not redoFA:
                altZCat = zcat
                
            # print(type(altZCat))
            # print(altZCat.dtype)
            
            # ADM update the appropriate ledger.
            print(f'Checkpoint L:  Updating the ledger at {althpdirname} with alternate ZCAT.')
            
            # Redshift catalog table with columns ``TARGETID``, ``NUMOBS``, ``Z``, ``ZWARN``, ``ZTILEID``, and ``msaddcols``
            # https://github.com/desihub/desitarget/blob/1c7edc091ba7b8c628914826abcd5ee9c7a8bf24/py/desitarget/mtl.py#L1777

            # ADM also ignore anything with NODATA set in ZWARN.
            # nodata = zcat["ZWARN"] & zwarn_mask["NODATA"] != 0
            '''
            # SB ignore targets that failed QA: ZWARN bits BAD_SPECQA|BAD_PETALQA
            # badqa = zcat["ZWARN"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
            # HACK
            update_ledger(althpdirname, altZCat, obscon=obscon.upper(),
                          numobs_from_ledger=numobs_from_ledger)
            '''
            '''
            if survey == "main":
                print('Adding TIMESTAMP after sleep(1).')
                
                sleep(1)
                
                tiles["TIMESTAMP"] = get_utc_date(survey=survey)
            '''
            print('Checkpoint M:  Finished with ledger update.')

            print(f'Checkpoint pre-N:  Ready to update {altmtltilefn}.')
                        
            io.write_mtl_tile_file(altmtltilefn, dateTiles)

            print('Checkpoint N:  Finished with writing mtl tile file.')
                
            # ADM for the main survey "holding pen" method, ensure the TIMESTAMP
            # ADM in the mtl-done-tiles file is always later than in the ledgers.
        
        # ADM write the processed tiles to the MTL tile file.
        # io.write_mtl_tile_file(altmtltilefn, tiles)

    return  althpdirname, altmtltilefn, ztilefn, tiles


if __name__ == '__main__':
    redoFA        = True
    singleDate    = 20210405
    altmtlbasedir = '/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/'
    
    loop_alt_mockledger('BRIGHT',\
                        survey='sv3',\
                        zcatdir=None,\
                        mtldir=None,
                        altmtlbasedir=altmtlbasedir,\
                        ndirs=1,\
                        numobs_from_ledger=True,\
                        secondary=False,\
                        singletile=None,\
                        singleDate=singleDate,\
                        debugOrig=False,\
                        getosubp=False,\
                        quickRestart=False,\
                        redoFA=redoFA,\
                        multiproc=False,\
                        nproc=None)
    
