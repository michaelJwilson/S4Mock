import  glob
 
from    astropy.table   import Table
from    desitarget.io   import read_targets_in_tiles, read_mtl_ledger
from    desitarget.mtl  import tiles_to_be_processed, inflate_ledger
from    altcreate_mtl   import altcreate_mtl
from    mock_zcat       import fba2zcat
from    update_ledger   import update_ledger


# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L293
# ADM grab the MTL directory (in case we're relying on $MTL_DIR)
# mtldir = get_mtl_dir(mtldir)
# '/global/cfs/cdirs/desi/target/mtl/mtl-done-tiles.ecsv'
survey       = 'sv3'
obscon       = 'BRIGHT'
altmtltilefn = '/global/cscratch1/sd/mjwilson/altmtls/mtl-done-tiles.ecsv'
zcatdir      = '/global/cscratch1/sd/mjwilson/altmtls/barbaro/'

# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L314
# ADM And contruct the associated ZTILE filename.
# ztilefn = os.path.join(zcatdir, get_ztile_file_name())
# /global/cfs/cdirs/desi/spectro/redux/daily/tiles-specstatus.ecsv

# Alt. mtl dir has its own mtl tiles file:
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L342
# altmtltilefn = os.path.join(altmtldir, get_mtl_tile_file_name(secondary=secondary))

# Build the name of an output target file (or directory)
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/io.py#L2502
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L344
# althpdirname = io.find_target_files(altmtldir, flavor="mtl", resolve=resolve, survey=survey, obscon=obscon, ender=form)

# Checkpoint 'A'
# Redshifts processed on date.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L353

# ZDONE (and ARCHIVE_DATE > 0 for main) in tiles-specstatus.
# Diffed against list of done tiles /global/cscratch1/sd/mjwilson/altmtls/mtl-done-tiles.ecsv  
# **  newtiles = tiles_to_be_processed(zcatdir, altmtltilefn, obscon, survey, reprocess=False)  **
# print(newtiles)

# Checkpoint 'B' - limit to single tileid
# Checkpoint 'C' - limit to single zdate.
mtldir  = '/global/cscratch1/sd/mjwilson/altmtls/ledger/'
tpath   = '/global/homes/m/mjwilson/desi/S4MOCK/LEAHFORK/S4Mock/test-tiles.fits'
tiles   = Table.read(tpath)
mtltime = None

'''
ledgers = sorted(glob.glob(mtldir + '*.ecsv'))

for led in ledgers:
    read_mtl_ledger(led, unique=True, isodate=None, initial=False, leq=False)

    print(f'Fetched {led}')
'''
'''
targ    = read_targets_in_tiles(mtldir, tiles, quick=False, mtl=True, unique=True, isodate=mtltime)
print(targ)

# altcreate_mtl                                                                                                                                                                                     
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/fatools.py#L375
# If mtl:  ``True`` then the `columns and `header` kwargs are ignored and the full ledger is returned.                                                                                               
'''
'''
altmtl = altcreate_mtl(tpath,\
                       mtldir,\
                       'dr2',\
                       'n',\
                       '/global/cscratch1/sd/mjwilson/altmtls/{:06d}-targ.fits'.format(39),\
                       None,\
                       survey='sv3',\
                       mtltime=None,\
                       pmtime_utc_str=None,\
                       add_plate_cols=True)
'''
'''
'''
'''                                                                                                                                                                                                  
isodate : :class:`str`, defaults to ``None``                                                                                                                                                         
  Only used if `mtl` is ``True`` An ISO date, such as returned by                                                                                                                                    
  :func:`desitarget.mtl.get_utc_date() `. The ledger is restricted                                                                                                                                   
  to entries strictly BEFORE `isodate` before being extracted.                                                                                                                                       
  If ``None`` is passed then no date restrictions are applied.                                                                                                                                                                                                                                                                                                                                            
unique : :class:`bool`, optional, defaults to ``True``                                                                                                                                              
        If ``True`` then only read targets with unique `TARGETID` from                                                                                                                               
        MTL ledgers. Only used if `mtl` is ``True``.                                                                                                                                                         Returns the last entry. 
'''
# Checkpoint 'G'
# Create targets files input to fiberassign from ledger.  
# get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)

# altcreate_mtl
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/fatools.py#L375
# Array of tiles in the desimodel format
# If mtl:  ``True`` then the `columns and `header` kwargs are ignored and the full ledger is returned.

# Checkpoint 'H'
# Run fiberassign with subprocess.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L410

# /global/cscratch1/sd/mjwilson/altmtls/fa-000037.sh
# /global/cscratch1/sd/mjwilson/altmtls/fba-000037.fits

# Make redshift catalog.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L417
# make_zcat
# Full path to the "daily" directory that hosts redshift catalogs.
# Numpy array of tiles to be processed. Must contain at least:
#        * TILEID - unique tile identifier.
#        * ZDATE - final night processed to complete the tile (YYYYMMDD).
#
# /zcatdir/tiles/cumulative/
#
# zbestfns = sorted(glob(os.path.join(ymdir, "zbest*")))
#
# ADM recover the information for unique targets based on the
# ADM first entry for each TARGETID.
#
# zcat datamodel
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L2383
# 
# RA, DEC, ZTILEID, NUMOBS, DELTACHI2, ZWARN.
#
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L2393
#
# Per-fiber info. swapping. 
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L447

# Update the ledger.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L453
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L1777
althpdirname = '/global/cscratch1/sd/mjwilson/altmtls/ledger/'
altZCat      =  fba2zcat()

'''
hpdirname : :class:`str`
        Full path to a directory containing an MTL ledger that has been
        partitioned by HEALPixel (i.e. as made by `make_ledger`).
'''

'''
cat : :class:`~astropy.table.Table`, optional
        Redshift catalog table with columns ``TARGETID``, ``NUMOBS``,
        ``Z``, ``ZWARN``, ``ZTILEID``, and ``msaddcols`` at the top of
        the code for the Main Survey.
'''

'''
targets : :class:`~numpy.array` or `~astropy.table.Table`, optional
        A numpy rec array or astropy Table with at least the columns
        ``RA``, ``DEC``, ``TARGETID``, ``DESI_TARGET``, ``NUMOBS_INIT``,
        and ``PRIORITY_INIT``. If ``None``, then assume the `zcat`
        includes ``RA`` and ``DEC`` and look up `targets` in the ledger.
'''

print(althpdirname)

# Returns:  Nothing, but relevant ledger files are updated.
update_ledger(althpdirname, altZCat, obscon='BRIGHT', numobs_from_ledger=True)

# Write ledger.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L459
