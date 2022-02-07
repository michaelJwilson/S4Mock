from astropy.table  import Table
from desitarget.io  import read_targets_in_tiles
from desitarget.mtl import tiles_to_be_processed, inflate_ledger


# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L293
# '/global/cfs/cdirs/desi/target/mtl/mtl-done-tiles.ecsv'
survey       = 'sv3'
obscon       = 'BRIGHT'
altmtltilefn = '/global/cscratch1/sd/mjwilson/altmtls/mtl-done-tiles.ecsv'
zcatdir      = '/global/cscratch1/sd/mjwilson/altmtls/barbaro/'

# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L314
# /global/cfs/cdirs/desi/spectro/redux/daily/tiles-specstatus.ecsv

# Alt. mtl dir has its own mtl tiles file:
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L342

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

# Checkpoint 'G'
# Create targets files input to fiberassign from ledger.  
# get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)

# altcreate_mtl
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/fatools.py#L375
# Array of tiles in the desimodel format
# If mtl:  ``True`` then the `columns and `header` kwargs are ignored and the full ledger is returned.

'''
isodate : :class:`str`, defaults to ``None``
  Only used if `mtl` is ``True`` An ISO date, such as returned by
  :func:`desitarget.mtl.get_utc_date() `. The ledger is restricted
  to entries strictly BEFORE `isodate` before being extracted.
  If ``None`` is passed then no date restrictions are applied.

unique : :class:`bool`, optional, defaults to ``True``
        If ``True`` then only read targets with unique `TARGETID` from
        MTL ledgers. Only used if `mtl` is ``True``.

        Returns the last entry. 
'''

mtldir  = '/global/cscratch1/sd/mjwilson/altmtls/iledger/'
tiles   = Table.read('/global/homes/m/mjwilson/desi/S4MOCK/LEAHFORK/S4Mock/test-tiles.fits')
mtltime = None

targ    = read_targets_in_tiles(mtldir, tiles, quick=False, mtl=True, unique=True, isodate=mtltime)

# Checkpoint 'H'
# Run fiberassign with subprocess.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L410

# Make redshift catalog.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L417

# Per-fiber info. swapping. 
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L447

# Update the ledger.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L453
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L1777
# update_ledger(althpdirname, altZCat, obscon=obscon.upper(), numobs_from_ledger=numobs_from_ledger)

# Write ledger.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L459
