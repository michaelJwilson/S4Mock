# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L293
# '/global/cfs/cdirs/desi/target/mtl/mtl-done-tiles.ecsv'

# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L314
# /global/cfs/cdirs/desi/spectro/redux/daily/tiles-specstatus.ecsv

# Alt. mtl dir has its own mtl tiles file:
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L342

# Build the name of an output target file (or directory)
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/io.py#L2502
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L344
althpdirname = io.find_target_files(altmtldir, flavor="mtl", resolve=resolve, survey=survey, obscon=obscon, ender=form)

# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L353
# Checkpoint 'A'
tiles = tiles_to_be_processed(zcatdir, altmtltilefn, obscon, survey)

# Checkpoint 'B' - limit to single tileid
# Checkpoint 'C' - limit to single zdate.

# Redshifts processed on date. 
tiles = tiles_to_be_processed(zcatdir, altmtltilefn, obscon, survey)

# Checkpoint 'G'
# Create targets files input to fiberassign from ledger.  
get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)

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
update_ledger(althpdirname, altZCat, obscon=obscon.upper(), numobs_from_ledger=numobs_from_ledger)

# Write ledger.
# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/altmtltools.py#L459
