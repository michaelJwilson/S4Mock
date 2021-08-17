import astropy.io.fits as fits

from   astropy.table import Table

# https://docs.astropy.org/en/stable/table/index.html

fpath = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/v0.1/LRGAlltiles_full.dat.fits'
# fpath = '/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/dark/LRG_oneperztrue_clus.dat.fits'

dat = Table.read(fpath)
dat.pprint()

for x in dat.dtype.names:
    print(x)
