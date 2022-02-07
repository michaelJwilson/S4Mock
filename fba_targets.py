import os
import sys

from   astropy.table import Table
from   desimodel.footprint import tiles2pix

sys.path.append(os.environ['HOME'] + '/LSS/py/')

from   LSS.SV3.cattools import tile2rosette
from   desimodel.footprint import tiles2pix


root  = os.environ['CSCRATCH'] + '/altmtl/'

# ---- Get SV3 tiles and assign rosette ---
sv3_tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
sv3_tiles = sv3_tiles[sv3_tiles['PROGRAM'] == 'BRIGHT']                                                                                                                                              

# https://github.com/desihub/LSS/blob/092c7d98ad50b3a70f8d4e1d2b148fb36e21bb14/py/LSS/SV3/cattools.py#L21
ros       = [tile2rosette(x) for x in sv3_tiles['TILEID'].data]

sv3_tiles['ROS'] = ros
sv3_tiles.pprint()

# ----  Limit tiles to rosettes: G12: [1,2]; G15: [8,9,10, 17] ----
sv3_tiles = sv3_tiles[sv3_tiles['ROS'] == 1]
sv3_tiles.pprint()

'''                                                                                                                                                                                                
----  For each tile remaining, calculating the pix list at the desitarget nside  ---- 
      pixarea2nside(7.) = 32? TBC.
 
TILEID    RA     DEC   OBSCONDITIONS IN_DESI PROGRAM                                                                                                                                                 
int32  float64 float64     int32      int16   bytes6                                                                                                                                                
------ ------- ------- ------------- ------- -------                                                                                                                                                
    37 179.481  -0.018            15       1     SV3                                                                                                                                                                                                                                                                                                                                                '''
pix = tiles2pix(32, tiles=sv3_tiles, radius=None, per_tile=True, fact=2**7)

print('\n\n')

for tt, pp in zip(sv3_tiles['TILEID'], pix):
    print(tt, pp)

sv3_tiles = sv3_tiles[:2]
sv3_tiles.pprint()

sv3_tiles.write('test-tiles.fits', format='fits', overwrite=True)
    
'''
For every BGS Bright galaxy, assign a hp at the same nside.   By doing so, 
get a target list for each tile in the above rosettes.

Limit to minimal columns:
['RA', 'DEC', 'OBSCONDITIONS', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET', 'TARGETID', 'NUMOBS_MORE', 'PRIORITY', 'SUBPRIORITY']

assigning modal values, as before. 

write to 
/global/cscratch1/sd/mjwilson/altmtl/targets/

starting with one tile to check, with tileid 000037.
'''
'''
tpath = root + '/original/000037-targ.fits'
targ  = Table.read(tpath)

for x in sorted(targ.dtype.names):
    print(x)

# Complete list of columns needed by fba_run (could be less still). 
# https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/io.py#L53
min_cols = ['RA', 'DEC', 'OBSCONDITIONS', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET', 'TARGETID', 'NUMOBS_MORE', 'PRIORITY', 'SUBPRIORITY']
targ     = targ[min_cols]
targ.pprint()

targ.write(root + '/000037-targ.fits', format='fits', overwrite=True)
'''
print('\n\nDone.\n\n')
