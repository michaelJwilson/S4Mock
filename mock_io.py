import time
import glob
import desitarget
import numpy as np

from   ros_tools import tile2rosette
from   astropy.table import Table, vstack

'''
def read_targets():
  # What's the original path.
  # select bright targets.
  # calculate mags and ebv correction.
  to_grab=glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/targets-bright-hp-*.fits') 

  # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
  to_grab = sorted(to_grab) 
  
  hp_stack = []

  #do timer as takes a while
  start = time.time() 

  #total number of pixels, not quite sure where this has come from as npix is less than this above 
  mmask = 'BGS_TARGET' # NOTE: SV3 targets.  
  ttype = 'BGS_BRIGHT' 

  #loop through pixels
  for i, x in enumerate(to_grab):
      x = fits.open(x)
      f = np.array(x[1].data)[['RA','DEC','TARGETID', 'BGS_TARGET', 'MWS_TARGET','PHOTSYS']]
      #mask for bgs objects
      
      is_bgs = (f[mmask] & bgs_mask[ttype]) != 0
      #idx = np.arange(len(x))[is_bgs]
      #x = x.iloc[idx] 
      hp_stack.append(f[is_bgs])
      
      #more timing stuff
      if (i % 20) == 0:
          runtime = (time.time() - start)
          
        print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

  data_stack = np.concatenate(hp_stack)      
        
  data_stack = Table(data_stack)
  mask,idx = np.unique(data_stack['TARGETID'],return_index=True)
  data_stack = data_stack[idx]

  # MW_TRANSMISSION_GRZ, EBV, 
  # data_stack['RMAG_DRED'] = 22.5 - 2.5 * np.log10(flux / mwtrans)
  
  
def read_ledgers(mock=True):
    # desitarget.io.function
    # altmtls/ledger/
    if mock:
        root = ''
    else:
        root = ''

    # '/global/cscratch1/sd/mjwilson/S4MOCK/SV3REAL/SV3REALLEDGER/bright/*'
    fpaths = sorted(glob.glob(root + '*.'))
   
    return vstack([Table.read(x) for x in fpaths])

def read_desitargetrandoms(number=1, mask=None):
    ns = np.arange(number)

    fpaths = [f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-{nn}-0.fits' for nn in ns]

    rand = vstack([Table.read(x) for x in fpaths])

    if mask is not None:
        # use mask.
        raise NotImplementedError('Stop.')
       
    return rand 
'''
def read_sv3tiles():
    tiles = Table.read('/global/cscratch1/sd/mjwilson/S4MOCK/tiles-sv3.ecsv')
    tiles = tiles[(tiles['STATUS'] == 'done') & (tiles['PROGRAM']=='BRIGHT')]
    tiles['ROSETTE'] = np.array([tile2rosette(x) for x in tiles['TILEID']])
    return tiles
