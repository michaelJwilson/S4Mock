import time
import glob
import desitarget
import numpy as np
from astropy.io import fits as fits
from desitarget.targets import desi_mask, bgs_mask, mws_mask 
import pandas as pd



from   astropy.table import Table, vstack


def read_mainsurvey_targets(pixlist=None):
    # What's the original path.
    # select bright targets.
    # calculate mags and ebv correction.
    root = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/'

    if pixlist == None:
        to_grab=glob.glob(root + '/targets-bright-hp-*.fits') 
    else:
        to_grab=[root + f'/targets-bright-hp-{pix}.fits' for pix in pixlist]
        
    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
    to_grab = sorted(to_grab) 

    hp_stack = []

    #do timer as takes a while
    start = time.time() 

    #total number of pixels, not quite sure where this has come from as npix is less than this above 
    mmask = 'BGS_TARGET' # NOTE: SV3 targets.  
    ttype = 'BGS_BRIGHT' 

    min_cols = ['RA','DEC','TARGETID', 'BGS_TARGET', 'MWS_TARGET','PHOTSYS']
    
    #loop through pixels
    for i, x in enumerate(to_grab):
        x = fits.open(x)
        # f = np.array(x[1].data)[]
        f = fitsio.read(x, columns=min_cols)
        
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
    # RMAG, GMAG, W1, ZMAG, RFIBMAG, 
    return data_stack
  
  
def read_mainsurvey_ledgers(mock=True):
    # desitarget.io.function
    # altmtls/ledger/
    if mock:
        to_grab=glob.glob('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/bright/mtl-bright-hp-*.ecsv') 

    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
        to_grab = sorted(to_grab) 

        hp_stack = []

        #do timer as takes a while
        start = time.time() 

        #total number of pixels, not quite sure where this has come from as npix is less than this above 
        npix_todo = 200000

        mmask = 'BGS_TARGET'
        ttype = 'BGS_BRIGHT'

        #loop through pixels
        for i, x in enumerate(to_grab):
            x = pd.read_csv(x, comment='#', delimiter='\s+', usecols=['RA', 'DEC', 'TARGETID', 'BGS_TARGET', 'MWS_TARGET'])

            #mask for bgs objects
            is_bgs = (x[mmask] & bgs_mask[ttype]) != 0
            idx = np.arange(len(x))[is_bgs]
            x = x.iloc[idx] 
            hp_stack.append(x)

            #more timing stuff
            if (i % 100) == 0:
                runtime = (time.time() - start)

                print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

            if i > npix_todo:
                break

        # Create a big table from the list of tables.  
        data_stack = pd.concat(hp_stack, ignore_index=True)

        #unique targets only and put it in right table format
        mask,idx = np.unique(data_stack['TARGETID'],return_index=True)
        data_stack = data_stack.iloc[idx]
        data_stack = Table.from_pandas(data_stack)

        #more timing stuff
        runtime = (time.time() - start)
        print('\n\nTotal runtime of {:.6f} seconds after {:d} pixels'.format(runtime, npix_todo))

        
    else:
        root = '/global/cscratch1/sd/mjwilson/S4MOCK/SV3REAL/SV3REALLEDGER/bright/'
        fpaths = sorted(glob.glob(root + '*.'))

        data_stack = vstack([Table.read(x) for x in fpaths])
    
    return data_stack

def read_desitargetrandoms(number=1, mask=None):
    ns = np.arange(number)

    fpaths = [f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-{nn}-0.fits' for nn in ns]

    rand = vstack([Table.read(x) for x in fpaths])

    if mask is not None:
        # use mask.
        raise NotImplementedError('Stop.')
       
    return rand 

def read_sv3tiles():
    tiles = Table.read('/global/cscratch1/sd/mjwilson/S4MOCK/tiles-sv3.ecsv')
    tiles = tiles[(tiles['STATUS'] == 'done') & (tiles['PROGRAM']=='BRIGHT')]
    tiles['ROSETTE'] = np.array([tile2rosette(x) for x in tiles['TILEID']])
    return tiles
