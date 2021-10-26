from LSS.SV3.altmtltools import initializeAlternateMTLs, loop_alt_ledger
from astropy.table import Table

obscon='BRIGHT'
survey='sv3'

initMTL_path = '/global/cscratch1/sd/lbigwood/S4MOCK/orig_mtls/bright/sv3mtl-bright-hp-0.ecsv'
# initMTL = Table.read(initMTL_path)
# initMTL.pprint()

outputMTL_path = '/global/cscratch1/sd/lbigwood/S4MOCK/alt_mtls/{:d}/'
# outputMTL_path = outputMTL_path.format(1)

# print(outputMTL_path)

# >> https://github.com/desihub/LSS/blob/b3e06635459771c622135a0090acd34f3a174ca3/py/LSS/SV3/altmtltools.py#L167 
#
#  Note:  writes to e.g.  /global/cscratch1/sd/lbigwood/S4MOCK/alt_mtls/0/sv3/bright/sv3mtl-bright-hp-0.ecsv
#

# initializeAlternateMTLs(initMTL_path, outputMTL_path, nAlt = 2, seed = 314159, obscon = 'BRIGHT')

alt_mtldir = '/global/cscratch1/sd/lbigwood/S4MOCK/alt_mtls/'

all_dates = ['20210405','20210406','20210407','20210408','20210409','20210410','20210411','20210412','20210413','20210414','20210415','20210416','20210417','20210418','20210419','20210420','20210428','20210429','20210430','20210501','20210502','20210503','20210504','20210505','20210506','20210507','20210508','20210509','20210510','20210511','20210512','20210513','20210517','20210529']

# print(all_dates)

single_date = all_dates[0]

althpdirname, altmtltilefn, ztilefn, tiles = loop_alt_ledger(obscon, survey=survey, zcatdir=None, mtldir=None,
                                                             altmtlbasedir=alt_mtldir, ndirs=2, numobs_from_ledger=True,
                                                             secondary=False, singletile = None, singleDate=single_date, debugOrig = False, getosubp = False)
