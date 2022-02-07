#!/bin/bash

source /global/project/projectdirs/desi/software/desi_environment.sh master

module swap fiberassign/2.2.0

fba_run --targets /global/cscratch1/sd/mjwilson/altmtls/000039-targ.fits /global/cscratch1/sd/mjwilson/altmtls//000039-scnd.fits --sky /global/cscratch1/sd/mjwilson/altmtls//000039-sky.fits --footprint /global/cscratch1/sd/mjwilson/altmtls//000039-tiles.fits --rundate 2021-04-06T00:39:37 --fieldrot 0.000560328685810459 --dir /global/cscratch1/sd/mjwilson/altmtls/ --sky_per_petal 40 --standards_per_petal 10

# --overwrite
