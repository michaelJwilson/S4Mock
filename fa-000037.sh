#!/bin/bash

source /global/project/projectdirs/desi/software/desi_environment.sh master

export SKYBRICKS_DIR=${DESI_ROOT}/target/skybricks/v2

module swap fiberassign/2.4.0

fba_run --targets /global/cscratch1/sd/mjwilson/altmtls/000037-targ.fits /global/cscratch1/sd/mjwilson/altmtls/000037-scnd.fits --sky /global/cscratch1/sd/mjwilson/altmtls/000037-sky.fits --footprint /global/cscratch1/sd/mjwilson/altmtls/000037-tiles.fits --rundate 2021-05-07T19:16:32+00:00 --fieldrot 0.000817750726657575 --dir /global/cscratch1/sd/mjwilson/altmtls/ --sky_per_petal 40 --standards_per_petal 10 --overwrite --sky_per_slitblock 1
