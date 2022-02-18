import os
import sys
import fitsio

from   astropy.table import Table
from   desitarget import io 
from   desitarget.mtl import inflate_ledger

sys.path.append(os.environ['HOME'] + '/LSS/py/')

from LSS.SV3.fatools import force_nonzero_refepoch, force_finite_pm

minimal_target_columns=['RA', 'DEC', 'OBSCONDITIONS', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET', 'TARGETID', 'NUMOBS_MORE', 'PRIORITY', 'SUBPRIORITY']

gaia_ref_epochs = {"dr2": 2015.5}

def altcreate_mtl(
    tilesfn,
    mtldir,            
    gaiadr,
    pmcorr,
    outfn,
    targdir,
    survey='sv3',
    mtltime=None,#I think we will just want this to be the latest for the re/alt runs    
    pmtime_utc_str=None,
    add_plate_cols=True#,
    #tmpoutdir=tempfile.mkdtemp(),
):
    """
    Mostly copied from fiberassign.fba_launch_io.create_mtl
    Create a (primary or secondary) target fits file, based on MTL ledgers (and complementary columns from desitarget targets files).
    
    Args:
        tilesfn: path to a tiles fits file (string)
        mtldir: folder with ledger files        
        targdir: desitarget targets folder (or file name if secondary) (string)        
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        survey: should just be sv3
        mtltime: MTL isodate (string formatted as yyyy-mm-ddThh:mm:ss+00:00); this needs be considered carefully for alt mtls
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_targets will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        add_plate_cols (optional, defaults to True): adds a PLATE_RA and PLATE_DEC columns (boolean)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
        for sv3-backup, we remove BACKUP_BRIGHT targets.
        TBD : if secondary targets, we currently disable the inflate_ledger(), as it
                seems to not currently work.
                hence if secondary and pmcorr="y", the code will crash, as the 
                GAIA_ASTROMETRIC_EXCESS_NOISE column will be missing; though we do not
                expect this configuration to happen, so it should be fine for now.
        TBD: the PLATE_{RA,DEC,REF_EPOCH} columns currently simply are copy of RA,DEC,REF_EPOCH
        TBD:    but it prepares e.g. to add chromatic offsets.
        20210526 : implementation of using subpriority=False in write_targets
                    to avoid an over-writting of the SUBPRIORITY; AJR changed to True reproduce SV3
    """
    tiles = fitsio.read(tilesfn)

    # AR mtl: read mtl
    d = io.read_targets_in_tiles(
        mtldir,
        tiles,
        quick=False,
        mtl=True,
        unique=True,
        isodate=mtltime,
    )
    
    # AR mtl: removing by hand BACKUP_BRIGHT for sv3/BACKUP
    # AR mtl: using an indirect way to find if program=backup,
    # AR mtl:   to avoid the need of an extra program argument
    # AR mtl:   for sv3, there is no secondary-backup, so no ambiguity
    if (survey == "sv3") & ("backup" in mtldir):
        from desitarget.sv3.sv3_targetmask import mws_mask

        keep = (d["SV3_MWS_TARGET"] & mws_mask["BACKUP_BRIGHT"]) == 0
        d = d[keep]

    #AJR added this in
    columns = [key for key in minimal_target_columns if key not in d.dtype.names]
    #tcol = ['SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET','SV3_SCND_TARGET']
    #for col in tcol:
    #    columns.append(col) 

    # Does not bail if columns already present
    if len(columns) > 0:
        print('Retrieving {}'.format(columns))

        # AR adding PLATE_RA, PLATE_DEC, PLATE_REF_EPOCH ?
        d = inflate_ledger(d, targdir, columns=columns, header=False, strictcols=False, quick=True)
        
    if add_plate_cols:
        d = Table(d)
        d["PLATE_RA"] = d["RA"]
        d["PLATE_DEC"] = d["DEC"]
        d["PLATE_REF_EPOCH"] = d["REF_EPOCH"]
        d = d.as_array()
    
    # AR mtl: PMRA, PMDEC: convert NaN to zeros
    d = force_finite_pm(d)
    # AR mtl: update RA, DEC, REF_EPOCH using proper motion?
    if pmcorr == "y":
        if pmtime_utc_str is None:
            sys.exit(1)
        d = update_nowradec(d, gaiadr, pmtime_utc_str)
    else:
        d = force_nonzero_refepoch(
            d, gaia_ref_epochs[gaiadr]
        )
    d = Table(d)

    outfndir = '/'.join(outfn.split('/')[:-1])

    print('Writing targ. to {}'.format(outfndir))
    
    if not os.path.exists(outfndir):
        os.makedirs(outfndir, exist_ok=True)

    d.write(outfn,format='fits', overwrite=True)
    del d
    return True
    # AR mtl: write fits
    #n, tmpfn = io.write_targets(tmpoutdir, d, indir=mtldir, indir2=targdir, survey=survey, subpriority=True)
    #_ = mv_write_targets_out(tmpfn, tmpoutdir, outfn)

    # AR mtl: update header if pmcorr = "y"
    if pmcorr == "y":
        fd = fitsio.FITS(outfn, "rw")
        fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
        fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
        fd.close()
