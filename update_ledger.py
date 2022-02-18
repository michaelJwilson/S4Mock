import os
import healpy as hp
import numpy as np

from   desitarget import io
from   desitarget.mtl import _get_mtl_nside, get_mtl_ledger_format, mtlformatdict
from   desitarget.geomask import match, match_to
from   astropy.io import ascii
from   make_mtl import make_mtl

def update_ledger(hpdirname, zcat, targets=None, obscon="DARK",
                  numobs_from_ledger=False):
    """
    Update relevant HEALPixel-split ledger files for some targets.
    Parameters
    ----------
    hpdirname : :class:`str`
        Full path to a directory containing an MTL ledger that has been
        partitioned by HEALPixel (i.e. as made by `make_ledger`).
    zcat : :class:`~astropy.table.Table`, optional
        Redshift catalog table with columns ``TARGETID``, ``NUMOBS``,
        ``Z``, ``ZWARN``, ``ZTILEID``, and ``msaddcols`` at the top of
        the code for the Main Survey.
    targets : :class:`~numpy.array` or `~astropy.table.Table`, optional
        A numpy rec array or astropy Table with at least the columns
        ``RA``, ``DEC``, ``TARGETID``, ``DESI_TARGET``, ``NUMOBS_INIT``,
        and ``PRIORITY_INIT``. If ``None``, then assume the `zcat`
        includes ``RA`` and ``DEC`` and look up `targets` in the ledger.
    obscon : :class:`str`, optional, defaults to "DARK"
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set using "obsconditions". Basically a
        check on whether the files in `hpdirname` are as expected.
    numobs_from_ledger : :class:`bool`, optional, defaults to ``True``
        If ``True`` then inherit the number of observations so far from
        the ledger rather than expecting it to have a reasonable value
        in the `zcat.`
    Returns
    -------
    Nothing, but relevant ledger files are updated.
    """
    # ADM find the general format for the ledger files in `hpdirname`.
    # ADM also returning the obsconditions.

    print('Loading: {}'.format(hpdirname))
    
    fileform, oc = io.find_mtl_file_format_from_header(hpdirname, returnoc=True)
    # ADM this is the format for any associated override ledgers.
    overrideff = io.find_mtl_file_format_from_header(hpdirname,
                                                     forceoverride=True)

    # ADM check the obscondition is as expected.
    if obscon != oc:
        msg = "File is type {} but requested behavior is {}".format(oc, obscon)
        log.critical(msg)
        raise RuntimeError(msg)

    # ADM if targets wasn't sent, that means the zcat includes
    # ADM coordinates and we can read relevant targets from the ledger.
    if targets is None:
        nside = _get_mtl_nside()
        theta, phi = np.radians(90-zcat["DEC"]), np.radians(zcat["RA"])

        # print(np.sort(theta))
        # print(np.sort(phi))
        
        pixnum = hp.ang2pix(nside, theta, phi, nest=True)
        pixnum = list(set(pixnum))
        # ADM we'll read in too many targets, here, but that's OK as
        # ADM make_mtl(trimtozcat=True) only returns the updated targets.
        targets = io.read_mtl_in_hp(hpdirname, nside, pixnum, unique=True)

    # ADM if requested, use the previous values in the ledger to set
    # ADM NUMOBS in the zcat.
    if numobs_from_ledger:
        # ADM match the zcat to the targets.
        tii, zii = match(targets["TARGETID"], zcat["TARGETID"])
        # ADM update NUMOBS in the zcat for matches.
        zcat["NUMOBS"][zii] = targets["NUMOBS"][tii] + 1
        
    # ADM run MTL, only returning the targets that are updated.
    # https://github.com/desihub/desitarget/blob/b3c58c89bbc5e07902154a0f0d890f62d4e29539/py/desitarget/mtl.py#L410
    
    mtl = make_mtl(targets, oc, zcat=zcat, trimtozcat=True, trimcols=True)

    print(mtl.dtype.names)
    
    # ADM this is redundant if targets wasn't sent, but it's quick.
    nside = _get_mtl_nside()
    theta, phi = np.radians(90-mtl["DEC"]), np.radians(mtl["RA"])
    pixnum = hp.ang2pix(nside, theta, phi, nest=True)

    # ADM loop through the pixels and update the ledger, depending
    # ADM on whether we're working with .fits or .ecsv files.
    ender = get_mtl_ledger_format()
    
    for pix in set(pixnum):
        # ADM grab the targets in the pixel.
        ii = pixnum == pix
        mtlpix = mtl[ii]

        # ADM sorting on TARGETID is neater (although not strictly
        # ADM necessary when using io.read_mtl_ledger(unique=True)).
        mtlpix = mtlpix[np.argsort(mtlpix["TARGETID"])]

        # ADM the correct filenames for this pixel number.
        fn = fileform.format(pix)
        overfn = overrideff.format(pix)

        print()
        print(fn)
        print(overfn)
        
        # ADM if an override ledger exists, update it and recover its
        # ADM relevant MTL entries.
        if os.path.exists(overfn):
            # overmtl = process_overrides(overfn)
            # ADM add any override entries TO THE END OF THE LEDGER.
            # mtlpix = vstack([mtlpix, overmtl])
            pass

        # ADM if we're working with .ecsv, simply append to the ledger.
        if ender == 'ecsv':
            f = open(fn, "a")
            ascii.write(mtlpix, f, format='no_header', formats=mtlformatdict)
            f.close()
        
    return
