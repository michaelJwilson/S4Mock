import os
import numpy as np
import healpy as hp
import numpy.lib.recfunctions as rfn
import sys
from astropy.table import Table, hstack, vstack
from astropy.io import ascii
import fitsio
from time import time, sleep
from datetime import datetime, timezone, date, timedelta
from glob import glob, iglob

from desitarget.targetmask import obsmask, obsconditions, zwarn_mask
from desitarget.targets import calc_priority, calc_numobs_more
from desitarget.targets import main_cmx_or_sv, switch_main_cmx_or_sv
from desitarget.targets import set_obsconditions, decode_targetid
from desitarget.geomask import match, match_to
from desitarget.internal import sharedmem
from desitarget import io
from desimodel.footprint import is_point_in_desi, tiles2pix
from desitarget.mtl import survey_data_model, mtldatamodel, get_utc_date
from desitarget import __version__ as dt_version

# ADM set up the DESI default logger.
from desiutil.log import get_logger

log = get_logger()

# HACK: zcat.masked throws error. 
def make_mtl(targets, obscon, zcat=None, scnd=None,
             trim=False, trimcols=False, trimtozcat=False):
    """
    Add zcat columns to a targets table, update priorities and NUMOBS.
    
    Parameters
    ----------
    targets : :class:`~numpy.array` or `~astropy.table.Table`
        A numpy rec array or astropy Table with at least the columns
        `TARGETID`, `DESI_TARGET`, `BGS_TARGET`, `MWS_TARGET` (or the
        corresponding columns for SV or commissioning) `NUMOBS_INIT` and
        `PRIORITY_INIT`. `targets` must also contain `PRIORITY` if `zcat`
        is not ``None`` (i.e. if this isn't the first time through MTL
        and/or if `targets` is itself an mtl array). `PRIORITY` is needed
        to "lock in" the state of Ly-Alpha QSOs. `targets` may also
        contain `SCND_TARGET` (or the corresponding columns for SV) if
        secondary targets are under consideration.
    obscon : :class:`str`
        A combination of strings that are in the desitarget bitmask yaml
        file (specifically in `desitarget.targetmask.obsconditions`), e.g.
        "DARK|GRAY". Governs the behavior of how priorities are set based
        on "obsconditions" in the desitarget bitmask yaml file.
    zcat : :class:`~astropy.table.Table`, optional
        Redshift catalog table with columns ``TARGETID``, ``NUMOBS``,
        ``Z``, ``ZWARN``, ``ZTILEID``, and possibly the extra columns in
        ``msaddcols`` at the top of the module.
    scnd : :class:`~numpy.array`, `~astropy.table.Table`, optional
        TYPICALLY, we have a separate secondary targets (they have their
        own "ledger"). So passing associated secondaries is DEPRECATED
        (but still works). `scnd` is kept for backwards compatibility.
        A set of secondary targets associated with the `targets`. As with
        the `target` must include at least ``TARGETID``, ``NUMOBS_INIT``,
        ``PRIORITY_INIT`` or the corresponding SV columns.
        The secondary targets will be padded to have the same columns
        as the targets, and concatenated with them.
    trim : :class:`bool`, optional
        If ``True`` (default), don't include targets that don't need
        any more observations.  If ``False``, include every input target.
    trimcols : :class:`bool`, optional, defaults to ``False``
        Only pass through columns in `targets` that are actually needed
        for fiberassign (see `desitarget.mtl.mtldatamodel`).
    trimtozcat : :class:`bool`, optional, defaults to ``False``
        Only return targets that have been UPDATED (i.e. the targets with
        a match in `zcat`). Returns all targets if `zcat` is ``None``.
        See important Notes about `trimtozcat`=``False``!!!
    Returns
    -------
    :class:`~astropy.table.Table`
        MTL Table with targets columns plus:
        * NUMOBS_MORE    - number of additional observations requested
        * PRIORITY       - target priority (larger number = higher priority)
        * TARGET_STATE   - the observing state that corresponds to PRIORITY
        * OBSCONDITIONS  - replaces old GRAYLAYER
        * TIMESTAMP      - time that (this) make_mtl() function was run
        * VERSION        - version of desitarget used to run make_mtl()
    Notes
    -----
    - Sources in the zcat with `ZWARN` of `BAD_SPECQA|BAD_PETALQA|NODATA`
      are always ignored.
    - As sources with `BAD_SPECQA|BAD_PETALQA|NODATA` are ignored they
      will have ridiculous z-entries if `trimtozcat`=``False`` is passed!
    - The input `zcat` MAY BE MODIFIED. If a desideratum is that `zcat`
      remains unaltered, make sure to copy `zcat` before passing it.
    """
    start = time()
    # ADM set up the default logger.
    from desiutil.log import get_logger
    log = get_logger()

    # ADM if trimcols was passed, reduce input target columns to minimal.
    if trimcols:
        mtldm = switch_main_cmx_or_sv(mtldatamodel, targets)
        # ADM the data model for mtl depends on the survey type.
        _, _, survey = main_cmx_or_sv(mtldm)
        mtldm = survey_data_model(mtldm, survey=survey)
        cullcols = list(set(targets.dtype.names) - set(mtldm.dtype.names))
        if isinstance(targets, Table):
            targets.remove_columns(cullcols)
        else:
            targets = rfn.drop_fields(targets, cullcols)

    # ADM determine whether the input targets are main survey, cmx or SV.
    colnames, masks, survey = main_cmx_or_sv(targets, scnd=True)
    # ADM set the first column to be the "desitarget" column
    desi_target, desi_mask = colnames[0], masks[0]
    scnd_target, scnd_mask = colnames[-1], masks[-1]

    # ADM if secondaries were passed, concatenate them with the targets.
    if scnd is not None:
        nrows = len(scnd)
        log.info('Pad {} primary targets with {} secondaries...t={:.1f}s'.format(
            len(targets), nrows, time()-start))
        padit = np.zeros(nrows, dtype=targets.dtype)
        sharedcols = set(targets.dtype.names).intersection(set(scnd.dtype.names))
        for col in sharedcols:
            padit[col] = scnd[col]
        targets = np.concatenate([targets, padit])
        # APC Propagate a flag on which targets came from scnd
        is_scnd = np.repeat(False, len(targets))
        is_scnd[-nrows:] = True
        log.info('Done with padding...t={:.1f}s'.format(time()-start))

    # Trim targets from zcat that aren't in original targets table.
    # ADM or that didn't actually obtain an observation.
    if zcat is not None:
        ok = np.in1d(zcat['TARGETID'], targets['TARGETID'])
        num_extra = np.count_nonzero(~ok)
        if num_extra > 0:
            log.info("Ignoring {} z entries that aren't in the input target list"
                     " (e.g. likely skies, secondaries-when-running-primary, "
                     "primaries-when-running-secondary, etc.)".format(num_extra))
            zcat = zcat[ok]
        # ADM also ignore anything with NODATA set in ZWARN.
        nodata = zcat["ZWARN"] & zwarn_mask["NODATA"] != 0
        num_nod = np.sum(nodata)
        if num_nod > 0:
            log.info("Ignoring a further {} zcat entries with NODATA set".format(
                num_nod))
            zcat = zcat[~nodata]
        # SB ignore targets that failed QA: ZWARN bits BAD_SPECQA|BAD_PETALQA
        badqa = zcat["ZWARN"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
        num_badqa = np.sum(badqa)
        if num_badqa > 0:
            log.info(f"Ignoring a further {num_badqa} zcat entries with BAD_SPECQA or BAD_PETALQA set")
            zcat = zcat[~badqa]
        # ADM simulations (I think) and some unit tests expect zcat to
        # ADM be modified by make_mtl().
        if num_extra > 0 or num_nod > 0 or num_badqa > 0:
            msg = "The size of the zcat has changed, so it won't be modified!"
            log.warning(msg)

    n = len(targets)
    # ADM if a redshift catalog was passed, order it to match the input targets
    # ADM catalog on 'TARGETID'.
    if zcat is not None:
        # ADM find where zcat matches target array.
        zmatcher = match_to(targets["TARGETID"], zcat["TARGETID"])
        ztargets = zcat

        masked   = hasattr(ztargets, 'masked')
        
        if ztargets.masked:
            unobs = ztargets['NUMOBS'].mask
            ztargets['NUMOBS'][unobs] = 0
            unobsz = ztargets['Z'].mask
            ztargets['Z'][unobsz] = -1
            unobszw = ztargets['ZWARN'].mask
            ztargets['ZWARN'][unobszw] = -1
    else:
        ztargets = Table()
        ztargets['TARGETID'] = targets['TARGETID']
        ztargets['NUMOBS'] = np.zeros(n, dtype=np.int32)
        ztargets['Z'] = -1 * np.ones(n, dtype=np.float32)
        ztargets['ZWARN'] = -1 * np.ones(n, dtype=np.int32)
        # ADM a catch all for added zcat columns.
        xtradescr = [dt for dt in mtldatamodel.dtype.descr if dt[0] == "ZTILEID"]
        if survey == 'main':
            xtradescr += msaddcols.dtype.descr
        for xctuple in xtradescr:
            xtracol, dt = xctuple
            # ADM default to "_" instead of -1 for strings.
            if isinstance(np.empty(1, dtype=dt).item(), str):
                ztargets[xtracol] = np.full(n, "-", dtype=dt)
            else:
                ztargets[xtracol] = np.full(n, -1, dtype=dt)
        # ADM if zcat wasn't passed, there is a one-to-one correspondence
        # ADM between the targets and the zcat.
        zmatcher = np.arange(n)
        
    # ADM extract just the targets that match the input zcat.
    targets_zmatcher = targets[zmatcher]

    # ADM special cases for SV3.
    if survey == "sv3":
        if zcat is not None:
            # ADM a necessary hack as we created ledgers for SV3 with
            # ADM NUMOBS_INIT==9 then later decided on NUMOBS_INIT==3.
            ii = targets_zmatcher["NUMOBS_INIT"] == 9
            targets_zmatcher["NUMOBS_INIT"][ii] = 3
            # ADM make sure to also force a permanent change of state for
            # ADM the actual *targets* that will be returned as the mtl.
            targets["NUMOBS_INIT"][zmatcher[ii]] = 3
        if (obsconditions.mask(obscon) & obsconditions.mask("DARK")) != 0:
            # ADM In dark time, if a QSO target is above feasible galaxy
            # ADM redshifts, NUMOBS should behave like a QSO, not an ELG.
            ii = targets_zmatcher[desi_target] & desi_mask["QSO"] != 0
            # ADM the secondary bit-names that correspond to primary QSOs.
            sns = [bn for bn in scnd_mask.names() if scnd_mask[bn].flavor == 'QSO']
            for sn in sns:
                ii |= targets_zmatcher[scnd_target] & scnd_mask[sn] != 0
            # ADM above feasible galaxy redshifts (with no warning).
            ii &= ztargets['Z'] > 1.6
            ii &= ztargets['ZWARN'] == 0
            targets_zmatcher["NUMOBS_INIT"][ii] = desi_mask["QSO"].numobs

    # ADM update the number of observations for the targets.
    ztargets['NUMOBS_MORE'] = calc_numobs_more(targets_zmatcher, ztargets, obscon)

    # ADM assign priorities. Only things in the zcat can have changed
    # ADM priorities. Anything else is assigned PRIORITY_INIT, below.
    priority, target_state = calc_priority(
        targets_zmatcher, ztargets, obscon, state=True)

    # If priority went to 0==DONOTOBSERVE or 1==OBS or 2==DONE, then
    # NUMOBS_MORE should also be 0.
    # ## mtl['NUMOBS_MORE'] = ztargets['NUMOBS_MORE']
    ii = (priority <= 2)
    log.info('{:d} of {:d} targets have priority <=2, setting N_obs=0.'.format(
        np.sum(ii), n))
    ztargets['NUMOBS_MORE'][ii] = 0

    # - Set the OBSCONDITIONS mask for each target bit.
    obsconmask = set_obsconditions(targets)

    # APC obsconmask will now be incorrect for secondary-only targets. Fix this
    # APC using the mask on secondary targets.
    if scnd is not None:
        obsconmask[is_scnd] = set_obsconditions(targets[is_scnd], scnd=True)

    # ADM set up the output mtl table.
    mtl = Table(targets)
    mtl.meta['EXTNAME'] = 'MTL'

    # ADM use the Main Survey data model, if appropriate.
    mtldm = survey_data_model(mtldatamodel, survey=survey)

    # ADM add a placeholder for the secondary bit-mask, if it isn't there.
    if scnd_target not in mtl.dtype.names:
        mtl[scnd_target] = np.zeros(len(mtl),
                                    dtype=mtldm["SCND_TARGET"].dtype)

    # ADM initialize columns to avoid zero-length/missing/format errors.
    zcols = ["NUMOBS_MORE", "NUMOBS", "Z", "ZWARN", "ZTILEID"]
    if survey == 'main':
        zcols += list(msaddcols.dtype.names)
    for col in zcols + ["TARGET_STATE", "TIMESTAMP", "VERSION"]:
        mtl[col] = np.empty(len(mtl), dtype=mtldm[col].dtype)

    # ADM any target that wasn't matched to the ZCAT should retain its
    # ADM original (INIT) value of PRIORITY and NUMOBS.
    mtl['NUMOBS_MORE'] = mtl['NUMOBS_INIT']
    mtl['PRIORITY'] = mtl['PRIORITY_INIT']
    mtl['TARGET_STATE'] = "UNOBS"
    # ADM add the time and version of the desitarget code that was run.
    mtl["TIMESTAMP"] = get_utc_date(survey=survey)
    mtl["VERSION"] = dt_version

    # ADM now populate the new mtl columns with the updated information.
    mtl['OBSCONDITIONS'] = obsconmask
    mtl['PRIORITY'][zmatcher] = priority
    mtl['TARGET_STATE'][zmatcher] = target_state
    for col in zcols:
        mtl[col][zmatcher] = ztargets[col]
    # ADM add ZTILEID, other columns, if passed, otherwise we're likely
    # ADM to be working with non-ledger-based mocks and can let it slide.
    xtradescr = [dt for dt in mtldatamodel.dtype.descr if dt[0] == "ZTILEID"]
    if survey == 'main':
        xtradescr += msaddcols.dtype.descr
    for xctuple in xtradescr:
        xtracol, dt = xctuple
        if xtracol in ztargets.dtype.names:
            mtl[xtracol][zmatcher] = ztargets[xtracol]
        else:
            # ADM default to "_" instead of -1 for strings.
            if isinstance(np.empty(1, dtype=dt).item(), str):
                mtl[xtracol] = "-"
            else:
                mtl[xtracol] = -1

    # Filter out any targets marked as done.
    if trim:
        notdone = mtl['NUMOBS_MORE'] > 0
        log.info('{:d} of {:d} targets are done, trimming these'.format(
            len(mtl) - np.sum(notdone), len(mtl))
        )
        mtl = mtl[notdone]

    # Filtering can reset the fill_value, which is just wrong wrong wrong
    # See https://github.com/astropy/astropy/issues/4707
    # and https://github.com/astropy/astropy/issues/4708
    mtl['NUMOBS_MORE'].fill_value = -1

    log.info('Done...t={:.1f}s'.format(time()-start))

    if trimtozcat:
        return mtl[zmatcher]
    return mtl
