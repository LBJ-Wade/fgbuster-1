import healpy as hp
import numpy as np


def healpix_multiresolution(map_in, condition, power=0,
                            nside_min=0, nside_max=None):
    '''
    Parameters
    ----------
    map_in: array
        Initial map to be used
    condition: function
        Takes the map (at any resolution) and returns a boolean. True means
        that a new index should be assigned to the pixel
    power: float
        The downgraded sup-pixel is the sum over the subpixel divided the number
        of subpixels to this power
    nside_min: int
        All pixels that are still unassigned at this nside are assigned to a
        single region
    nside_max: int
        Do not trigger ne regions above this nside

    Returns
    -------
    idx:
        Map of the index of the region to which each pixel in map_in belongs to.
    nsides:
        Map of the index of the region to which each pixel in map_in belongs to.
    '''
    map_in = map_in.copy()
    mask_bad = hp.ma(map_in).mask
    idx = np.zeros_like(map_in, dtype=int)
    nside = hp.get_nside(map_in)
    nsides = [0]
    #healpix_ids = [0]
    while True:
        if nside_max is None or nside <= nside_max:
            trigger_new_pix = condition(map_in)
            map_in[trigger_new_pix] = hp.UNSEEN
            n_pix_triggered = _update_idx(idx, trigger_new_pix, mask_bad)
            #healpix_id = _update_idx(idx, trigger_new_pix, mask_bad)
            #healpix_ids += list(healpix_id)
            #nsides += [nside] * len(healpix_id)
            nsides += [nside] * n_pix_triggered
        nside //= 2
        if nside < max(nside_min, 1):
            break
        map_in = downgrade(map_in, nside, power)

    assert len(set(idx)) == idx.max()+1, '%i %i' % (len(set(idx)), idx.max()+1)
    return idx, np.array(nsides)#, np.array(healpix_ids)


def _update_idx(map_idx, trigger_new_pix, mask_bad=None):
    if hp.get_nside(trigger_new_pix) == 8:
        trigger_new_pix[341] = True
        trigger_new_pix[180] = True
    '''
    if hp.get_nside(trigger_new_pix) == 16:
        trigger_new_pix[1330] = True
    '''


    # NOTE index equal zero -> pixel still unassigned
    # Checks and preliminaries
    assert len(trigger_new_pix) <= len(map_idx)
    if hp.pixelfunc.is_ma(trigger_new_pix):
        trigger_new_pix = trigger_new_pix.data * hp.mask_good(trigger_new_pix)

    # NOTE trigger_new_pix can be true for pixels already assigned to an index
    # Build a mask of the pixels that are already assigned
    mask = np.zeros_like(map_idx, dtype=float)
    mask[map_idx != 0] = hp.UNSEEN
    mask = hp.ud_grade(mask, hp.get_nside(trigger_new_pix))
    # Intersect with trigger_new_pix and make mask bool, non-ma ndarray
    mask = hp.mask_good(mask) * trigger_new_pix
    if mask_bad is not None:
        mask *= hp.ud_grade(mask_bad.astype(float),
                            hp.get_nside(trigger_new_pix)) != 1.
    # Now mask is true if and only if the pixel needs a new index assigned to it
    idx = np.zeros_like(mask, dtype=float)
    idx[mask] = np.arange(mask.sum()) + (1 + map_idx.max())
    '''
    print('Filling holes')
    for i in np.arange(idx.size)[~mask]:
        neighbours = hp.get_all_neighbours(hp.get_nside(trigger_new_pix), i)
        if np.all(mask[neighbours]):
            idx[i] = idx.max() + 1
    print('Holes filled')
    '''

    idx = hp.ud_grade(idx, hp.get_nside(map_idx)).astype(int)

    if mask_bad is None:
        idx[map_idx != 0] = 0  # update only unassigned pixels
    else:
        idx[(map_idx != 0) | mask_bad] = 0  # update only unassigned pixels



    map_idx += idx
    return mask.sum()


def downgrade(map_in, nside, power=0):
    nside_in = hp.get_nside(map_in)
    assert nside_in >= nside
    assert nside > 0
    idx = hp.ud_grade(np.arange(hp.nside2npix(nside)), nside_in).astype(int)
    map_in = hp.ma(map_in)
    map_out = np.bincount(idx[~map_in.mask], map_in[~map_in.mask],
                          minlength=hp.nside2npix(nside))
    nhits = np.bincount(idx[~map_in.mask], minlength=hp.nside2npix(nside))
    map_out[nhits == 0] = hp.UNSEEN
    map_out = hp.ma(map_out)
    if power:
        return map_out / nhits**power
    return map_out

'''
import pylab as pl
m = np.arange(12*8**2, dtype=float)
m[568] = hp.UNSEEN
condition = lambda x: x > m[-1] / 1.001
#
m = hp.read_map('../../../analyses/map.fits')
condition = lambda x: x**2 > 1000
idx = healpix_multiresolution(hp.ma(m), condition, power=0.5)
hp.mollview(idx)
pl.show()
'''
