import healpy as hp
import numpy as np


def radec2pix(ra, dec, nside=32):
    theta = np.pi / 2. - np.radians(dec)
    phi   = np.radians(ra)
    
    all_pixel_indices = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)

    return all_pixel_indices

def hp_props(nside):
    npix = hp.nside2npix(nside)
    pixel_area = hp.nside2pixarea(nside,degrees=True)

    return npix, pixel_area

def targ_hpmap(targs, norm=None, nside=32):
    pix = radec2pix(ra, dec, nside=nside)
    
    #indice of filled pixels and corrosponding targets in pixel
    filled_pixel_index, filled_targets_per_pixel = np.unique(pix, return_counts=True) 

    #no. targets per pixel, initially 0 
    target_pixel_density = np.zeros(hp.nside2npix(nside))

    #update no. targets per pixel 
    target_pixel_density[filled_pixel_index] = filled_targets_per_pixel 
    target_pixel_density[target_pixel_density == 0] = np.NaN

    if norm != None:
        target_pixel_density /= norm
    
    return  target_pixel_density

def rand_inrect(ra_lower, ra_upper, dec_lower,dec_upper,nside=32):
    ra_rand = np.random.uniform(ra_lower,ra_upper,100)
    dec_rand = np.random.uniform(dec_lower,dec_upper,100)
    
    theta = np.pi / 2. - np.radians(dec_rand.data)
    phi = np.radians(ra_rand.data)

    #indices of pixels with non-zero density, unorganised list
    all_pixel_indices = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)
    
    return np.unique(all_pixel_indices)
    