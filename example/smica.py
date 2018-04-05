from fgbuster.pysm_helpers import get_instrument, get_sky
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from fgbuster.angular_spectrum_estimation import TEB_spectra

# simulate sky
nside = 32
sky = get_sky(nside, 'c1d0s0')
instr = get_instrument(nside, 'litebird')
freq_maps = instr.observe(sky, write_outputs=False)[0]

# build an isolat binary mask
nside_hr = 1024
mask = np.zeros((12*nside_hr**2))
theta,phi = hp.pix2ang(nside, np.arange(len(mask)))
lat_cut = 15
mask[theta<np.pi/2-lat_cut*np.pi/180] = 1
mask[theta>np.pi/2+lat_cut*np.pi/180] = 1
hp.mollview(mask)

# apodize it
amask = hp.alm2map(hp.almxfl(hp.map2alm(mask, lmax=2*nside, iter=0), hp.gauss_beam(10*np.pi/180.)), nside)
amask[amask<0.01] = 0
amask[amask>0.99] = 1
hp.mollview(amask)

# compute spectral covariance matrices
nmaps = freq_maps.shape[0]
lmax = nside*2
hRTT = np.zeros((nmaps, nmaps, lmax+1))
for i in range(nmaps):
    for j in range(i, nmaps):
        hRTT[i,j,:] = TEB_spectra(freq_maps[i,:,:], IQU_map_2=freq_maps[j,:,:],
                                  ell_max=lmax, n_iter=0,
                                  mask=mask,
                                  fwhm_beam=instr.Beams[i], fwhm_beam_2=instr.Beams[j])[0]
        hRTT[j,i,:] = hRTT[i,j,:]

plt.figure()
for i in range(nmaps):
    for j in range(i, nmaps):
        plt.semilogy(hRTT[i,j,:])
plt.show()

# binning

