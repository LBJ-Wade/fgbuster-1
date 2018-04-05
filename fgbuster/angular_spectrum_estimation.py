""" Toolbox in order to transform maps to harmonic domain
The toolbox should eventually include the possibility to estimate power spectra with
    - Anafast
    - Xpol
    - Xpure
    - NaMaster
    - PolSpice?
"""

import numpy as np
import healpy as hp
import sys

arcmin2rad = np.pi/(60.*180.)

def TEB_spectra( IQU_map, IQU_map_2=None, ell_max=0.0, estimator=None, n_iter=10, *args, **kwargs ):
    """ Get a pre-defined PySM sky

    Parameters
    ----------
    IQU_map: float, array-like shape (Npix,) or (3, Npix)
             Either an array representing a map, or a sequence of 3 arrays
             representing I, Q, U maps
    IQU_map_2: float, array-like shape (Npix,) or (3, Npix)
               Either an array representing a map, or a sequence of 3 arrays
               representing I, Q, U maps
    ell_max: int, scalar, optional
             Maximum l of the power spectrum (default: 3*nside-1)
    n_iter: int, number of iteration 
    fwhm_beam: float, FWHM of the Gaussian beam to be used to deconvolve spectra 
    fwhm_beam_2: float, FWHM of the Gaussian beam to be used to deconvolve spectra     
    estimator: string
               choice of the power spectrum estimator, among 'anafast' (default) and NaMaster 
               the user can provide the necessary *args and **kwargs for each method 
               - nothing required for Anafast
               - w (NaMaster work space), mask (the binary mask) and mask_apo (the apodized mask)
                 apotype ('C1' by default), apodization_size (30 arcmin by default), nlb (number of ell bins)
                 if not provided, we build the necessary quantities in the code

    Returns
    -------
    sky: tuple containing the 6 angular spectra (ClTT, ClEE, ClBB, ClTE, ClTB, ClEB )
    """

    if estimator=='NaMaster':

        # XXX: to be debugged and validated
        # print("power spectrum estimation with NaMaster not ready yet", file=sys.stderr)
        print >> sys.stderr, "power spectrum estimation with NaMaster not ready yet"
        #--------------- 
        # import pymaster as nmt
        # if not w: 
        #     print('building mode coupling matrix as it is not provided')
        #     b=nmt.NmtBin( hp.npix2nside(IQU_map.shape[1]), nlb=nlb )
        #     w=nmt.NmtWorkspace()
        #     w.compute_coupling_matrix( Q, U, b )
        # if not mask_apo:
        #     if not mask:
        #         mask = np.ones(IQU_map[0].shape[0])
        #         mask[np.where(IQU_map[0] == 0)[0]] = 0.
        #     mask_apo = nmt.mask_apodization(mask, apodization_size, apotype=apotype)
        # f2 = nmt.NmtField(mask_apo, [Q,U], purify_e=True, purify_b=True)
        # field_0 = f2
        # if ((len(Q2)>1) and (len(U2)>1)):
        #     f2_2 = nmt.NmtField(mask_apo,[Q2,U2], purify_e=True, purify_b=True)
        #     field_1 = f2_2
        # else: field_1 = f2
        # cl_coupled=nmt.compute_coupled_cell( field_0, field_1 )
        # cl_decoupled=w.decouple_cell( cl_coupled )
        # return cl_decoupled

    elif (estimator=='Xpol' or estimator=='Xpure'):

        # XXX: to be written, using something close to what Julien Peloton
        #      did for S4CMB
        # print("power spectrum estimation with Xpol/Xpure not ready yet", file=sys.stderr)
        print >> sys.stderr, "power spectrum estimation with Xpol/Xpure not ready yet"
        #--------------- 
        # write_maps_a_la_xpure(sky_out_tot, name_out=name_out,
        #                           output_path='xpure/maps')
        # write_weights_a_la_xpure(sky_out_tot, name_out=name_out,
        #                              output_path='xpure/masks')
        
        # params_xpure = import_string_as_module(args.inifile_xpure)

        # batch_file = 'sim{:03d}_{}_{}_{}.batch'.format(
        #             args.sim_number,
        #             params.tag,
        #             params.name_instrument,
        #             params.name_strategy)
        
        # create_batch(batch_file, name_out, params, params_xpure)

        # qsub = commands.getoutput('sbatch ' + batch_file)
        # print(qsub)

    else:

        #if isinstance(IQU_map, tuple):
        #    print >> sys.stderr, "anafast requires input map to be an array, not a tuple"

        #if IQU_map.ndim > 1:
        #    nside_input_map = hp.npix2nside(IQU_map[0].shape[0])
        #else:
        #    nside_input_map = hp.npix2nside(len(IQU_map))

        #if ell_max <= 3*nside_input_map-1:
        # 	ell_max = 3*nside_input_map-1
        pol = True if IQU_map.ndim>1 else False
        
        Cl = hp.sphtfunc.anafast(map1=IQU_map, map2=IQU_map_2, iter=n_iter, lmax=ell_max)
        
        if "fwhm_beam" in kwargs.keys():
            Bl = hp.gauss_beam(kwargs["fwhm_beam"]*arcmin2rad, lmax=ell_max, pol=pol)
        else:
            Bl = np.ones((ell_max+1,3))
        if "fwhm_beam_2" in kwargs.keys():
            Bl2 = hp.gauss_beam(kwargs["fwhm_beam_2"]*arcmin2rad, lmax=ell_max, pol=pol)
        else:
            Bl2 = Bl
        inv_Bl = 1/(Bl*Bl2);  max_inv = 1e4 # maximum inverse value
        inv_Bl[inv_Bl>max_inv] = max_inv
        if pol: #TT, EE, BB, TE, EB, TB
            Cl = [Cl[0]*inv_Bl[:,0], Cl[1]*inv_Bl[:,1], Cl[2]*inv_Bl[:,2],
                  Cl[3]*inv_Bl[:,3], Cl[4]*inv_Bl[:,1], Cl[5]*inv_Bl[:,3]]
        else:
            Cl = Cl*inv_Bl
            
            
        return Cl
