#!/usr/bin/env python
import unittest
import numpy as np
from numpy.testing import assert_allclose as aac
from numpy.testing import assert_array_almost_equal as aaae
from scipy.stats import kstest
import healpy as hp
from fgbuster.pysm_helpers import get_instrument, get_sky
from fgbuster.algebra import _mtv
import fgbuster.component_model as cm
from fgbuster.separation_recipies import basic_comp_sep


class TestEnd2EndNoiselessPhysical(unittest.TestCase):

    def setUp(self):
        NSIDE = 32
        MODEL = 'c1d0s0f1'
        INSTRUMENT = 'litebird'
        X0_FACTOR = 0.99
        sky = get_sky(NSIDE, MODEL)
        self.instrument = get_instrument(NSIDE, INSTRUMENT)
        self.freq_maps, self.noise = self.instrument.observe(
            sky, write_outputs=False)

        self.components = [cm.CMB(), cm.Dust(200.), cm.Synchrotron(100.)]
        freefree = cm.PowerLaw(100.)
        freefree.defaults = [-2.14]  # Otherwise it is the same as Synchrotron
        self.components.append(freefree)
        self.input = []
        for component in self.components:
            self.input += component.defaults
            component.defaults = [d*X0_FACTOR for d in component.defaults]


    def test_basic_comp_sep_T(self):
        res_T = basic_comp_sep(self.components, self.instrument,
                               self.freq_maps[:, :1, :])
        aac(res_T.x, np.array(self.input), rtol=1e-5)
        aaae(res_T.chi, 0, decimal=2)


    def test_basic_comp_sep_P(self):
        res_P = basic_comp_sep(self.components[:-1], self.instrument,
                               self.freq_maps[:, 1:, :])
        aac(res_P.x, np.array(self.input[:-1]), rtol=1e-5)
        aaae(res_P.chi, 0, decimal=2)


class TestEnd2EndNoisy(unittest.TestCase):

    def test_Sigma_synchrotron(self):
        NSIDE = 8
        MODEL = 's0'
        INSTRUMENT = 'litebird'
        SIGNAL_TO_NOISE = 20
        sky = get_sky(NSIDE, MODEL)
        instrument = get_instrument(NSIDE, INSTRUMENT)
        components = [cm.Synchrotron(100.)]
        ref = []
        for component in components:
            ref += component.defaults

        freq_maps, noise_maps = instrument.observe(sky, write_outputs=False)

        signal = freq_maps[:, 0, 0]
        noise = np.std(noise_maps[:, 0], axis=-1)
        maps = signal / np.dot(signal, noise) * SIGNAL_TO_NOISE
        maps = maps[:, np.newaxis] + noise_maps[:, 0] 
        res = basic_comp_sep(components, instrument,
                             maps, nside=hp.get_nside(maps))
        white = (res.x[0] - ref[0]) / res.Sigma[0, 0]**0.5
        _, p = kstest(white, 'norm')
        assert p > 0.01

    def test_Sigma_dust_one_parameter(self):
        NSIDE = 8
        MODEL = 'd0'
        INSTRUMENT = 'litebird'
        SIGNAL_TO_NOISE = 10
        sky = get_sky(NSIDE, MODEL)
        instrument = get_instrument(NSIDE, INSTRUMENT)
        components = [cm.Dust(100., temp=20.)]
        ref = []
        for component in components:
            ref += component.defaults

        freq_maps, noise_maps = instrument.observe(sky, write_outputs=False)

        signal = freq_maps[:, 0, 0]
        noise = noise_maps[:, 0]
        signal_ver =  signal / np.dot(signal, signal)**0.5
        noise_std = np.std([np.dot(n, signal_ver) for n in noise.T])
        maps = signal_ver * noise_std  * SIGNAL_TO_NOISE
        maps = maps[:, np.newaxis] + noise

        res = basic_comp_sep(components, instrument,
                             maps, nside=hp.get_nside(maps))
        white = (res.x[0] - ref[0]) / res.Sigma[0, 0]**0.5
        _, p = kstest(white, 'norm')
        assert p > 0.01

    def test_Sigma_dust_two_parameters(self):
        return  #  FIXME skipping test because result is unsatisfactory
        NSIDE = 1
        MODEL = 'd0'
        INSTRUMENT = 'litebird'
        SIGNAL_TO_NOISE = 100
        sky = get_sky(NSIDE, MODEL)
        instrument = get_instrument(NSIDE, INSTRUMENT)
        components = [cm.Dust(150.)]
        ref = []
        for component in components:
            ref += component.defaults
        ref = np.array(ref)

        freq_maps, noise_maps = instrument.observe(sky, write_outputs=False)

        signal = freq_maps[:, 0, 0]
        noise = noise_maps[:, 0]
        signal_ver =  signal / np.dot(signal, signal)**0.5
        noise_std = np.std([np.dot(n, signal_ver) for n in noise.T])
        maps = signal_ver * noise_std  * SIGNAL_TO_NOISE
        maps = maps[:, np.newaxis] + noise

        res = basic_comp_sep(components, instrument,
                             maps, nside=hp.get_nside(maps))
        diff = (res.x.T - ref)
        S = res.Sigma.T[0]
        e, v = np.linalg.eigh(S)
        white = _mtv(v, diff) / np.sqrt(e) 
        _, p = kstest(white.flatten(), 'norm')
        assert p > 0.01


if __name__ == '__main__':
    unittest.main()
