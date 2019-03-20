# FGBuster
# Copyright (C) 2019 Davide Poletti, Josquin Errard and the FGBuster developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

""" All the routines for making all Josquin's lovely plots
"""
from corner import corner
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import healpy as hp


def corner_norm(mean, cov, *args, **kwargs):
    ''' Corner plot for multivariate gaussian

    Just like corner.corner, but you privide mean and covariance instead of `xs`
    '''
    xs = np.random.multivariate_normal(mean, cov, 100000)  # TODO: not hardcoded
    corner(xs, *args, **kwargs)


def plot_component(component, nu_min, nu_max):
    nus = np.logspace(np.log10(nu_min), np.log10(nu_max), 1000)
    emission = component.eval(nus, *(component.defaults))
    plt.loglog(nus, emission, label=type(component).__name__)


def check_fit(**kwargs):
    #chi = hp.ma(kwargs['chi'])
    chi = kwargs['chi']
    try:
        output = kwargs['output']
        del kwargs['output']
    except KeyError:
        output = False

    mask = np.any(hp.ma(chi).mask, axis=tuple(range(chi.ndim-1)))

    if 'mask' in kwargs:
        mask |= kwargs['mask']

    chi = chi[..., ~mask]

    if 'patch_ids' in kwargs:
        ids = kwargs['patch_ids'][~mask]
        del kwargs['patch_ids']

        max_id = ids.max()
        digit = int(np.log10(max_id)) + 1

        for i_patch in range(max_id):
            kwargs['chi'] = chi[..., ids == i_patch]

            if output:
                tag = ('_%0'+str(digit)+'d') % i_patch
                i_dot = output.rfind('.')
                kwargs['output'] = output[:i_dot] + tag + output[i_dot:]
                print(kwargs['output'])
            else:
                assert max_id < 10, 'Too many plot to visualize, use "output"'

            check_fit(**kwargs)

        return

    assert chi.ndim > 1

    if not chi.size:
        return

    n_freq = chi.shape[0]
    truths = [0.] * n_freq
    if 'frequencies' in kwargs:
        labels = ['%i GHz' % f for f in kwargs['frequencies']]
    else:
        labels = None

    fig = corner(chi.reshape(n_freq, -1).T, truths=truths, color='k',
                 hist_kwargs=dict(density=True), labels=labels)

    x = np.linspace(-3, 3, 100)
    axes = np.array(fig.axes).reshape((n_freq, n_freq))
    for i_freq in range(n_freq):
        axes[i_freq, i_freq].plot(x, stats.norm.pdf(x))

    if n_freq > 2:
        plt.sca(axes[0,2])
        plt.cla()
        chi2 = np.sum(chi**2, axis=0).flatten()
        n, _, _ = plt.hist(chi2, fill=False, histtype='step', density=True, color='k')
        n_min = min([i for i in n if i > 0])
        plt.yscale('log')
        plt.ylim(n_min, None)
        plt.title('$\\chi^2$')
        if 'components' in kwargs:
            xmin, xmax = plt.gca().get_xlim()
            x = np.linspace(xmin, xmax, 100)
            dof = n_freq - len(kwargs['components'])
            plt.plot(x, stats.chi2.pdf(x, dof))

    if output:
        plt.savefig(output)
        plt.close()
