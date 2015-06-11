import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from bsfh import sps_basis
from bsfh.gp import GaussianProcess

from plotting import *

sps = sps_basis.StellarPopBasis()
#gp = GaussianProcess(None, None)
import george
kernel = (george.kernels.WhiteKernel(0.0) +
          0.0 * george.kernels.ExpSquaredKernel(0.0))
gp = george.GP(kernel, solver=george.HODLRSolver)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        resfile=sys.argv[1]
    else:
        resfile = 'results/ggc_mock.u0.t1.0_z0.0_a0.5_1426268715_mcmc'
    model_file = resfile.replace('_mcmc','_model')
    nsample = 10
    
    res, pr, mod = bread.read_pickles(resfile, model_file=model_file)
    obsdat = res['obs']
    
    fsamples = np.random.uniform(0,1,nsample)
    thetas, start, samples = theta_samples(res, samples=fsamples, start=0.75, thin=1)
    fwave, mospec, mounc, specvecs = comp_samples_fullspec(thetas, mod, obsdat, sps=sps,
                                                  gp=None)
    pwave, mosed, mosed_unc, pvecs = comp_samples_phot(thetas, mod, obsdat, sps=sps)

    tspec, tphot, ttheta = true_sed(mod, obsdat, sps=sps, fullspec=True)
    twave = sps.ssp.wavelengths
    
    #calvec = np.copy(obsdat['calibration'])
    #if np.size(calvec) == 1:
    #    calvec = np.zeros(len(fwave)) + calvec
    #else:
    #    calvec = calvec[obsdat['mask']]
    #norm = obsdat['normalization_guess'] * obsdat['rescale']

    obsdat['normalization_guess'] = 1.0
    fig = pl.figure(figsize=(6,8))
    #gs1 = gridspec.GridSpec(3, 1)
    #gs1.update(left=0.05, right=0.48, wspace=0.05)
    #oax = pl.subplot(gs1[:-1,0])
    #cax = pl.subplot(gs1[-1,0])
    
    gs2 = gridspec.GridSpec(4, 2)
    #gs2.update(left=0.55, right=0.98, hspace=0.5)
    gs2.update(left = 0.05, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])
    
    sfig, sax = sedfig(fwave, specvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax),
                       peraa=True)
    sax.plot(twave, tspec, color='black', label='True spectrum', alpha = 0.75)

    sax.set_xscale('log')
    sax.set_xlim(2.5e3, 1.6e4)
    sax.legend(loc=0, prop={'size':12})
    
    #sfig.show()
    #sfig.savefig('example_sed.png')
    
    pnames = ['mass', 'tage', 'dust2', 'zmet']
    samples, pnames_ord = hist_samples(res, mod, pnames, start=0.5, thin=10)
    truths = [obsdat['mock_params'][k] for k in pnames_ord]
    hfig, haxes = histfig(samples, pnames_ord, truths = truths,
                          color='green', fax=(None, haxes))
    haxes[0].legend(loc=0, prop={'size':8})


    fig.savefig(resfile.replace('_mcmc','.dashboard.pdf'))
    #gs1.tight_layout(fig)
    #hfig.show()
    #hfig.savefig('example_hist.png')

    #[fig.add_axes(ax) for ax in [oax, cax, sax] + haxes.flatten()]
    #fig.show()
