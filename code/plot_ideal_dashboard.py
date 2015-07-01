import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec
import matplotlib


from bsfh import read_results as bread
from bsfh import sps_basis
from bsfh.gp import ExpSquared

from plotting import *
matplotlib.colors.cnames.update(newcolors)


sps = sps_basis.StellarPopBasis()
import george
kernel = (george.kernels.WhiteKernel(0.0) +
          0.0 * george.kernels.ExpSquaredKernel(0.0))
#gp = george.GP(kernel)#, solver=george.HODLRSolver)
gp = ExpSquared(None, None)


param_name_map = {'tage':'Age (Gyr)',
                  'mass': '$M_*$ $(M_\odot)$',
                  'dust2':'$A_V$ (mag)',
                  'zmet': '$\log Z/Z_\odot$',
                  'sigma_smooth': '$\sigma_{{LSF}}$',
                  'zred': '${\it z}$',
                  }

pnmap = param_name_map

def format_sedax(sax):
    sax.set_xscale('log')
    sax.set_xlim(3e3, 1.6e4)
    ticks = list(sax.get_xlim()) + [4e3, 6e3, 10e3]
    sax.set_xticks(ticks)
    sax.set_xticklabels(['{:4.0f}'.format(t) for t in ticks], fontsize=8)
    sax.set_xlabel('$\lambda (\AA)$', fontsize=12)
    sax.set_yscale('log')
    sax.tick_params(axis='both', which='major', labelsize=8)
    #sax.set_yticklabels(sax.get_yticklabels(), fontsize=8)
    sax.set_ylabel(sax.get_ylabel(), fontsize=12)
    sax.set_ylim(3e-14, 1e-12)
    sax.legend(loc=0, prop={'size':12})
    return sax

if __name__ == "__main__":

    if len(sys.argv) > 1:
        resfile=sys.argv[1]
    else:
        resfile = 'results/ggc_mock_ideal.c0.t9.0_z0.0_a0.5_1430261146_mcmc'
    model_file = resfile.replace('_mcmc','_model')
    basecolor = 'green'
    pointcolor='orange'
    nsample = 10

    #read the results file
    res, pr, mod = bread.read_pickles(resfile, model_file=model_file)
    obsdat = res['obs']

    # Get samples and component spectra
    fsamples = np.random.uniform(0,1,nsample)
    thetas, start, samples = theta_samples(res, samples=fsamples, start=0.75, thin=1)
    fwave, _, _, fspecvecs = comp_samples_fullspec(thetas, mod, obsdat,
                                                   sps=sps, gp=None)
    mwave, mospec, mounc, specvecs = comp_samples(thetas, mod, obsdat, sps=sps,
                                                  gp=gp)
    pwave, mosed, mosed_unc, pvecs = comp_samples_phot(thetas, mod, obsdat, sps=sps)
    # Get truth spectrum and sed
    tspec, tphot, ttheta = true_sed(mod, obsdat, sps=sps, fullspec=True)
    twave = sps.ssp.wavelengths
    tspec /= obsdat['normalization_guess']
    
    norm = obsdat['normalization_guess'] * obsdat['rescale']
    #obsdat['normalization_guess'] = 1.0

    #### Build Figure ###
    fig = pl.figure(figsize=(6,8))
    # Set up right hand side
    gs2 = gridspec.GridSpec(4, 2)
    gs2.update(left=0.05, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])

    # Intrinsic SED figure
    sfig, sax = sedfig(fwave, fspecvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax),
                       peraa=True, basecolor=basecolor, pointcolor=pointcolor)
    sax.plot(twave, tspec, color='black', label='True spectrum', alpha=0.6)
    sax = format_sedax(sax)
    
    # Posterior parameter histograms
    pnames = ['mass', 'tage', 'dust2', 'zmet']
    samples, pnames_ord = hist_samples(res, mod, pnames, start=0.75, thin=5)
    truths = [obsdat['mock_params'][k] for k in pnames_ord]
    hfig, haxes = histfig(samples, pnames_ord, truths=truths, pname_map=pnmap,
                          basecolor=basecolor, fax=(None, haxes))
    haxes[0].legend(loc=0, prop={'size':8})

    # Save
    fig.suptitle('Ideal Case')
    fig.savefig(resfile.replace('_mcmc','.dashboard.pdf'))
