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


def calfig(wave, calvec, specvecs, norm=1.0, fax=None):
    """Plot the calibration and posterior samples of it
    """
    if fax is None:
        cfig, cax = pl.subplots()
    else:
        cfig, cax = fax
    #plot the calibration vector 
    cax.plot(wave, calvec, color='black', label='Mock Truth',
             linewidth=3.0)
    # and posterior samples of it
    for i, specs in enumerate(specvecs):
        if i==0:
            label = 'Posterior sample'
        else:
            label = None
        cax.plot(wave, norm * np.exp(np.log(specs[1]) + specs[2]),
                 color='green', alpha=0.3, label=label)
    
    return cfig, cax

def obsfig(wave, obsvec, specvecs, unc=None, fax=None):
    """Plot the observed spectrum and posterior samples of it
    """
    if fax is None:
        ofig, oax = pl.subplots()
    else:
        ofig, oax = fax
    # Plot posterior samples of the observed spectrum
    for i, specs in enumerate(specvecs):
        if i==0:
            label = 'Posterior samples'
        else:
            label = None
        oax.plot(wave, specs[3],
                 color='green', alpha = 0.3, label=label)
    #plot the observation itself
    if unc is not None:
        x, y, e = wave, obsvec, unc
        oax.fill_between(x, y-e, y+e, facecolor='grey', alpha=0.3)
    oax.plot(wave, obsvec, color='black', label='Mock Observed Spectrum',
             linewidth=1.0, alpha=1.0)
    return ofig, oax

def sedfig(wave, specvecs, phot, photvecs, norm = 1.0, fax=None):
    """Plot the photometric SED, posterior samples of it, and
    posterior samples of the intrinsic spectrum.
    """
    if fax is None:
        sfig, sax = pl.subplots()
    else:
        sfig, sax = fax
    pwave, sed, sed_unc = phot
    # to convert from f_lambda cgs/AA to lambda*f_lambda cgs
    sconv = wave * norm
    # to convert from maggies to nu * f_nu cgs
    pconv = 3631e-23 * 2.998e18/pwave

    for i, (specs, seds) in enumerate(zip(specvecs, photvecs)):
        if i==0:
            label = 'Posterior samples'
        else:
            label = None
        sax.plot(wave, specs[0] * sconv,
                 color='green', alpha=0.3, label=label)
        sax.plot(pwave, seds[0] * pconv, markersize=8.0, linestyle='',
                 marker='o', color='magenta', label=label)

    sax.errorbar(pwave, sed * pconv, yerr=sed_unc * pconv,
                 marker='o', markersize=8.0,
                 color='black', linestyle='', label='Mock Photometry')

    return sfig, sax

def hist_samples(res, model, showpars, start=0, thin=1, **kwargs):
    
    nw, niter = res['chain'].shape[:-1]
    parnames = np.array(model.theta_labels())
    start_index = np.floor(start * (niter-1)).astype(int)
    flatchain = res['chain'][:,start_index::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    ind_show = np.array([p in showpars for p in parnames], dtype= bool)
    flatchain = flatchain[:,ind_show]
    return flatchain, parnames[ind_show]

def histfig(samples, parnames, truths=None, fax=None, truth_color='k', **kwargs):
    npar = len(parnames)

    if fax is None:
        nx = int(np.floor(np.sqrt(npar)))
        ny = int(np.ceil(npar*1.0/nx))
        hfig, haxes = pl.subplots(nx, ny)
    else:
        hfig, haxes = fax
        
    for i, (ax, name) in enumerate(zip(haxes.flatten(), parnames)):
        ax.hist(samples[:,i], bins=kwargs.get("bins", 50),
                histtype="stepfilled",
                color=kwargs.get("color", "k"),
                alpha = 0.5,
                label = 'posterior PDF')
        if truths is not None:
            ax.axvline(truths[i], color=truth_color, label='Mock Truth')
        ax.set_xlabel(name, fontsize=6)
        ax.set_yticklabels([])
        pl.setp(ax.get_xticklabels(), fontsize=6)
    return hfig, haxes

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
    fwave, mospec, mounc, fspecvecs = comp_samples_fullspec(thetas, mod, obsdat, sps=sps,
                                                  gp=gp)

    mwave, mospec, mounc, specvecs = comp_samples(thetas, mod, obsdat, sps=sps,
                                                  gp=gp)
    pwave, mosed, mosed_unc, pvecs = comp_samples_phot(thetas, mod, obsdat, sps=sps)

    calvec = obsdat['calibration']
    if np.size(calvec) == 1:
        calvec = np.zeros(len(mwave)) + calvec
    else:
        calvec = calvec[obsdat['mask']]
    norm = obsdat['normalization_guess'] * obsdat['rescale']

    fig = pl.figure(figsize=(10,8))
    gs1 = gridspec.GridSpec(3, 1)
    gs1.update(left=0.05, right=0.48, wspace=0.05)
    oax = pl.subplot(gs1[:-1,0])
    cax = pl.subplot(gs1[-1,0])
    
    gs2 = gridspec.GridSpec(4, 2)
    gs2.update(left=0.55, right=0.98, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])

    cfig, cax = calfig(mwave, calvec, specvecs, norm=norm, fax=(None, cax))
    cax.set_ylabel(r'Calibration $F_{{obs}}/F_{{\lambda, intrinsic}}$')
    cax.legend(loc=0, prop={'size':8})
    #cfig.show()
    #cfig.savefig('example_cal.png')
    
    sfig, sax = sedfig(mwave, specvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax))
    sax.set_xscale('log')
    sax.set_xlim(2.5e3, 1.6e4)
    sax.legend(loc=0, prop={'size':12})
    sax.set_ylabel(r'$\lambda F_{{\lambda}}$ (intrinsic, cgs)')
    #sfig.show()
    #sfig.savefig('example_sed.png')
    
    ofig, oax = obsfig(mwave, mospec, specvecs, unc=mounc, fax=(None,oax))
    oax.set_ylabel(r'Observed Spectrum $F_{{obs}}$ (unknown units)')
    oax.legend(loc=0, prop={'size':8})
    #ofig.show()
    #ofig.savefig('example_obs.png')

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
