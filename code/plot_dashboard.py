import sys, pickle, matplotlib
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from bsfh import sps_basis
from bsfh.gp import ExpSquared

from plotting import *
matplotlib.colors.cnames.update(newcolors)

sps = sps_basis.StellarPopBasis()
gp = ExpSquared(None, None)

param_name_map = {'tage':'Age (Gyr)',
                  'mass': '$M_*$ $(M_\odot/10^{5})$',
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
    sax.set_ylabel(sax.get_ylabel(), fontsize=12)
    sax.set_ylim(3e-14, 1e-12)
    sax.legend(loc=0, prop={'size':12})
    return sax

def format_calax(cax, norm):
    cax.set_ylabel(r'Calibration $F_{{obs}}/F_{{\lambda, intrinsic}}$',
                   fontsize=12)
    cax.legend(loc=0, prop={'size':8})
    cax.set_ylim(0.2*norm, 1.6*norm)

    #cax.set_xlim(3e3, 1.6e4)
    #ticks = list(cax.get_xlim()) + [4e3, 6e3, 10e3]
    #cax.set_xticks(ticks)
    #cax.set_xticklabels(['{:4.0f}'.format(t) for t in ticks], fontsize=8)
    cax.set_xlabel('$\lambda (\AA)$', fontsize=12)
    cax.tick_params(axis='both', which='major', labelsize=8)
    cax.set_ylabel(cax.get_ylabel(), fontsize=12)
    cax.legend(loc=0, prop={'size':12})
    return cax

if __name__ == "__main__":

    bcolor='magenta'
    if len(sys.argv) > 1:
        resfile=sys.argv[1]
        try:
            bcolor = sys.argv[2]
        except:
            pass
    else:
        resfile = 'results/ggc_mock.u0.t1.0_z0.0_a0.5_1426268715_mcmc'
    model_file = resfile.replace('_mcmc','_model')
    pcolor='orange'
    nsample = 10
    
    res, pr, mod = bread.read_pickles(resfile, model_file=model_file)
    obsdat = res['obs']
    
    # Get samples and component spectra
    fsamples = np.random.uniform(0,1,nsample)
    thetas, start, samples = theta_samples(res, samples=fsamples, start=0.95, thin=1)
    fwave, _, _, fspecvecs = comp_samples_fullspec(thetas, mod, obsdat,
                                                   sps=sps, gp=None)
    mwave, mospec, mounc, specvecs = comp_samples(thetas, mod, obsdat,
                                                  sps=sps, gp=gp)
    pwave, mosed, mosed_unc, pvecs = comp_samples_phot(thetas, mod, obsdat, sps=sps)

    # Get truth spectrum and sed
    tspec, tphot, ttheta = true_sed(mod, obsdat, sps=sps, fullspec=True)
    twave = sps.ssp.wavelengths
    tspec /= obsdat['normalization_guess']
    
    # Get calibration vector
    calvec = obsdat['calibration']
    if np.size(calvec) == 1:
        calvec = np.zeros(len(mwave)) + calvec
    else:
        calvec = calvec[obsdat['mask']]
    norm = obsdat['normalization_guess'] * obsdat['rescale']

    ### Build Figure ###
    fig = pl.figure(figsize=(10,8))
    # Set up left hand side
    gs1 = gridspec.GridSpec(3, 1)
    gs1.update(left=0.05, right=0.45, wspace=0.05)
    oax = pl.subplot(gs1[0,0])
    cax = pl.subplot(gs1[1:,0])

    # Set up right hand side
    gs2 = gridspec.GridSpec(4, 2)
    gs2.update(left=0.55, right=0.98, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])

    # Calibration figure
    cfig, cax = calfig(mwave, calvec, specvecs, norm=norm, fax=(None, cax),
                       basecolor=bcolor)
    cax = format_calax(cax, norm)
    
    # Intrinsic SED figure
    sfig, sax = sedfig(fwave, fspecvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax),
                       peraa=True, basecolor=bcolor, pointcolor=pcolor)
    sax.plot(twave, tspec, color='black', label='True spectrum', alpha=0.6)
    sax = format_sedax(sax)
    
    # Observed spectrum figure
    ofig, oax = residualfig(mwave, mospec, specvecs, unc=mounc,
                            fax=(None,oax), chi=True, basecolor=bcolor)
    oax.set_ylabel(r'$\chi$')
    oax.legend(loc=0, prop={'size':8})
    oax.set_ylim(-1,1)

    # Posterior parameter histograms
    pnames = ['mass', 'tage', 'dust2', 'zmet']
    samples, pnames_ord = hist_samples(res, mod, pnames, start=0.5, thin=10)
    truths = [obsdat['mock_params'][k] for k in pnames_ord]
    hfig, haxes = histfig(samples, pnames_ord, truths = truths, pname_map=pnmap,
                          basecolor=bcolor, fax=(None, haxes))
    haxes[0].legend(loc=0, prop={'size':8})

    # Save
    fig.savefig(resfile.replace('_mcmc','.dashboard.pdf'))
