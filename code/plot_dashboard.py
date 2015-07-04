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

if __name__ == "__main__":

    bcolor='magenta'
    if len(sys.argv) > 1:
        resfile=sys.argv[1]
        try:
            bcolor = sys.argv[2]
        except:
            pass
        try:
            suptitle = sys.argv[3]
        except:
            suptitle=''
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
    # Set up right hand side
    gs1 = gridspec.GridSpec(3, 1)
    gs1.update(left=0.55, right=0.98, wspace=0.05, hspace=0.01)
    
    oax = pl.subplot(gs1[0,0])
    cax = pl.subplot(gs1[1:,0])

    # Set up left hand side
    gs2 = gridspec.GridSpec(4, 2)
    gs2.update(left=0.05, right=0.45, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])

    # Calibration figure
    if norm < 1e10:
        rescale = 1
    else:
        rescale = 1e18
    cfig, cax = calfig(mwave, calvec, specvecs, norm=norm, rescale=rescale,
                       fax=(None, cax), caltype='full', basecolor=bcolor)
    cax = format_calax(cax, norm, rescale=rescale)

    # Residual Figure    
    ofig, oax = residualfig(mwave, mospec, specvecs, unc=mounc,
                            basecolor=bcolor,fax=(None, oax), chi=True)
    oax.set_ylabel(r'$\chi$')
    oax.legend(loc=0, prop={'size':8})
    oax.tick_params(axis='both', which='major', labelsize=8)
    oax.set_xticklabels([''])
    oax.set_ylabel(cax.get_ylabel(), fontsize=12)
    oax.set_ylim(-3,3)
    
        
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
    fig.suptitle(suptitle)
    fig.savefig(resfile.replace('_mcmc','.dashboard.pdf'))
