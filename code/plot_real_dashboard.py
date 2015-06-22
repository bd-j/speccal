import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from bsfh import sps_basis
from bsfh.gp import ExpSquared

from plotting import *

sps = sps_basis.StellarPopBasis()
#gp = GaussianProcess(None, None)
try:
    import george
    kernel = (george.kernels.WhiteKernel(0.0) +
            0.0 * george.kernels.ExpSquaredKernel(0.0))
    gp = george.GP(kernel)#, solver=george.HODLRSolver)
except:
    pass

gp = ExpSquared(None, None)

if __name__ == "__main__":

    if len(sys.argv) > 1:
        resfile=sys.argv[1]
    else:
        resfile = 'results/ggc_mock.u0.t1.0_z0.0_a0.5_1426268715_mcmc'
    model_file = resfile.replace('_mcmc','_model')
    bcolor='green'
    pcolor='orange'
    nsample = 10
    
    res, pr, mod = bread.read_pickles(resfile, model_file=model_file)
    obsdat = res['obs']
    
    # Get samples and component spectra
    fsamples = np.random.uniform(0,1,nsample)
    thetas, start, samples = theta_samples(res, samples=fsamples, start=0.75, thin=1)
    fwave, _, _, fspecvecs = comp_samples_fullspec(thetas, mod, obsdat,
                                                   sps=sps, gp=None)
    mwave, mospec, mounc, specvecs = comp_samples(thetas, mod, obsdat,
                                                  sps=sps, gp=gp)
    pwave, mosed, mosed_unc, pvecs = comp_samples_phot(thetas, mod, obsdat, sps=sps)

    # Get truth spectrum and sed
    # THERE IS NO TRUTH
    
    # Get calibration vector    
    if obsdat.get('spec_calibrated', True):
        norm = 1
        calvec = np.ones(len(mwave))
    else:
        import ggcdata
        cal = ggcdata.ggc_spec(obsdat['object_name'], 'a', '1',
                               fluxtype=None)
        calvec = cal['calibration'].copy()
        calvec = calvec[obsdat['mask']]
        norm = 1e18#obsdat['normalization_guess'] * obsdat['rescale']
        
    ### Build Figure ###
    fig = pl.figure(figsize=(10,8))
    # Set up left hand side
    gs1 = gridspec.GridSpec(3, 1)
    gs1.update(left=0.05, right=0.48, wspace=0.05)
    oax = pl.subplot(gs1[0,0])
    cax = pl.subplot(gs1[1:,0])    
    # Set up right hand side
    gs2 = gridspec.GridSpec(4, 2)
    gs2.update(left=0.55, right=0.98, hspace=0.5)
    sax = pl.subplot(gs2[0:2,:])
    haxes = np.array([pl.subplot(gs2[i, j]) for i in [2,3] for j in [0,1]])

    # Calibration figure
    cfig, cax = calfig(mwave, calvec, specvecs, norm=norm,
                       labelprefix='Schiavon ', fax=(None, cax), basecolor=bcolor)
    cax.set_ylabel(r'Calibration $F_{{obs}}/F_{{\lambda, intrinsic}}$')
    cax.legend(loc=0, prop={'size':8})
    #cax.set_ylim(0.2, 1.6)
    # Intrinsic SED figure
    sfig, sax = sedfig(fwave, fspecvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax),
                       labelprefix='Observed', peraa=True, basecolor=bcolor, pointcolor=pcolor)
    sax.set_xscale('log')
    sax.set_xlim(2.5e3, 1.6e4)
    sax.legend(loc=0, prop={'size':12})

    # Residual Figure    
    ofig, oax = residualfig(mwave, mospec, specvecs, unc=mounc,
                            basecolor=bcolor,fax=(None,oax), chi=True)
    oax.set_ylabel(r'$\chi$')
    oax.legend(loc=0, prop={'size':8})
    oax.set_ylim(-1,1)

    # Posterior parameter histograms
    pnames = ['mass', 'tage', 'dust2', 'zmet']
    samples, pnames_ord = hist_samples(res, mod, pnames, start=0.75, thin=10)
    hfig, haxes = histfig(samples, pnames_ord, truths=None,
                          basecolor=bcolor, fax=(None, haxes))
    haxes[0].legend(loc=0, prop={'size':8})


    fig.savefig(resfile.replace('_mcmc','.dashboard.pdf'))
