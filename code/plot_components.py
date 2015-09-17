import sys, pickle, matplotlib
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from bsfh import sps_basis
import bsfh.gp

from plotting import *
matplotlib.colors.cnames.update(newcolors)

sps = sps_basis.StellarPopBasis()
gp = bsfh.gp.ExpSquared(None, None)


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
    mod._has_parameter_dependencies = False
    obsdat = res['obs']
    if 'gp_jitter_add' not in obsdat['mock_params']:
        obsdat['mock_params'].update({'gp_jitter_add': np.array([0.0])})
    if res['run_params'].get('gp_type', 'ExpSquared') == 'Matern':
        gp = bsfh.gp.Matern(None, None)
    
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
    if norm < 1e10:
        rescale = 1
    else:
        rescale = 1e18

    
    ### Build Figure ###
    fig = pl.figure(figsize=(10,8))

    # Set up left hand side
    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(left=0.05, right=0.45, wspace=0.1)
    pax = pl.subplot(gs1[0,0])
    gax = pl.subplot(gs1[1,0])
    
    # Set up right hand side
    gs2 = gridspec.GridSpec(2, 1)
    gs2.update(left=0.55, right=0.98, wspace=0.1)
    fax = pl.subplot(gs2[0,0])
    rax = pl.subplot(gs2[1,0])


    # Polynomial
    pfig, pax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       mlabel='Input calibration curve', rescale=rescale,
                       fax=(None, pax), basecolor=bcolor, caltype='poly')
    pax.set_ylabel('$F_{obs}/F_{intrinsic}$')
    pax.legend(loc='lower right', prop={'size':12})
    pax = format_calax(pax, norm=norm, rescale=rescale)
    pax.text(0.1, 0.8, r'Smooth ($e^{f(\alpha)}$)', transform=pax.transAxes,
             fontsize=12, color='red')
        
    # GP
    gfig, gax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       rescale=rescale,
                       fax=(None, gax), basecolor=bcolor, caltype='gp')
    gax.text(0.1, 0.2, r'GP ($\tilde{\Delta}/\mu$)', transform=gax.transAxes,
             fontsize=12, color='red')
    gax.set_ylabel('$\delta(F_{obs}/F_{intrinsic})$')
    #gax = format_calax(gax, norm=norm, rescale=rescale)
    gax.set_xlabel('$\lambda (\AA)$', fontsize=12)
    gax.tick_params(axis='both', which='major', labelsize=8)
    gax.set_ylabel(gax.get_ylabel(), fontsize=12)
    
    # Full
    ffig, fax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       mlabel='Input calibration curve', rescale=rescale,
                       fax=(None, fax), basecolor=bcolor, caltype='full')
    fax = format_calax(fax, norm=norm, rescale=rescale)
    fax.text(0.1, 0.8, r'Sum ($e^{f(\alpha)} + \tilde{\Delta}/\mu$)',
             transform=fax.transAxes, fontsize=12, color='red')
    #fax.set_ylabel('$F_{obs}/F_{intrinsic}$')
    
    # Residual
    for spec in specvecs:
        sfull = norm * spec[1]
        rax.plot(mwave, calvec / sfull, color=bcolor, alpha=0.3)
    rax.axhline(1.0, linestyle=':', color='black')
    rax.text(0.1, 0.2, r'$[e^{f(\alpha)} + \tilde{\Delta}/\mu]/Truth$',
             transform=rax.transAxes, fontsize=12, color='red')
    rax.set_ylim(0.8, 1.2)
    rax.set_xlabel('$\lambda (\AA)$', fontsize=12)
    rax.tick_params(axis='both', which='major', labelsize=8)
    rax.set_ylabel(rax.get_ylabel(), fontsize=12)


        #fig.show()
    fig.savefig(resfile.replace('_mcmc','.components.pdf'))
