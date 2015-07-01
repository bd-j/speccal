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
    gs1.update(left=0.05, right=0.45, wspace=0.1)
    sax = pl.subplot(gs1[0,0])
    pax = pl.subplot(gs1[1,0])
    gax = pl.subplot(gs1[2,0])
    
    # Set up right hand side
    gs2 = gridspec.GridSpec(3, 1)
    gs2.update(left=0.55, right=0.98, wspace=0.1)
    fax = pl.subplot(gs2[0,0])
    tax = pl.subplot(gs2[1,0])
    rax = pl.subplot(gs2[2,0])

    #Instrinsic
    sfig, sax = sedfig(fwave, fspecvecs, [pwave, mosed, mosed_unc], pvecs,
                       norm=1/obsdat['normalization_guess'], fax=(None,sax),
                       peraa=True, basecolor=bcolor, pointcolor=pcolor)
    sax.plot(twave, tspec, color='black', label='True spectrum', alpha=0.6)
    sax = format_sedax(sax)
    sax.text(0.1, 0.8, r'Intrinsic ($\mu$)', transform=sax.transAxes,
             fontsize=12, color='red')

    # Polynomial
    pfig, pax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       mlabel='Input calibration curve',
                       fax=(None, pax), basecolor=bcolor, caltype='poly')
    pax.text(0.1, 0.8, r'Polynomial ($e^{f(\alpha)}$)', transform=pax.transAxes,
             fontsize=12, color='red')
    pax.set_ylabel('$F_{obs}/F_{intrinsic}$')
    sax.set_xlim(pax.get_xlim())
    pax.legend(loc='lower right', prop={'size':12})
    
    # GP
    gfig, gax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       fax=(None, gax), basecolor=bcolor, caltype='gp')
    gax.text(0.1, 0.8, r'GP ($\tilde{\Delta}/\mu$)', transform=gax.transAxes,
             fontsize=12, color='red')
    gax.set_ylabel('$\delta(F_{obs}/F_{intrinsic})$')
    
    # Full
    ffig, fax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       mlabel='Input calibration curve',
                       fax=(None, fax), basecolor=bcolor, caltype='full')
    fax.text(0.1, 0.8, r'Full ($e^{f(\alpha)} + \tilde{\Delta}/\mu$)',
             transform=fax.transAxes, fontsize=12, color='red')
    fax.set_ylabel('$F_{obs}/F_{intrinsic}$')
    
    # Total
    tfig, tax = calfig(mwave, calvec, specvecs, norm=norm, obsvec=mospec,
                       mlabel='Input calibration curve',
                       fax=(None, tax), basecolor=bcolor, caltype='total')
    #cax = format_calax(cax, norm)
    tax.text(0.1, 0.8, r'Total ($F_{obs}/\mu$)',
             transform=tax.transAxes, fontsize=12, color='red')
    tax.set_ylabel('$F_{obs}/F_{intrinsic}$')

    # Residual
    rfig, rax = residualfig(mwave, mospec, specvecs, unc=mounc,
                            fax=(None, rax), chi=True, basecolor=bcolor)
    rax.set_ylabel(r'$\chi$')
    rax.legend(loc=0, prop={'size':8})
    rax.set_ylim(-1,1)
    rax.text(0.1, 0.8, r'Residual $(F_{obs} - [\mu e^{f(\alpha)} + \tilde{\Delta}])/\sigma$',
             transform=rax.transAxes, fontsize=12, color='red')

    #fig.show()
    fig.savefig(resfile.replace('_mcmc','.components.pdf'))
