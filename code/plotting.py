import numpy as np
import matplotlib.pyplot as pl
import bsfh.read_results as bread
from copy import deepcopy
import matplotlib

newcolors = {'blue':u'#5DA5DA',
             'orange': u'#FAA43A',
             'green': u'#60BD68',
             'red': u'#F15854',
             'magenta':u'#FF70B1', #actually pink
             'cyan': u'#39CCCC', #actually teal
             'yellow': u'#FFDC00',
             'maroon': u'#85144B',
            }
matplotlib.colors.cnames.update(newcolors)

def comp_samples(thetas, model, obs, sps=None, gp=None):
    """Different components of the model for a given set of thetas.

    :param thetas:
        A list or iterable of theta vectors for which model components
        are desired.
    :returns wave:
        The full wavelength array
    :returns mospec:
        The observed spectrum, linear units
    :returns mounc:
        The observational errors, linear units
    :returns specvecs:
        A list length len(theta) where each element is a list of model
        components corresponding to theta.  The model components are:
    """
    logarithmic = obs.get('logify_spectrum', True)
    specvecs = []
    wave, ospec, mask = obs['wavelength'], obs['spectrum'], obs['mask']
    mwave, mospec = wave[mask], ospec[mask]
    mounc = obs['unc'][mask]
    unc = mounc
    gp.wave, gp.sigma = mwave, obs['unc'][mask]
    if logarithmic:
        mospec = np.exp(mospec)
        unc = mounc * mospec

    for theta in thetas:
        mu, cal, delta, mask, wave = bread.model_comp(theta, model, obs, sps,
                                                      gp=gp, photflag=0)
        
        if logarithmic:
            full_cal = np.exp(np.log(cal) + delta)
            mod = np.exp(np.log(mu) + np.log(cal) + delta)
            residual = np.exp(np.log(mospec) - np.log(mod))
            chi =  (np.log(mospec)-np.log(mod)) / mounc
        else:
            full_cal = cal + delta/mu
            mod = cal * mu + delta
            residual  = mospec - mod
            chi = residual / mounc
        specvecs += [ [mu, full_cal, delta, mod,
                       residual, chi, cal]]
            
    return wave, mospec, unc, specvecs

def comp_samples_phot(thetas, model, obs, sps=None):
    """Different components of the model for a given set of thetas,
    for the photometry """
    specvecs = []
    wave = np.array([f.wave_effective for f in obs['filters']])
    mask = obs['phot_mask']
    mospec = obs['maggies'][mask]
    mounc = obs['maggies_unc'][mask]
    zero = np.zeros(mask.sum())
    
    for theta in thetas:
        mu = model.mean_model(theta, obs, sps=sps)[1][mask]
        specvecs += [ [mu, zero, zero, mu, mospec - mu, (mospec - mu)/mounc] ]
    return wave[mask], mospec, mounc, specvecs

def comp_samples_fullspec(thetas, model, obs, sps=None, gp=None):
    """Different components of the model for a given set of thetas,
    for the full spectrum (not just the wavelengths in obs).
    """
    
    specvecs = []
    mospec, mounc = None, None
    fullobs = deepcopy(obs)
    fullobs['wavelength'] = sps.ssp.wavelengths.copy()
    fullobs['mask'] = np.ones( len(fullobs['wavelength']), dtype= bool)
    for theta in thetas:
        mu, cal, delta, mask, wave = bread.model_comp(theta, model, fullobs, sps,
                                                      gp=None, photflag=0)
        cal = np.exp(cal)
        full_cal = np.exp(np.log(cal) + delta)
        mod = np.exp(np.log(mu) + np.log(cal) + delta)
        specvecs += [[mu, cal, delta, mod]]
        
    return wave, mospec, mounc, specvecs

def true_sed(model, obs, sps=None, fullspec=False):
    """The true SED for something with known parameters (stored in
    obs['mock_params'])
    """
    fullobs = deepcopy(obs)
    mockpars = deepcopy(obs['mock_params'])
    model.params = mockpars
    theta = model.theta.copy()
    if fullspec:
        fullobs['wavelength'] = sps.ssp.wavelengths.copy()
        fullobs['mask'] = np.ones( len(fullobs['wavelength']), dtype= bool)
        
    mu, phot, x = model.sed(theta, fullobs, sps=sps)
    return mu, phot, theta
        
def theta_samples(res, samples=[1.0], start=0.0, thin=1):

    nw, niter = res['chain'].shape[:-1]
    start_index = np.floor(start * (niter-1)).astype(int)
    flatchain = res['chain'][:,start_index::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    ns = flatchain.shape[0]
    thetas = [flatchain[s,:] for s in np.floor(np.array(samples) * (ns-1)).astype(int)]
    return thetas, start_index, np.floor(np.array(samples) * (ns-1)).astype(int)

def calfig(wave, calvec, specvecs, obsvec=None, norm=1.0,
           mlabel='Mock Truth', fax=None, rescale=1,
           basecolor='green', caltype='full', **extras):
    """Plot the calibration and posterior samples of it.

    :param caltype: (default: 'full')
        The calibration type to plot.  One of 'full'|'gp'|'poly'|'total'
           
    """
    cfig, cax = fax
    if caltype == 'full':
        samples = [norm * spec[1] for spec in specvecs]
        cplot = True
    elif caltype == 'gp':
        samples = [norm * spec[2]/spec[0] for spec in specvecs]
        cplot = False
    elif caltype == 'poly':
        samples = [norm * spec[6] for spec in specvecs]
        cplot = True
    elif caltype == 'total':
        samples = [norm * obsvec/spec[0] for spec in specvecs]
        cplot = True

    if cplot:
        #plot the calibration vector 
        cax.plot(wave, calvec/rescale, color='black', label=mlabel,
                 linewidth=3.0)
    # and posterior samples of it
    for i, specs in enumerate(samples):
        if i==0:
            label = 'Posterior sample'
        else:
            label = None
        cax.plot(wave, specs/rescale, label=label,
                 color=basecolor, alpha=0.3)
    
    return cfig, cax

def obsfig(wave, obsvec, specvecs, unc=None,
           labelprefix='Mock Observed', fax=None):
    """Plot the observed spectrum and posterior samples of it.
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
    oax.plot(wave, obsvec, color='black', label=labelprefix +' Spectrum',
             linewidth=1.0, alpha=1.0)
    return ofig, oax

def residualfig(wave, obsvec, specvecs, unc=None, chi=False,
                basecolor=None, fax=None):
    if fax is None:
        rfig, rax = pl.subplots()
    else:
        rfig, rax = fax

    if chi:
        #normalize residuals by the uncertainty
        chi_unc = unc.copy()
    else:
        chi_unc = 1.0
        if unc is not None:
            x, y, e = wave, obsvec, unc
            rax.fill_between(x, -e, e, facecolor='grey', alpha=0.3)
    # Plot posterior samples of the observed spectrum as a residual
    for i, specs in enumerate(specvecs):
        if i==0:
            label = 'Posterior samples'
        else:
            label = None
        rax.plot(wave, (specs[3] - obsvec) / chi_unc, linewidth=0.5,
                 color=basecolor, alpha=0.3, label=label)
    rax.axhline(0.0, linestyle=':', color='k')
    return rfig, rax
                
        
def sedfig(wave, specvecs, phot, photvecs, norm = 1.0,
            labelprefix='Mock', fax=None, peraa=False,
            basecolor='green', pointcolor='magenta', **kwargs):
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
    ylabel = r'$\lambda F_{{\lambda}}$ (intrinsic, cgs)'
    if peraa:
        sconv /= wave
        pconv /= pwave
        ylabel = r'$F_{{\lambda}}$ (intrinsic, cgs)'
    for i, (specs, seds) in enumerate(zip(specvecs, photvecs)):
        if i==0:
            label = 'Posterior samples'
        else:
            label = None
        sax.plot(wave, specs[0] * sconv, linewidth=0.5,
                 color=basecolor, alpha=0.3, label=label)
        sax.plot(pwave, seds[0] * pconv, markersize=8.0, linestyle='',
                 marker='o', alpha=0.5, mec=pointcolor, color=pointcolor,
                 label=label)

    sax.errorbar(pwave, sed * pconv, yerr=sed_unc * pconv,
                 marker='o', markersize=5.0, ecolor='black',
                 color='white', alpha=1.0, mec='black', mew=2,
                 linestyle='', label=labelprefix+' Photometry')
    sax.set_ylabel(ylabel)
    return sfig, sax

def hist_samples(res, model, showpars, start=0, thin=1,
                 return_lnprob=False, **extras):
    """Get samples for the parameters listed in showpars.
    """
    nw, niter = res['chain'].shape[:-1]
    parnames = np.array(model.theta_labels())
    start_index = np.floor(start * (niter-1)).astype(int)
    flatchain = res['chain'][:,start_index::thin,:]
    dims = flatchain.shape[0], flatchain.shape[1], flatchain.shape[2]
    flatchain = flatchain.reshape(dims[0]*dims[1], dims[2])
    ind_show = np.array([p in showpars for p in parnames], dtype= bool)
    flatchain = flatchain[:, ind_show]
    if return_lnprob:
        flatlnprob = res['lnprobability'][:, start_index::thin].reshape(dims[0]*dims[1])
        return flatchain, parnames[ind_show], flatlnprob
    
    return flatchain, parnames[ind_show]

def histfig(samples, parnames, truths=None, fax=None,
            truth_color='k', basecolor='green', pname_map={},
            **kwargs):
    """Plot a histogram of the given samples in each parameter.
    """
    npar = len(parnames)

    if fax is None:
        nx = int(np.floor(np.sqrt(npar)))
        ny = int(np.ceil(npar*1.0/nx))
        hfig, haxes = pl.subplots(nx, ny)
    else:
        hfig, haxes = fax
        
    for i, (ax, name) in enumerate(zip(haxes.flatten(), parnames)):
        if name == 'mass':
            rescale = 1e5
        else:
            rescale = 1
        ax.hist(samples[:,i]/rescale, bins=kwargs.get("bins", 50),
                histtype="stepfilled",
                color=basecolor,
                alpha = 0.5,
                label = 'posterior PDF')
        if truths is not None:
            ax.axvline(truths[i]/rescale, color=truth_color, label='Mock Truth', linestyle='--')
        ax.set_xlabel(pname_map.get(name, name), fontsize=6)
        ax.set_yticklabels([])
        pl.setp(ax.get_xticklabels(), fontsize=6)
    return hfig, haxes

def deltafig(results, models, pnames, pmap={}, plabel_map={},
             sfraction=0.9, thin=10, ptile=[16, 50, 84],
             fax=None, colors=None, labels=None, verbose=False,
             fractional=False, **kwargs):
    """Make plots of the fractional uncertainty in parameters for
    different runs.
    """
    if fax is None:
        pfig, pax = pl.subplots()
    else:
        pfig, pax = fax
    if labels is None:
        labels = len(results) * ['']
    if colors is None:
        colors = len(results) * ['k']

    def identity(x):
        return x
        
    for res, mod, label, clr in zip(results, models, labels, colors):
        samples, pnames_ord = hist_samples(res, mod, pnames,
                                           start=sfraction, thin=thin)
        truths = np.array([res['obs']['mock_params'][k][0] for k in pnames_ord])
        pct = np.percentile(samples, ptile, axis=0)
        
        pct = np.array([pmap.get(k, identity)(pct[:,i])
                        for i,k in enumerate(pnames_ord)])
        truths = np.array([pmap.get(k, identity)(truths[i])
                           for i,k in enumerate(pnames_ord)])
            
        delta = (pct - truths[:, None])
        if fractional:
            delta /= truths[:, None]
        pax.plot(delta[:,1], '-o', color=clr, label=label)
        pax.fill_between(np.arange(len(pnames)), delta[:,0], delta[:,2],
                         alpha=0.3, color=clr)

    pretty_pname = [plabel_map.get(k, k) for k in pnames_ord]
    pax.set_xticks(np.arange(len(pnames_ord)))
    pax.set_xticklabels(pretty_pname)
    pax.set_ylabel(r'(Post-True) / True')
    return pfig, pax
    
def deltafig_vspar(results, models, pname, pvary, pmap={},
                   sfraction=0.9, thin=10, ptile=[16, 50, 84],
                   fax=None, color='k', label='', xlims=None,
                   fractional=False, verbose=False, **kwargs):
    """Make plots of the (fractional) uncertainty in a certain
    parameter for different runs.
    """
    if fax is None:
        pfig, pax = pl.subplots()
    else:
        pfig, pax = fax
    
    def identity(x):
        return x

    pct = np.zeros([ len(results), len(ptile)])
    truths = np.zeros( len(results) )
    vary_param = np.zeros( len(results) )
    for ir, (res, mod) in enumerate(zip(results, models)):
        samples, pname_o = hist_samples(res, mod, pname,
                                        start=sfraction, thin=thin)
        if verbose:
            print(pname, ir)
        truths[ir] = res['obs']['mock_params'][pname][0]
        vary_param[ir] = res['obs']['mock_params'][pvary][0]
        pct[ir,:] = np.percentile(samples, ptile)
        # Transform the values if pmap supplied for this param
        pct[ir,:] = pmap.get(pname, identity)(pct[ir,:])
        truths[ir] = pmap.get(pname, identity)(truths[ir])
        
    delta = (pct - truths[:, None])
    if fractional:
        delta /= truths[:, None]
        #title += '$/$Truth'
    pax.plot(vary_param, delta[:,1], '-o', color=color, label=label)
    pax.fill_between(vary_param, delta[:,0], delta[:,2],
                     alpha=0.3, color=color)
    pax.set_xlim(xlims)
    
    return pfig, pax


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

def format_calax(cax, norm, rescale=1):
    if rescale == 1:
        label = r'Calibration [Counts/Flux (cgs)]'
    else:
        label = r'Calibration [Counts/Flux (cgs)]' + '($\\times 10^{{{0:2.0f}}}$)'.format(np.log10(1/rescale))
    cax.set_ylabel(label, fontsize=12)
    cax.legend(loc=0, prop={'size':8})
    cax.set_ylim(0.2*norm/rescale, 1.6*norm/rescale)

    #cax.set_xlim(3e3, 1.6e4)
    #ticks = list(cax.get_xlim()) + [4e3, 6e3, 10e3]
    #cax.set_xticks(ticks)
    #cax.set_xticklabels(['{:4.0f}'.format(t) for t in ticks], fontsize=8)
    cax.set_xlabel('$\lambda (\AA)$', fontsize=12)
    cax.tick_params(axis='both', which='major', labelsize=8)
    cax.set_ylabel(cax.get_ylabel(), fontsize=12)
    cax.legend(loc=0, prop={'size':12})
    return cax
