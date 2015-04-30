import numpy as np
import matplotlib.pyplot as pl
import bsfh.read_results as bread
from copy import deepcopy

def comp_samples(thetas, model, obs, sps=None, gp=None):
    """Different components of the model for a given set of thetas.

    :returns wave:
        The full wavelength array
    :return 
    """
    specvecs = []
    wave, ospec, mask = obs['wavelength'], obs['spectrum'], obs['mask']
    mwave, mospec = wave[mask], ospec[mask]
    mounc = obs['unc'][mask]
    gp.wave, gp.sigma = mwave, obs['unc'][mask]
    mospec = np.exp(mospec)
         #mounc *= mospec

    for theta in thetas:
        mu, cal, delta, mask, wave = bread.model_comp(theta, model, obs, sps,
                                                      gp=gp, photflag=0)
        cal = np.exp(cal)
        full_cal = np.exp(np.log(cal) + delta)
        mod = np.exp(np.log(mu) + np.log(cal) + delta)
        #mu = np.exp(mu)
            
        specvecs += [ [mu, cal, delta, mod,
                       np.exp(np.log(mospec)-np.log(mod)),
                       (np.log(mospec)-np.log(mod)) / mounc] ]
            
    return wave, mospec, mounc, specvecs

def comp_samples_phot(thetas, model, obs, sps=None):
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
    fullobs = deepcopy(obs)
    mockpars = deepcopy(obs['mock_params'])
    model.params = mockpars
    theta = model.theta.copy()
    if fullspec:
        fullobs['wavelength'] = sps.ssp.wavelengths.copy()
        fullobs['mask'] = np.ones( len(fullobs['wavelength']), dtype= bool)
        
    mu, phot, x = model.sed(theta, fullobs, sps=sps)
    return mu, phot
        
def theta_samples(res, samples=[1.0], start=0.0, thin=1):

    nw, niter = res['chain'].shape[:-1]
    start_index = np.floor(start * (niter-1)).astype(int)
    flatchain = res['chain'][:,start_index::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    ns = flatchain.shape[0]
    thetas = [flatchain[s,:] for s in np.floor(np.array(samples) * (ns-1)).astype(int)]
    return thetas, start_index, np.floor(np.array(samples) * (ns-1)).astype(int)

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

def sedfig(wave, specvecs, phot, photvecs, norm = 1.0, fax=None, peraa=False):
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
        sax.plot(wave, specs[0] * sconv,
                 color='green', alpha=0.3, label=label)
        sax.plot(pwave, seds[0] * pconv, markersize=8.0, linestyle='',
                 marker='o', color='magenta', label=label)

    sax.errorbar(pwave, sed * pconv, yerr=sed_unc * pconv,
                 marker='o', markersize=8.0,
                 color='black', linestyle='', label='Mock Photometry')
    sax.set_ylabel(ylabel)
    return sfig, sax

def hist_samples(res, model, showpars, start=0, thin=1, **extras):
    
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
