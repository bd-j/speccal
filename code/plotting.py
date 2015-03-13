import numpy as np
import matplotlib.pyplot as pl
import bsfh.read_results as bread

def comp_samples(thetas, model, obs, sps=None, gp=None):
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
    mospec = obs['maggies']
    mounc = obs['maggies_unc']
    mask = obs['phot_mask']
    zero = np.zeros(mask.sum())
    
    for theta in thetas:
        mu = model.mean_model(theta, obs, sps=sps)[1][mask]
        specvecs += [ [mu, zero, zero, mu, mospec - mu, (mospec - mu)/mounc] ]
    return wave, mospec, mounc, specvecs

def theta_samples(res, samples=[1.0], start=0.0, thin=1):

    nw, niter = res['chain'].shape[:-1]
    start_index = np.floor(start * (niter-1)).astype(int)
    flatchain = res['chain'][:,start_index::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    ns = flatchain.shape[0]
    thetas = [flatchain[s,:] for s in np.floor(np.array(samples) * (ns-1)).astype(int)]
    return thetas, start_index, np.floor(np.array(samples) * (ns-1)).astype(int)
