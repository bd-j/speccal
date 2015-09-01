import numpy as np
import pickle, os
from bsfh import priors, sedmodel
from sedpy import attenuation
from ggcdata import *

#The speccal directory
sdir = os.path.join(os.environ['PROJECTS'], 'speccal')

run_params = {'verbose':True,
              'outfile':'results/ggc_ngc1851',
              'do_powell': False,
              'ftol':0.5e-4, 'maxfev':5000,
              'nwalkers':64,
              'nburn':[64, 64, 128, 128, 256], 'niter':512,
              'initial_disp':0.3,
              'debug':False,
              'logify_spectrum':False,
              'normalize_spectrum':True,
              'norm_band_name':'sdss_g0',
              'rescale':True,
              'objname':'NGC1851',
              'wlo':3350.,
              'whi':6500.,
              'noisefactor':1.0,
              'mask':True,
              'broaden_obs':True,
              'calibrated':True
              }

##### OBSERVATIONAL DATA ######

def load_obs(objname=None, noisefactor=1.0, calibrated=True,
             mask=True, broaden_obs=False, **extras):
    # get the spectrum
    if calibrated:
        fluxtype=None
    else:
        fluxtype=1
    obs = ggc_spec(objname, 'a', '1', fluxtype=fluxtype,
                   datadir=os.path.join(sdir, 'data/ggclib/spectra'), **extras)
    #mask the spectrum
    if mask:
        obs = ggc_mask(obs, lsf=6, pad=6, thresh=3, width=20)
    #broaden to constant FWHM in wave space.
    if broaden_obs:
        smask = obs.get('mask', np.ones(len(obs['wavelength']), dtype=bool))
        bflux = broaden_ggcspec(obs['wavelength'][smask], obs['spectrum'][smask])
        obs['spectrum'][smask] = bflux
    #add the photometric data
    obs.update(ggc_phot(objname, datadir=os.path.join(sdir, 'data/ggclib/photometry')))
    obs['phot_mask'] = np.array(['sdss' in filt.name for filt in obs['filters']])
    #adjust uncertainties
    obs['unc'] *= noisefactor
    obs['noisefactor'] = noisefactor
    obs['spec_calibrated'] = calibrated
    return obs

obs = load_obs(**run_params)

###### MODEL ###########
model_type = sedmodel.SedModel

model_params = []

###### Distance ##########
model_params.append({'name': 'lumdist', 'N': 1,
                     'isfree': False,
                     'init': 0.011,
                     'units': 'Mpc',
                     'prior_function': None,
                     'prior_args': None})

###### SFH ################

model_params.append({'name': 'mass', 'N': 1,
                     'isfree': True,
                     'init': 2e5,
                     'units': r'M$_\odot$',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini':1e4, 'maxi': 1e7}})

model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init': 10.0,
                        'units': 'Gyr',
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':5.0, 'maxi':15.0}})

model_params.append({'name': 'zmet', 'N': 1,
                        'isfree': True,
                        'init': -0.5,
                        'units': r'$\log (Z/Z_\odot)$',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-2, 'maxi':0.19}})

model_params.append({'name': 'sfh', 'N':1,
                        'isfree': False,
                        'init': 0,
                        'units': None})

###### DUST ##################

model_params.append({'name': 'dust_curve', 'N': 1,
                        'isfree': False,
                        'init': attenuation.cardelli,
                        'units': None})

model_params.append({'name': 'dust2', 'N': 1,
                        'isfree': True,
                        'init': 1.0,
                        'units': r'$\tau_V$',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':2.5}})

model_params.append({'name': 'dust1', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': r'$\tau_V$',
                        'prior_function': None,
                        'prior_args': None})

model_params.append({'name': 'dust_tesc', 'N': 1,
                        'isfree': False,
                        'init': 0.01,
                        'units': 'Gyr',
                        'prior_function': None,
                        'prior_args': None})

model_params.append({'name': 'dust_type', 'N': 1,
                        'isfree': False,
                        'init': 1,
                        'units': 'NOT USED'})

###### IMF ###################

model_params.append({'name': 'imf_type', 'N': 1,
                        'isfree': False,
                        'init': 2, #2 = kroupa, #0=salpeter
                        'units': None})

model_params.append({'name': 'imf3', 'N':1,
                        'isfree': False,
                        'init': 2.35,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':1.3, 'maxi':3.3}})

###### WAVELENGTH SCALE ######

model_params.append({'name': 'zred', 'N':1,
                        'isfree': True,
                        'init': 0.000001,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-0.001, 'maxi':0.001}})

model_params.append({'name': 'sigma_smooth', 'N': 1,
                        'isfree': True,
                        'init': 1.0,
                        'units': r'$\AA$',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':3}})
                        #'prior_function': priors.lognormal,
                        #'prior_args': {'log_mean':np.log(2.2)+0.05**2, 'sigma':0.05}})

model_params.append({'name': 'smooth_velocity', 'N': 1,
                        'isfree': False,
                        'init': False,
                        'units': None})

model_params.append({'name': 'min_wave_smooth', 'N': 1,
                        'isfree': False,
                        'init': 3200.0,
                        'units': r'$\AA$'})

model_params.append({'name': 'max_wave_smooth', 'N': 1,
                        'isfree': False,
                        'init': 7000.0,
                        'units': r'$\AA$'})

#model.params.append({'name': 'lsf', 'N':1,
#                         'isfree':False,
#                         'init': line_spread_function,
#                         'units': None})
                         
###### CALIBRATION ###########

polyorder = 2
polymin = [-5e1, -1e2]
polymax = [5e1, 1e2]
polyinit = [0.01, 0.01]

model_params.append({'name': 'poly_coeffs', 'N': polyorder,
                        'isfree': True,
                        'init': polyinit,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':polymin, 'maxi':polymax}})
    
model_params.append({'name': 'spec_norm', 'N':1,
                        'isfree': True,
                        'init':0.0001,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-1.0, 'maxi':1.0}})

# This is for use with the matern branch for bsfh, where sqrt(jitter)
# now multiplies the noise (instead of adding)
model_params.append({'name': 'gp_jitter', 'N':1,
                        'isfree': True,
                        'init': 1.0,
                        'units': 'fractional noise squared',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.01, 'maxi':100}})

model_params.append({'name': 'gp_amplitude', 'N':1,
                        'isfree': True,
                        'init': 0.04,
                        'units': 'spec units',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.02, 'maxi':0.1}})

model_params.append({'name': 'gp_length', 'N':1,
                        'isfree': True,
                        'init': 300.0,
                        'units': r'$\AA$',
#                        'prior_function': priors.lognormal,
#                        'prior_args': {'log_mean':np.log(200.0)+0.5**2, 'sigma':0.5}})
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':250.0, 'maxi':350.0}})

model_params.append({'name': 'phot_jitter', 'N':1,
                        'isfree': True,
                        'init': 0.01,
                        'units': 'mags',
                        'prior_function': priors.logarithmic,
                        'prior_args': {'mini':0.0, 'maxi':0.05}})
