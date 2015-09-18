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
    #add the photometric data
    obs.update(ggc_phot(objname, datadir=os.path.join(sdir, 'data/ggclib/photometry')))
    obs['phot_mask'] = np.array(['sdss' in filt.name for filt in obs['filters']])

    obs['spectrum'] = None
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
                        'isfree': False,
                        'init': 0.000000,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-0.001, 'maxi':0.001}})

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


model_params.append({'name': 'phot_jitter', 'N':1,
                        'isfree': True,
                        'init': 0.01,
                        'units': 'mags',
                        'prior_function': priors.logarithmic,
                        'prior_args': {'mini':0.001, 'maxi':0.05}})