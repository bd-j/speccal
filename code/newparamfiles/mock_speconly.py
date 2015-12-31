import numpy as np
import pickle, os
from bsfh.source_basis import StellarPopBasis
from bsfh.models import SedModel, priors
from sedpy import attenuation

# ------------------
# Global parameters
# ------------------
sdir = os.path.join(os.environ['PROJECTS'], 'speccal')

run_params = {'verbose':True,
              'debug':False,
              'outfile':'results/ggc_mock_u0_t9.0_z0.0_a0.5',
              # Sampling parameters
              'do_powell': False,
              'ftol':0.5e-4, 'maxfev':5000,
              'nwalkers':64, 
              'nburn':[64, 64, 128, 128, 256], 'niter':512,
              'initial_disp':0.1,
              # Observational Data Manipulation
              'logify_spectrum':False,
              'normalize_spectrum':True,
              'norm_band_name':'sdss_g0',
              'rescale':True,
              # mock parameters
              'mass': 1e5,
              'tage': 9.0,
              'logzsol': 0.0,
              'zred': 0.0,
              'sigma_smooth':1.0,
              # mock noise parameters
              # Object to use for wavelength, SNR, and cal vectors
              'objname': 'NGC1851',
              'phot_snr': 20.0,
              'spec_snr': None, # Use GGC example
              'noisefactor': 1.0, # scale GGC SNR
              'add_noise': False,
              'calibrated': True,
              'wlo':3990.,
              'whi':6150.,
              # SPS basis parameters
              'zcontinuous': 1,
              }

# ------------------
# (MOCK) OBSERVATIONAL DATA
# ------------------

mockpars = ['mass', 'tage', 'logzsol', 'zred', 'sigma_smooth', 'lumdist']

def load_obs(mockpars=mockpars, **kwargs):

    filters = [kwargs['norm_band_name']]
    
    import mocks
    params = {}
    for p in mockpars:
        try:
            params[p] = np.array(kwargs[p])
        except:
            pass
    
    sps = load_sps(**kwargs)
    mod = load_model(**kwargs)
    mock = mocks.make_mock(mockparams=params, model=mod, sps=sps, filters=filters,
                           datadir=os.path.join(sdir, 'data/ggclib/spectra'), **kwargs)

    return mock

# ------------------
# MODEL
# ------------------
    
def load_sps(**kwargs):
    sps = StellarPopBasis(**kwargs)
    return sps

def load_gp(**kwargs):
    from bsfh.likelihood import gp
    spec_gp = gp.ExpSquared(None, None)
    phot_gp = gp.PhotOutlier(None, None)
    return spec_gp, phot_gp

model_params = []

# ------- Distance ---------------
model_params.append({'name': 'lumdist', 'N': 1,
                     'isfree': False,
                     'init': 0.01,
                     'units': 'Mpc',
                     'prior_function': None,
                     'prior_args': None})

# --------- SFH ------------------

model_params.append({'name': 'mass', 'N': 1,
                     'isfree': True,
                     'init': 1e5,
                     'init_disp': 1e4,
                     'units': r'M$_\odot$',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini':1e4, 'maxi': 1e6}})

model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init': 5.0,
                        'init_disp': 2.0,
                        'units': 'Gyr',
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':0.1, 'maxi':15.0}})

model_params.append({'name': 'zmet', 'N': 1,
                        'isfree': True,
                        'init': -0.5,
                        'init_disp': 0.2,
                        'units': r'$\log (Z/Z_\odot)$',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-2, 'maxi':0.19}})

model_params.append({'name': 'sfh', 'N':1,
                        'isfree': False,
                        'init': 0,
                        'units': None})

# ---------- Dust -----------------------

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

# --------- IMF -------------------

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

# ------- Wavelength Scale ---------

model_params.append({'name': 'zred', 'N':1,
                        'isfree': True,
                        'init': 1e-6,
                        'init_disp': 1e-5,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-0.001, 'maxi':0.001}})

model_params.append({'name': 'sigma_smooth', 'N': 1,
                        'isfree': True,
                        'init': 1.0,
                        'init_disp': 0.1,
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

# --------- Calibration -------------

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

model_params.append({'name': 'gp_jitter', 'N':1,
                        'isfree': True,
                        'init': 1.3,
                        'units': 'fractional noise squared',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':1.0, 'maxi':10.0}})

model_params.append({'name': 'gp_jitter_add', 'N':1,
                        'isfree': True,
                        'init': 1e-4,
                        'units': 'noise squared',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':1e-8, 'maxi':2.5e-3}})

model_params.append({'name': 'gp_amplitude', 'N':1,
                        'isfree': True,
                        'init': 0.01,
                        'init_disp': 0.005,
                        'units': 'spec units',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':0.25}})

model_params.append({'name': 'gp_length', 'N':1,
                        'isfree': True,
                        'init': 100.0,
                        'init_disp': 50.0,
                        'units': r'$\AA$',
                        'prior_function': priors.lognormal,
                        'prior_args': {'log_mean':np.log(100.0)+0.75**2, 'sigma':0.75}})
                        #'prior_function':priors.tophat,
                        #'prior_args': {'mini':10.0, 'maxi':1000}})

model_params.append({'name': 'phot_jitter', 'N':1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'mags',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':0.1}})


def load_model(**kwargs):
    # adjust model parameters based on command line arguments
    # initialize the model
    model = SedModel(model_params)
    return model
