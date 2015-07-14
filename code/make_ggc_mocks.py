import numpy as np
import pickle, os, sys
from bsfh import sps_basis, priors, sedmodel, elines
from sedpy import attenuation
import ggcdata

#The speccal directory
sdir = os.path.join(os.environ['PROJECTS'], 'speccal')

sps = sps_basis.StellarPopBasis()
nw = len(sps.ssp.wavelengths)
if nw > 7000:
    lib = 'ckc'
elif nw > 2000:
    lib = 'miles'
    sigma_library = 1.08
else:
    lib = 'basel'

model_params = []

###### Distance ##########
model_params.append({'name': 'lumdist', 'N': 1,
                     'isfree': False,
                     'init': 0.01,
                     'units': 'Mpc',
                     'prior_function': None,
                     'prior_args': None})

###### SFH ################

model_params.append({'name': 'mass', 'N': 1,
                     'isfree': False,
                     'init': 1e5,
                     'units': r'M$_\odot$',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini':1e4, 'maxi': 1e6}})

model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init': 12.0,
                        'units': 'Gyr',
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':0.1, 'maxi':15.0}})

model_params.append({'name': 'zmet', 'N': 1,
                        'isfree': True,
                        'init': 0.0,
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
                        'init': 0.5,
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
                        'init': 0.0000,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':-0.001, 'maxi':0.001}})

model_params.append({'name': 'sigma_smooth', 'N': 1,
                        'isfree': False,
                        'init': 1.3,
                        'units': r'$\AA$',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':3}})
                        #'prior_function': priors.lognormal,
                        #'prior_args': {'log_mean':np.log(2.2)+0.05**2, 'sigma':0.05}})

model_params.append({'name': 'lsf', 'N': 1,
                        'isfree': False,
                        'init': None,
                        #'init': ggcdata.spec_lsf,
                        'units': r'$\AA (dispersion)$',
                        'prior_function': None,
                        'prior_args': None})

                        
model_params.append({'name': 'smooth_velocity', 'N': 1,
                        'isfree': False,
                        'init': False,
                        'units': None})

model_params.append({'name': 'min_wave_smooth', 'N': 1,
                        'isfree': False,
                        'init': 3000.0,
                        'units': r'$\AA$'})

model_params.append({'name': 'max_wave_smooth', 'N': 1,
                        'isfree': False,
                        'init': 7500.0,
                        'units': r'$\AA$'})

###### CALIBRATION ###########

polyorder = 2
polymin = [-1e1, -1e1]
polymax = [1e1, 1e1]
polyinit = [0.0, 0.0]

model_params.append({'name': 'poly_coeffs', 'N': polyorder,
                        'isfree': False,
                        'init': polyinit,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':polymin, 'maxi':polymax}})
    
model_params.append({'name': 'spec_norm', 'N':1,
                        'isfree': False,
                        'init':1.00000,
                        'units': None,
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.2, 'maxi':5}})

model_params.append({'name': 'gp_jitter', 'N':1,
                        'isfree': False,
                        'init': 0.0000,
                        'units': 'spec units',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':0.2}})

model_params.append({'name': 'gp_amplitude', 'N':1,
                        'isfree': False,
                        'init': 0.0000,
                        'units': 'spec units',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':0.5}})

model_params.append({'name': 'gp_length', 'N':1,
                        'isfree': False,
                        'init': 60.0,
                        'units': r'$\AA$',
                        'prior_function': priors.lognormal,
                        'prior_args': {'log_mean':np.log(60.0)+1.0**2, 'sigma':1.0}})

model_params.append({'name': 'phot_jitter', 'N':1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'mags',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini':0.0, 'maxi':0.1}})

fd = {'tage':12.0, 'zmet':0.0, 'dust2':0.5}

def main(versiondir='', do_vary=False, do_noise=False,
         fiducial_params=fd, **kwargs):

        
    
    vary_params = {'tage': [0.3, 1.1, 3.0, 6.0, 9.0, 12.0],
                   'zmet': [-1.5, -1.0, -0.5, 0.0, 0.15],
                   'dust2': [0, 0.5, 1.0, 2.0]}
        
    info = {'objname': 'NGC1851',
            'datadir': os.path.join(sdir, 'data/ggclib/spectra/'),
            'outdir': os.path.join(sdir, 'data/ggclib/mocks/'),
            'apply_cal': True,
            'add_noise': False,
            'mask': True,
            'phot_snr':10}
        
    info.update(kwargs)
    caltype = ['c', 'u']
    noisetype = ['0','1']
    name_template = os.path.join(info['outdir'], lib,
                                 'ggc_mock.{0}{1}.t{2:3.1f}_z{3:3.1f}_a{4:3.1f}.pkl')
    model = sedmodel.SedModel(model_params)
    for k, v in fiducial_params.iteritems():
        model.params[k] = np.atleast_1d(np.array(v))
    theta_default = model.theta.copy()

    # Make calibrated mock
    info['apply_cal'] = False
    theta = theta_default.copy()
    mock = ggcdata.ggc_mock(model, theta, sps, **info)
    filename = name_template.format(caltype[info['apply_cal']], 0, *theta)
    print('writing to {0}'.format(filename))
    with open(filename, 'w') as f:
        pickle.dump(mock, f)

    # Make uncalibrated mock
    info['apply_cal'] = True
    theta = theta_default.copy()
    mock = ggcdata.ggc_mock(model, theta, sps, **info)
    filename = name_template.format(caltype[info['apply_cal']], 0, *theta)
    print('writing to {0}'.format(filename))
    with open(filename, 'w') as f:
        pickle.dump(mock, f)
   
    
    # Make uncalibrated mocks with varied parameters
    info['apply_cal'] = True
    if do_vary:
        for p, vals in vary_params.iteritems():
            ind = model.theta_labels().index(p)
            theta = theta_default.copy()
            for v in vals:
                theta[ind] = v
                mock = ggcdata.ggc_mock(model, theta, sps, **info)
                filename = name_template.format(caltype[info['apply_cal']], 0,
                                                *theta)
                print('writing to {0}'.format(filename))
                with open(filename, 'w') as f:
                    pickle.dump(mock, f)
    
    # Make mocks with different noise realizations
    info['apply_cal'] = False
    if do_noise:
        info['add_noise'] = True
        theta = theta_default.copy()
        nnoise = 10 #number of noise realizations
        for i in range(nnoise):
            mock = ggcdata.ggc_mock(model, theta, sps, **info)
            filename = name_template.format(caltype[info['apply_cal']], i+1,
                                            *theta)
            print('writing to {0}'.format(filename))
            with open(filename, 'w') as f:
                pickle.dump(mock, f)
        
    return fiducial_params['tage'], fiducial_params['zmet'], fiducial_params['dust2']

if __name__ == "__main__":
    fp = main()
    print(fp)
