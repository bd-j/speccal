import os
import numpy as np
import astropy.io.fits as pyfits
from sedpy import observate

def make_mock(objname='', spec_snr=None, phot_snr=20, add_noise=False,
              noisefactor=1.0, filters=[], calibrated=True,
              mockparams={}, model=None, sps=None, **kwargs):

    mock = {}
    
    # Read in a calibrated GGC object to use its SNR, wavelength vector, and
    # calibration vector
    ggc = ggc_spec(objname, 'a', '1', fluxtype=None, **kwargs)
    if calibrated:
        calvector = 1.0
    else:
        calvector = ggc['calibration']
    if spec_snr is None:
        spec_snr = ggc['spectrum'] / (ggc['unc'] * noisefactor)
        mock['noisefator'] = noisefactor

    # Set up the output mock observations dictionary
    mock['wavelength'] = ggc['wavelength']
    mock['filters'] = observate.load_filters(filters)

    # Build and store the intrinsic SED
    model.params.update(mockparams)
    spec, phot, _ = model.sed(model.theta, mock, sps=sps)
    mock['intrinsic_true_spectrum'] = spec.copy()
    mock['intrinsic_true_maggies'] = phot.copy()
    mock['mock_params'] = model.params.copy()

    # Apply calibration
    spec *= calvector
    mock['calibration'] = np.array(calvector).copy()
    
    # Add noise
    noise_sigma = spec / spec_snr
    pnoise_sigma = phot / phot_snr
    noise = np.random.normal(0, 1, len(spec)) * noise_sigma
    pnoise = np.random.normal(0, 1, len(phot)) * pnoise_sigma
    mock['add_noise'] = add_noise
    if add_noise:
        spec += noise.copy()
        mock['added_noise'] = noise.copy()
        phot += pnoise.copy()
        mock['added_phot_noise'] = pnoise.copy()
    else:
        mock['added_noise'] = noise*0.0
        mock['added_phot_noise'] = pnoise*0.0
        
    mock['spectrum'] = spec
    mock['maggies'] = phot
    mock['unc'] = noise_sigma 
    mock['maggies_unc'] = pnoise_sigma

    return mock

def ggc_spec(objname, exp='a', ext='1', fluxtype=None,
             datadir='./spectra/', **extras):
    """
    :param fluxtype: (default: None)
        Flag describing the the type of flux calibration:
        
        * None: Use the flux calibrated spectrum
        * 0: Use the variance-weighted, CR-cleaned, sky-subtracted,
          spectrum.
        * 1: Use the unweighted, uncleaned, sky-subtracted, spectrum.
          Don't ever use this.
    """
    name = objname.upper().strip().replace(' ','')
    obs = {'object_name':name}
    sfile = '{0}_{1}_{2}.fits'.format(name.upper(), exp, ext)
    sfile = os.path.join(datadir, sfile)
    auxfile = sfile.replace('.fits','.aux.fits')
    if not os.path.exists(sfile):
        raise ValueError('{0} does not exist!'.format(sfile))
    
    shdr, sdata = pyfits.getheader(sfile), pyfits.getdata(sfile)
    ahdr, adata = pyfits.getheader(auxfile), pyfits.getdata(auxfile)
    crpix = (shdr['CRPIX1'] -1) #convert from FITS to numpy indexing
    assert sdata.shape[0] == adata.shape[1]
    assert sdata.ndim == 1
    
    try:
        cd = shdr['CDELT1']
    except (KeyError):
        cd = shdr['CD1_1']

    obs['wavelength'] = (np.arange(len(sdata)) - crpix) * cd + shdr['CRVAL1']
    if fluxtype is None:
        obs['spectrum'] = sdata.copy()
        obs['spec_units'] = shdr['BUNIT']
        obs['calibration'] = adata[0,:] / obs['spectrum']
    else:
        obs['spectrum'] = adata[fluxtype, :].copy()
        obs['spec_units'] = 'Counts'
        obs['calibration'] = 1.0
    obs['unc'] = obs['spectrum'] / adata[3,:]
    obs['sky'] = adata[2,:] * (obs['spectrum'] / adata[0,:])
    return obs
