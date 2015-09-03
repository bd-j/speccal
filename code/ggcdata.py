import os
import numpy as np
try:
    import astropy.io.fits as pyfits
except(ImportError):
    import pyfits
try:
    from sedpy import observate
except:
    pass

ggcdir = os.path.join(os.environ['PROJECTS'], 'speccal',
                      'data', 'ggclib')


def ggc_spec(objname, exp='a', ext='1', fluxtype=None,
             datadir=ggcdir+'/spectra/', **extras):
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


def spec_lsf(wave, sigma_smooth=1.31,
             sigma_library=1.08, lsf_ref_wave=4838,
             **extras):
    """Line spread function of the GGC spectroscopy from Schiavon.

    :param wave:
        Wavelengths in angstroms
        
    :param sigma_smooth:
        Dispersion at ref_wave, in AA.
        
    :param sigma_lib:
        Resolution of the spectral library.  It will be subtracted in
        quadrature from the desired total resolution, so when the
        library is convolved with the result the desired LSF is
        acheived

    :param lsf_ref_wave:
        The wavelength in AA at which sigma_smooth is defined.
        
    :returns disp:
        The dispersion (or fwhm if sigma_smooth was gien as such) at
        each wavelength given by wave.
    """
    coeffs = np.array([15.290, -6.079e-3, 9.472e-7, -4.395e-11])
    powers = np.arange(len(coeffs))
    fwhm = np.dot(coeffs, wave[None, :] ** powers[:, None])
    fwhm_ref = (coeffs * lsf_ref_wave**powers).sum()
    target = (fwhm/fwhm_ref) * sigma_smooth
    delta_sigma_sq = np.clip(target**2 - sigma_library**2, 0, np.inf)
    return np.sqrt(delta_sigma_sq)

def broaden_ggcspec(wave, flux, minsig=1e-6):
    """Convolve the observed spectrum to have have constant wavelength
    resolution across the spectral range.  This is not kosher, but
    makes run times tractable.
    """
    sigma_data = spec_lsf(wave, sigma_library=0.0, sigma_smooth=1.31)
    delta_sigma = np.sqrt(sigma_data.max()**2 - sigma_data**2)
    delta_sigma[delta_sigma <= 0.0] = minsig
    def lsf(x):
        return delta_sigma
    newflux = observate.lsf_broaden(wave, flux, lsf=lsf)
    return newflux

def ggc_phot(objname, datadir=ggcdir+'/photometry/', **extras):
    
    name = objname.upper().strip().replace(' ','')
    obs = {'object_name':name}
    
    # Optical
    bands = ['g','r','i','z']
    obs['filternames'] = ['sdss_{0}0'.format(b) for b in bands]
    maggies, mags_unc, tel = optical_maggies(name, datadir=datadir, bands=bands)
    obs['maggies'] = maggies
    obs['maggies_unc'] =  mags_unc * maggies / 1.086
    obs['maggies_unc'] = np.sqrt(obs['maggies_unc']**2 + (0.05*maggies)**2)
    obs['optical_telescope'] = tel
    
    # NIR.  These are in Vega!!!
    obs['filternames'] += ['twomass_J', 'twomass_H', 'twomass_Ks']
    maggies, mags_unc = twomass_maggies(name, datadir=datadir)
    obs['maggies'] = np.concatenate([obs['maggies'], maggies])
    obs['maggies_unc'] = np.concatenate([obs['maggies_unc'], mags_unc/1.086*maggies])


    #Galex
    obs['filternames'] =  ['galex_FUV', 'galex_NUV'] + obs['filternames']
    maggies, mags_unc = np.zeros(2)+ np.nan, np.zeros(2)+ np.nan
    obs['maggies'] = np.concatenate([maggies, obs['maggies']])
    obs['maggies_unc'] = np.concatenate([mags_unc/1.086*maggies, obs['maggies_unc']])
    
    try:
        obs['filters'] = observate.load_filters(obs['filternames'])
        #convert 2mass mags to AB
        ab2v =  np.array([p.ab_to_vega for p in obs['filters']])
        obs['maggies'][-3:] /= 10**(-0.4 * ab2v[-3:])
    except(NameError):
        print('Warning:  Filter objects not loaded, 2MASS magnitudes are still Vega!!!!')
    
    return obs
    
def optical_maggies(name,  bands=['g','r','i','z'],
                    datadir=ggcdir+'/photometry/', **extras):
    """
    :returns maggies:
        g,r,i,z GC maggies within the half-light radius, on the AB
        system.
    :returns mags_unc:
        g,r,i,z *magnitude* uncertainties
    """
    tel = None
    # make more DRY
    try:
        # ctio catalog
        ctio = pyfits.getdata(os.path.join(datadir,'ctio.fits'))
        cnames = [n.upper().strip().replace(' ','') for n in ctio['Name']]
        opt = ctio[cnames.index(name)]
        mags = np.array([opt[b+'_mag'] for b in bands]).flatten()
        mags_unc = np.array([opt['e_{0}_mag'.format(b)] for b in bands]).flatten()
        tel = 'ctio'
    except(ValueError):
        pass
    try:
        # sdss catalog
        sdss = pyfits.getdata(os.path.join(datadir,'sdss.fits'))
        cnames = [n.upper().strip().replace(' ','') for n in sdss['Name']]
        opt = sdss[cnames.index(name)]
        mags = np.array([opt[b+'mag'] for b in bands]).flatten()
        mags_unc = np.array([opt['e_{0}mag'.format(b)] for b in bands]).flatten()
        tel = 'sdss'
    except(ValueError):
        pass
        
    if tel is None:
        mags = np.zeros(4) + np.nan
        mags_unc = np.zeros(4) + np.nan
        print('Warning: no optical photometry found for {0}'.format(name))
    return 10**(-0.4*mags), mags_unc, tel

def twomass_maggies(name, datadir=ggcdir+'/photometry/', **extras):
    """
    :returns maggies:
        J,H,K integrated GC maggies, on the Vega system.
    :returns mags_unc:
        J,H,K *magnitude* uncertainties
    """
    r_c, r_t, r_h = gc_structural_params(name, datadir=datadir)
    tmass = pyfits.getdata(os.path.join(datadir,'twomass.fits'))
    tnames = [n.upper().strip().replace(' ','') for n in tmass['Name']]
    try:
        nir = tmass[tnames.index(name)]
    except(ValueError):
        print('Warning: {0} not found in 2mass catalog, setting NIR mags to NaN'.format(name))
        return np.zeros(3)+ np.nan, np.zeros(3)+ np.nan
    integral = king(r_h, r_c=nir['Jrc'], r_t=r_t, f_0=1.0)[1]
    #print(dn0[0]*integral)
    magzp = np.array([20.45, 20.90, 19.93])
    dn0 = np.array([nir[col] for col in ['Ja0', 'Ha0', 'Ka0']])
    dn_unc = np.array([nir['e_'+col] for col in ['Ja0', 'Ha0', 'Ka0']])
    return integral * dn0 * 10**(-0.4 * magzp), 1.086 * dn_unc/dn0
    
def gc_structural_params(name, datadir=''):
    """
    :returns r_c:
        core radius in ''
    :returns r_t
        tidal radisu in arcsec
    :returns r_h:
        half_light radius in arcsec
    """
    harris = pyfits.getdata(os.path.join(datadir,'harris97.fits'))
    hnames = [n.upper().strip().replace(' ','') for n in harris['ID']]
    st = harris[hnames.index(name)]
    return 60.0 * st['Rc'], 60.0 * st['Rc'] * 10**st['c'], 60.0 * st['Rh']

def king(r, r_c=1., r_t=30., f_0=1.):
    x, xt = (r/r_c)**2, (r/r_t)**2
    xtc = (r_t/r_c)**2
    t1, t2 = 1 + x, 1 + xt
    k = f_0 * (1 - 1/(1 + xtc)**(0.5))**(-2)
    value = k * ((1 + x)**(-0.5) - (1 + xtc)**(-0.5))**(2)
    integral =  ( np.log(1 + x) + x/(1 + xt) -
                  4 * ((1+x)**(0.5) - 1) / (1 + xt)**(0.5) )
    integral *= np.pi * r_c**2 * k
                                       
    return value, integral
    
def integrated_flux_king(outer, A, r_c=1, r_t=30):
    """
    :param A:
        surface brightness (mag/sq'')
    :param outer:
        Outer radius ('')
    :param r_c:
        Core radius ('')
    """
    def func(x, c, t, k):
        return x * king(x, r_c=c, r_t=t, f_0=k)[0]
    
    f0 = 10**(-0.4 * A)
    total_flux = 2*np.pi * quad(func, 0, outer, args=(r_c, r_t, f0))
    return -2.5*np.log10(total_flux)


def make_skymask(wave, flux, width=10, thresh=3., lsf=6., **extras):
    """Generate a mask based on peaks in the sky spectrum, found using
    median filtering.

    :param width:
        Width of the median filter, in pixels
        
    :param thresh:
        Detection threshold, in sigma

    :param lsf:
        Distance from detected peaks to include in the mask of that
        peak, in wavelength units.

    :returns mask:
        Boolean array with ``False`` for pixels to mask due to bright
        sky emission.
    """
    #median filter
    shifts = np.arange(width) - int(width/2.0)
    medspec = np.median(np.vstack([np.roll(flux, s) for s in shifts]), axis=0) 
    sigma = (flux-medspec).std()
    #find outliers
    peak_ind = abs(flux-medspec) > (thresh * sigma)
    peak_wave = wave[peak_ind]
    #grow outliers based on lsf
    dwave = abs(peak_wave[:, None] - wave[None, :])/lsf
    mask = (dwave < 1).sum(axis=0)
    return mask == 0


def ggc_mock(model, theta, sps, objname='', apply_cal=True, mask=True,
             add_noise=False, phot_snr=30, spec_snr=None, **extras):
    """Generate a mock spectrum
    """
    mock = {}
    fnames = ['galex_FUV', 'galex_NUV',
              'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
              'twomass_J', 'twomass_H', 'twomass_Ks']

    # load the calibrated spectrum
    cal = ggc_spec(objname, 'a', '1', fluxtype=None, **extras)
    if spec_snr is None:
        spec_snr = cal['spectrum']/cal['unc']
    mock['filters'] = observate.load_filters(fnames)
    mock['wavelength'] = cal['wavelength']
    mock['sky'] = cal['sky']
    s, p, x = model.sed(theta, mock, sps=sps)
    mock['intrinsic_true_spectrum'] = s.copy()
    mock['intrinsic_true_maggies'] = p.copy()
    mock['mock_params'] = model.params
    mock['mock_theta'] = model.theta.copy()

    mock['calibration'], noise = 1.0, 0.0
    mock['added_noise'] = None
    if apply_cal:
        s *= cal['calibration']
        mock['calibration'] = cal['calibration'].copy()
    noise_sigma = s/spec_snr
    pnoise_sigma = p/phot_snr
    if add_noise:
        noise = np.random.normal(0, 1, len(s)) * noise_sigma
        s += noise.copy()
        mock['added_noise'] = noise.copy()
        pnoise = np.random.normal(0, 1, len(p)) * pnoise_sigma
        p += pnoise.copy()
        mock['added_phot_noise'] = pnoise.copy()
        
    mock['spectrum'] = s
    mock['maggies'] = p
    # Should use s/spec_snr??? no, that does not actually give the
    # likelihood of the mock data given the model and the uncertainty
    mock['unc'] = noise_sigma 
    mock['maggies_unc'] = pnoise_sigma

    if mask:
        mock = ggc_mask(mock)
    
    return mock

def ggc_mask(obs, minwave=3602, maxwave=1e4, pad=10.0, **kwargs):
    """Generate a mask based on peaks in the sky spectrum, min and max
    wavelengths, and regions around ISM and sky lines.
    """
    from bsfh import elines
    wave = obs['wavelength']
    skymask = make_skymask(wave, obs['sky'], **kwargs)
    mask = (skymask & (wave > minwave) &
            (wave < maxwave))
        
    for line in elines.ism_lines + elines.sky_lines:
        inline  = ((wave > (elines.wavelength[line] - pad)) &
                   (wave < (elines.wavelength[line] + pad)))
        mask = mask & ~inline

    # mask weird remaining lines
    inline = (wave > 6050.0) & (wave < 6075.0)
    inline = inline | ((wave > 5036) & (wave < 5055.0))
    inline = inline | ((wave > 4536) & (wave < 4556.0))
    inline = inline | ((wave > 5878) & (wave < 5890.0))
    inline = inline | ((wave > 6212) & (wave < 6232.0))
    inline = inline | ((wave > 6266) & (wave < 6300.0))
    mask = mask & ~inline
    
    obs['mask'] = obs.get('mask', True) & mask
    return obs
