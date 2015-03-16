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

def ggc_spec(objname, exp='a', ext='1',
             fluxtype=None, datadir='', **extras):
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

def ggc_phot(objname, datadir='', **extras):
    
    name = objname.upper().strip().replace(' ','')
    obs = {'object_name':name}
    
    # Optical
    bands = ['g','r','i','z']
    obs['filternames'] = ['sdss_{0}0'.format(b) for b in bands]
    maggies, mags_unc = optical_maggies(name, datadir=datadir, bands=bands)
    obs['maggies'] = maggies
    obs['maggies_unc'] =  mags_unc * maggies / 1.086
    
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
    
def optical_maggies(name, datadir='', bands=['g','r','i','z'], **extras):
    """
    :returns maggies:
        g,r,i,z GC maggies within the half-light radius, on the AB
        system.
    :returns mags_unc:
        g,r,i,z *magnitude* uncertainties
    """
    found = False
    try:
        # ctio catalog
        ctio = pyfits.getdata(os.path.join(datadir,'ctio.fits'))
        cnames = [n.upper().strip().replace(' ','') for n in ctio['Name']]
        opt = ctio[cnames.index(name)]
        mags = np.array([opt[b+'_mag'] for b in bands]).flatten()
        mags_unc = np.array([opt['e_{0}_mag'.format(b)] for b in bands]).flatten()
        found = True
    except(ValueError):
        found = False
    try:
        # sdss catalog
        sdss = pyfits.getdata(os.path.join(datadir,'sdss.fits'))
        cnames = [n.upper().strip().replace(' ','') for n in sdss['Name']]
        opt = sdss[cnames.index(name)]
        mags = np.array([opt[b+'mag'] for b in bands]).flatten()
        mags_unc = np.array([opt['e_{0}mag'.format(b)] for b in bands]).flatten()
        found = found or True
    except(ValueError):
        found = found or False
        
    if not found:
        mags = np.zeros(4) + np.nan
        mags_unc = np.zeros(4) + np.nan
        print('Warning: no optical photometry found for {0}'.format(name))
    return 10**(-0.4*mags), mags_unc

def twomass_maggies(name, datadir='', **extras):
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


def skymask(wave, flux, width=10, thresh=5., lsf=6.):
    """Generate a mask based on peaks in the sky spectrum.
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
    
def ggc_mock(model, theta, sps, objname='', apply_cal=True,
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
    
    s, p, x = model.sed(theta, mock, sps=sps)
    mock['intrinsic_true_spectrum'] = s.copy()
    mock['intrinsic_true_maggies'] = p.copy()
    mock['mock_params'] = model.params
    mock['mock_theta'] = model.theta.copy()

    mock['calibration'], noise = 1.0, 0.0
    mock['added_noise'] = add_noise
    if apply_cal:
        s *= cal['calibration']
        mock['calibration'] = cal['calibration'].copy()
    noise_sigma = s/spec_snr
    pnoise_sigma = p/phot_snr
    if add_noise:
        noise = np.random.normal(0, 1, len(s)) * noise_sigma
        s += noise
        pnoise = np.random.normal(0, 1, len(p)) * pnoise_sigma
        p += pnoise
        
    mock['spectrum'] = s
    mock['maggies'] = p
    # Should use s/spec_snr??? no, that does not actually give the
    # likelihood of the mock data given the model and the uncertainty
    mock['unc'] = noise_sigma 
    mock['maggies_unc'] = pnoise_sigma

    return mock
