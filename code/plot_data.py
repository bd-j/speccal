import numpy as np
import matplotlib.pyplot as pl
import ggcdata, os

#The speccal directory
sdir = os.path.join(os.environ['PROJECTS'], 'speccal')

##### OBSERVATIONAL DATA ######

if __name__ == "__main__":

    try:
        objname = sys.argv[1]
    except:
        objname = 'NGC1851'
    try:
        noisefactor = float(sys.argv[2])
    except:
        noisefactor = 10.0
        print('noisefactor = 1')

    clr = 'black'
    
    datadir = os.path.join(sdir, 'data/ggclib/spectra/')
    dat = ggcdata.ggc_spec(datadir+objname, 'a', '1', fluxtype=None)
    mdat = ggcdata.ggc_mask(dat, thresh=3, lsf=6)
    mask = mdat['mask']
    #skymask = ggcdata.make_skymask(dat['wavelength'],
    #                               dat['sky'], thresh=3, lsf=6)
    fig, ax = pl.subplots(3, 1, sharex=True, figsize=(5, 8))
    ax[0].plot(dat['wavelength'][mask],
               dat['spectrum'][mask]*1e14, label='Spectrum',
               color = clr)
    ax[0].set_ylabel(r'F$_\lambda$ ($10^{-14} erg\, s^{-1} cm^{-2} \AA^{-1}$)')
    ax[1].plot(dat['wavelength'][mask],
               dat['spectrum'][mask]/(dat['unc'][mask] * noisefactor), label='SNR',
               color='black')
    ax[1].set_ylabel('S/N')
    ax[2].plot(dat['wavelength'][mask], 1/dat['calibration'][mask]*1e18, label='cal',
               color=clr)
    ax[2].set_ylabel(r'Calibration (F$_\lambda/$counts)')
    ax[2].set_xlabel(r'$\lambda (\AA)$')
    fig.show()
    fig.savefig('../tex/figures/data.pdf')
