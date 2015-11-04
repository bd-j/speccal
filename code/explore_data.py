import os, glob
import numpy as np
import matplotlib.pyplot as pl
import ggcdata

names = glob.glob('../data/ggclib/spectra/*.fits')
names = [os.path.basename(n).split('_')[0] for n in names]
unique_names = np.unique(np.array(names))
names = unique_names


def objdata(n):
    spec = ggcdata.ggc_spec(name, datadir='../data/ggclib/spectra')
    phot = ggcdata.ggc_phot(name, datadir='../data/ggclib/photometry')
    rc, rt, rh = ggcdata.gc_structural_params(name, datadir='../data/ggclib/photometry')
    nphot = np.isfinite(phot['maggies'][2:6]).sum()
    spec_snr = (spec['spectrum']/spec['unc']).max()
    spec_mags = np.array([f.ab_mag(spec['wavelength'], spec['spectrum'])
                          for f in phot['filters']])

    return n, nphot, spec_snr, -2.5 * np.log10(phot['maggies'][2]), rh


if __name__ == "__main__":
    ngood = 0
    out = open('ggc_summary.dat','w')
    out.write('Name   phot_flag   SNR_spec_max  g_AB  r_h \n')

    for name in unique_names:
        values = objdata(name)
        print(len(values), values)
        out.write('{0} {1} {2:4.0f}  {3:4.2f}  {4:6.1f}\n'.format(*values))
        #ngood += int(status)

    out.close()

#status = np.array([(np.isfinite(p['maggies'])).sum() == 7 for p in phot])
#status = status & (names != 'NGC6553') & (names != 'NGC2808') #& (names != 'NGC7089')
#gobs = np.array([p['maggies'][2] for p in phot])
#ind = gobs.tolist().index(gobs[status].max())
#pl.figure()    
#pl.plot([s[2] for s in struct], -2.5*np.log10(np.array(gobs)) - np.array(gmag), 'o')
        
#pl.figure()

#pl.plot(spec[ind]['wavelength'], spec[ind]['wavelength']*spec[ind]['spectrum']/fratio[ind])
#pconv = 3631e-23 * 2.998e18/w
#pl.plot(w, phot[ind]['maggies'] * pconv)
#pl.title(names[ind])

