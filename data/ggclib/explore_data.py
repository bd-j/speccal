import os, glob
import numpy as np
import matplotlib.pyplot as pl
import ggc

names = glob.glob('spectra/*.fits')
names = [os.path.basename(n).split('_')[0] for n in names]

unique_names = np.unique(np.array(names))
spec, phot = [], []
struct, gmag, fratio = [], [], []

names = unique_names

pfig = pl.figure()
for name in unique_names:
    spec.append(ggc.ggc_spec(name, datadir='spectra'))
    phot.append(ggc.ggc_phot(name, datadir='photometry'))
    w = np.array([p.wave_effective for p in phot[-1]['filters']])
    f = phot[-1]['maggies']
    m = np.isfinite(f)
    pl.plot(w[m], f[m], '-o', label=name)
    tmp = phot[-1]['filters'][2].ab_mag(spec[-1]['wavelength'], spec[-1]['spectrum'])
    gmag.append(tmp)
    s = ggc.gc_structural_params(name, datadir='photometry')
    struct.append(s)
    int1 = ggc.king(s[0], r_c=s[0], r_t=s[1])[1]
    int2 = ggc.king(s[2], r_c=s[0], r_t=s[1])[1]
    fratio.append( int1/int2)
    print(name, int1/int2)


status = np.array([(np.isfinite(p['maggies'])).sum() == 7 for p in phot])
status = status & (names != 'NGC6553') & (names != 'NGC2808') #& (names != 'NGC7089')
gobs = np.array([p['maggies'][2] for p in phot])
ind = gobs.tolist().index(gobs[status].max())
#pl.figure()    
#pl.plot([s[2] for s in struct], -2.5*np.log10(np.array(gobs)) - np.array(gmag), 'o')
        
pl.figure()

pl.plot(spec[ind]['wavelength'], spec[ind]['wavelength']*spec[ind]['spectrum']/fratio[ind])
pconv = 3631e-23 * 2.998e18/w
pl.plot(w, phot[ind]['maggies'] * pconv)
pl.title(names[ind])

