import numpy as np
import matplotlib.pyplot as pl
import sys, os, glob
import astropy.io.fits as pyfits

def rfits(filename, fluxconv =1):
    dat = pyfits.getdata(filename)
    hdr = pyfits.getheader(filename)
    crpix = (hdr['CRPIX1'] -1) #convert from FITS to numpy indexing
    try:
        cd = hdr['CDELT1']
    except (KeyError):
        cd = hdr['CD1_1']
    wave = (np.arange(dat.shape[1]) - crpix) * cd + hdr['CRVAL1']
    spec = dat[0,:] * fluxconv
    unc = np.sqrt(dat[1,:]) * fluxconv
    return wave, spec, unc

specdir = os.getenv('HOME')+'/Projects/m31spectra/data/mmt/nocal/'
files = glob.glob(specdir + '*.s.fits')

expnum = np.array([int(os.path.basename(f).split('.')[0]) for f in files])
objname = np.array([os.path.basename(f).split('.')[1] for f in files])

uobj, uind = np.unique(objname, return_inverse =True)
uobj = uobj[:-1] #drop M020
nobj = len(uobj)

calcolor = ['blue', 'cyan', 'magenta', 'red']
#rawcolor = ['red', , 'orange', 'grey']
caltype = ['s','v']
calname = ['calib.', 'raw']
labpos = ['left','right']
ylab = [r'F$_\lambda \times $ C', 'counts/s/pix']
fig, axes = pl.subplots(6,2, figsize= (11, 13))
axes =axes.flatten()
iax = 0


for i,obj in enumerate(uobj):
    
    exps = expnum[objname == obj]
    
    maxf = 0
    for j, ct in enumerate(caltype):
        ax =axes[iax]
        for k,e in enumerate(exps):
            #print(obj,e,ct)
            cw, cs, cu = rfits('{0}{1:03.0f}.{2}.{3}.fits'.format(specdir,e,obj, ct))
            #rw, rs, ru = rfits('{0}{1:03.0f}.{2}.v.fits'.format(specdir,e,obj))
            ax.plot(cw,cs, label = 'Obs. # {0}'.format(e), color = calcolor[k])
            #ax.plot(rw,rs, label = 'Exp. # {0}, raw'.format(e), color = rawcolor[j])
            maxf = max([maxf, cs.max()])
            
        ax.yaxis.set_label_position(labpos[j])
        ax.set_ylabel(ylab[j])
        if j == 1:
            ax.set_yticklabels([])
            ax.legend(prop = {'size':6})
        ax.set_ylim(0,maxf*1.05)
        ax.set_xlim(3500,9500)
        ax.annotate('{0} {1}'.format(obj, calname[j]), xy = (5500,maxf*0.9))
        if i < (nobj-1):
                ax.set_xticklabels([])
        else:
            ax.set_xlabel(r'$\AA$')
        iax += 1
fig.subplots_adjust(hspace =0, wspace =0, left =0.10, right =0.9, bottom = 0.1, top = 0.95)
pl.savefig('mmt_cluster_data.pdf')
#pl.show()
pl.close()
