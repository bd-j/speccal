import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits as pyfits

from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u


def read_reddening():
    abands = ['u', 'g', 'r', 'i', 'z']
    dt = np.dtype([('name', 'S10')] + zip(['A_{}'.format(b) for b in abands], len(abands) * [np.float]))
    with open('Schlafly_Finkbeiner_NED_reddening.txt', 'r') as f:
        lines = f.readlines()
    count = 0
    name, att, comment = [], [], []
    for l in lines:
        if l[0] == '#':
            continue
        cols = l.split()
        if len(cols) == 0:
            continue
        name.append(cols[0])
        tline = tuple([cols[0]] + [float(a) for a in cols[1:6]])
        att.append(tline)
        comment.append(' '.join(cols[6:]))
    return np.array(name), np.array(att, dtype=dt), comment

_, vred, _ = read_reddening()
sdss = pyfits.getdata('sdss.fits')
ctio = pyfits.getdata('ctio.fits')
harris = pyfits.getdata('harris97.fits')
hnames = [n.upper().strip().replace(' ','') for n in harris['ID']]
cnames = [n.upper().strip().replace(' ','') for n in ctio['Name']]
snames = [n.upper().strip().replace(' ','') for n in sdss['Name']]
vnames = [n.upper().strip().replace(' ','') for n in vred['name']]

hcols = {'HBR': 'HBR', 'met': '__Fe_H_', 'ebv': 'E_B-V_',
         'vel': 'Vr', 'dist': 'Rsun', #km/s, kpc
         'Rc': 'Rc', 'Rhalf': 'Rh', 'conc': 'c', # arcmin
         'RA': 'RA2000', 'Dec': 'De2000'}

bands = ['g', 'r', 'i', 'z']
abands = ['A_'+b for b in bands]
def make_cat():
    cols = bands + abands + hcols.keys()
    mydt = [('name', 'S10'), ('tel', 'S4')] + zip(cols, len(cols) * [np.float])
    mycat = np.zeros(len(harris), dtype=np.dtype(mydt))

    for i, n in enumerate(hnames):
        print(n)
        mycat[i]['name'] = n

        for m, h in hcols.items():
            mycat[i][m] = harris[h][i]

        if n in vnames:
            vcat = vred[vnames.index(n)]
            for b in bands:
                mycat[i]['A_'+b] = vcat['A_'+b]

        if n in cnames:
            cat = ctio[cnames.index(n)]
            tail = '_mag'
            tel = 'ctio'
        elif n in snames:
            cat = sdss[snames.index(n)]
            tail = 'mag'
            tel = 'sdss'
        else:
            continue

        mycat[i]['tel'] = tel
        for b in bands:
            mycat[i][b] = cat[b+tail]
            
            
    return mycat

if __name__ == "__main__":
    mycat = make_cat()
    good = (mycat['g'] > 1.0) & np.isfinite(mycat['g'])
    mycat = mycat[good]


    fig, axes = pl.subplots(1, 2)
    eax = axes[0]
    zax = axes[1]
    eax.plot( mycat['ebv'], mycat['g'] - mycat['z'], 'o')
    eax.set_ylim(-0.5, 2.0)
    zax.plot( mycat['met'], mycat['g'] - mycat['z'], 'o')
    zax.set_ylim(-0.5, 2.0)
    fig.show()
