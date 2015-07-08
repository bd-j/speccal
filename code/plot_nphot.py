import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from plotting import *

pnames = ['mass','tage','zmet','dust2']
pmap = {'mass': lambda x: np.log10(x),
        'tage': lambda x: np.log10(x)
        }
plabel_map = {'mass':'log M$_*$',
              'tage': 'log Age',
              'zmet': 'log $Z/Z_\odot$',
              'dust2':'$\tau_V$'
              }
pvlabel_map = {'mass':'M$_*$',
              'tage': 'Age (Gyr)',
              'zmet': 'log $Z/Z_\odot$',
              'dust2':'$\tau_V$'
              }

sfraction = 0.90
thin = 10

def read_results(runs):
    results, models = [], []
    for run in runs:
        res, pr, mod = bread.read_pickles(run, model_file=run.replace('_mcmc','_model'))
        results.append(res)
        models.append(mod)
    return results, models

def identity(x):
    return x


if __name__ == "__main__":

    kwargs = {}
    
    as_hist = False
    
    runs = ['results/ggc_mock_speconly.u0.t9.0_z0.0_a0.5_1431313829_mcmc', #g
            'results/ggc_mock_gi_spec_linear.u0.t12.0_z0.0_a0.5_5486144_1436374818_mcmc', #gi
            'results/ggc_mock_griz_spec_linear.u0.t12.0_z0.0_a0.5_5429625_1435022619_mcmc', #griz
            #'results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_5280432_1431898211_mcmc', #grizJHK
            'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.5_5313542_1432675345_mcmc', #grizJHK 
            ]
    runlabel = ['$g$','$gi$','$griz$', '$griz$\n$JHK_s$']
    pvary = [1, 2, 3, 4]
    nreal = len(runs)
    results, models = read_results(runs)
    ytitle=r'$\Delta$'
    
    pfig, paxes = pl.subplots(2,2, figsize=(8, 5))
    for ip, pname in enumerate(pnames):
        dax = paxes.flat[ip]
        dfig, dax = deltafig_vspar(results, models, pname, pvary,
                                   sfraction=sfraction, thin=thin,
                                   fax=(None, dax), pmap=pmap,
                                   xlims=(min(pvary)-0.5, max(pvary)+0.5),
                                   verbose=False, fractional=False)
        pretty_pname = plabel_map.get(pname, pname)
        #pretty_vpname = pvlabel_map.get(pvary, pvary)
        dax.axhline(0.0, linestyle=':', color='k')
        dax.text(0.6, 0.95, ytitle+pretty_pname,
                 transform=dax.transAxes, fontsize=14,
                 verticalalignment='top', bbox=None)#props)
        dax.set_ylim(-0.5, 0.5)
        if (ip<=1):
            dax.set_xticks([])
        else:
            dax.set_xlabel('N$_{phot}$')
            dax.set_xticks(pvary)
            dax.set_xticklabels(runlabel)
        if (ip==0) or (ip==2):
            dax.set_ylabel(ytitle, size=12)
        else:
            dax.set_yticks([])
        #fig.add_subplot(dax)
        
    pfig.savefig('../tex/figures/vary_nphot.pdf')
    #pl.close(fig)

