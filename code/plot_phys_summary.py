import sys, pickle, matplotlib
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec
from matplotlib import rcParams
#rcParams['xtick.direction'] = 'in'
#rcParams['ytick.direction'] = 'in'
        
from bsfh import read_results as bread
from plotting import *
matplotlib.colors.cnames.update(newcolors)

pnames = ['mass','tage','zmet','dust2']
pmap = {'mass': lambda x: np.log10(x),
        'tage': lambda x: np.log10(x)
        }
plabel_map = {'mass':r'log M$_*$',
              'tage': 'log Age',
              'zmet': r'log $Z/Z_\odot$',
              'dust2':r'A$_V$'
              }
pvlabel_map = {'mass':r'M$_*$',
              'tage': 'Age (Gyr)',
              'zmet': r'log $Z/Z_\odot$',
              'dust2':r'A$_V$'
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

if __name__ == "__main__":

    oldstyle=False
    
    truns = ['results/ggc_mock_specphot_linear.u0.t0.3_z0.0_a0.5_5305122_1432333491_mcmc',
             'results/ggc_mock_specphot_linear.u0.t1.1_z0.0_a0.5_5305126_1432333542_mcmc',
             'results/ggc_mock_specphot_linear.u0.t3.0_z0.0_a0.5_5305127_1432333568_mcmc',
             'results/ggc_mock_specphot_linear.u0.t6.0_z0.0_a0.5_5313538_1432675118_mcmc',
             'results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_5280432_1431898211_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.5_5313542_1432675345_mcmc'
             ]
    aruns = ['results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.0_5313557_1432675536_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.5_5313542_1432675345_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a1.0_5313614_1432675422_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a2.0_5391265_1434182123_mcmc'
             ]
    zruns = ['results/ggc_mock_specphot_linear.u0.t12.0_z-1.5_a0.5_5391259_1434126258_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z-1.0_a0.5_5391258_1434122880_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z-0.5_a0.5_5391256_1434100928_mcmc',
             'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.5_5313542_1432675345_mcmc'
             ]
        
    runs = [truns, aruns, zruns]
    vtype = ['Vary Age', 'Vary A$_V$', 'Vary Z']
    vary_params = ['tage', 'dust2', 'zmet']
    vary_lims = [(0.1, 12.1), (-0.1, 2.1), (-1.59, 0.2)]
    props = dict(boxstyle='round', facecolor='w', alpha=0.5)

    if not oldstyle:
        ytitle=r'$\Delta$'
        fig = pl.figure(figsize=(15,5))
        gs = [gridspec.GridSpec(2, 2) for i in range(len(runs))]
        leftpad, rightpad, margin = 0.06, 0.00, 0.05
        width = (1.0 - 2*margin)/len(gs)
        for ivp, (pruns, pvary, tgs) in enumerate(zip(runs, vary_params, gs)):
            results, models = read_results(pruns)
            left = margin + width * ivp + leftpad
            right =  margin + width * (ivp+1) - rightpad 
            print(left, right)
            tgs.update(left=left, right=right, hspace=0.0, wspace=0.0)
            paxes = [pl.Subplot(fig, spec) for spec in tgs]
            for ip, pname in enumerate(pnames):
                dax = paxes[ip]
                dfig, dax = deltafig_vspar(results, models, pname, pvary,
                                           sfraction=sfraction, thin=thin,
                                           fax=(None, dax), pmap=pmap,
                                           xlims=vary_lims[ivp],
                                           verbose=False, fractional=False)
                pretty_pname = plabel_map.get(pname, pname)
                pretty_vpname = pvlabel_map.get(pvary, pvary)
                dax.axhline(0.0, linestyle=':', color='k')
                dax.text(0.1, 0.95, ytitle+pretty_pname,
                         transform=dax.transAxes, fontsize=14,
                         verticalalignment='top', bbox=None)#props)
                dax.set_ylim(-0.1, 0.1)
                if (ip<=1):
                    dax.set_xticks([])
                else:
                    dax.set_xlabel(pretty_vpname)
                if (ip==0) or (ip==2):
                    dax.set_ylabel(ytitle, size=12)
                else:
                    dax.set_yticks([])
                fig.add_subplot(dax)
        fig.savefig('../tex/figures/vary_params.pdf')
        pl.close(fig)
                
    if oldstyle:
        labels = [['0.3 Gyr', '1.1  Gyr', '3 Gyr', '6 Gyr', '9 Gyr'],
                  ['0.0', '1.0']]

        colors = pl.rcParams['axes.color_cycle'][:3] + ['red', 'purple']
        props = dict(boxstyle='round', facecolor='w', alpha=0.5)
        pfig, paxes = pl.subplots(len(runs), 1)


        for ip, pruns in enumerate(runs):
            ax = np.array(paxes).flat[ip]
            results, models = read_results(pruns)
            pfig, ax = deltafig(results, models, pnames,
                                sfraction=sfraction, thin=thin,
                                fax=(pfig, ax),
                                pmap=pmap, plabel_map=plabel_map,
                                colors=colors, labels=labels[ip],
                                verbose=False)

            ax.plot(np.arange(4), np.zeros(4), ':k')
            ax.legend(loc=0)
            ax.text(0.1, 0.95, vtype[ip],
                transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

            ax.set_ylim(-0.2, 0.2)
        pfig.show()
