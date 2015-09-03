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
plabel_map = {'mass':r'log M$_*$',
              'tage': 'log Age',
              'zmet': r'log $Z/Z_\odot$',
              'dust2':r'$\tau_V$'
              }
pvlabel_map = {'mass':r'M$_*$',
              'tage': 'Age (Gyr)',
              'zmet': r'log $Z/Z_\odot$',
              'dust2':r'$\tau_V$'
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
    
    runs = ['results/ggc_mock_specphot_linear.u1.t12.0_z0.0_a0.5_1433337568_mcmc',
            'results/ggc_mock_specphot_linear.u2.t12.0_z0.0_a0.5_5362177_1433529210_mcmc',
            'results/ggc_mock_specphot_linear.u3.t12.0_z0.0_a0.5_5362178_1433535984_mcmc',
            'results/ggc_mock_specphot_linear.u4.t12.0_z0.0_a0.5_5363683_1433543728_mcmc',
            'results/ggc_mock_specphot_linear.u5.t12.0_z0.0_a0.5_5363685_1433556968_mcmc',
            'results/ggc_mock_specphot_linear.u6.t12.0_z0.0_a0.5_5387917_1433979023_mcmc',
            'results/ggc_mock_specphot_linear.u7.t12.0_z0.0_a0.5_5387918_1433980001_mcmc',
            'results/ggc_mock_specphot_linear.u8.t12.0_z0.0_a0.5_5387925_1433980286_mcmc',
            'results/ggc_mock_specphot_linear.u9.t12.0_z0.0_a0.5_5387929_1433982909_mcmc',
            'results/ggc_mock_specphot_linear.u10.t12.0_z0.0_a0.5_5388209_1433996870_mcmc'
            ]
    nreal = len(runs)
    results, models = read_results(runs)
    noiseless_run = 'results/ggc_mock_specphot_linear.u0.t12.0_z0.0_a0.5_5313542_1432675345_mcmc'
    bcolor = 'maroon'

    if as_hist:
        hfig, haxes = pl.subplots(2,2)
        for j, (res, mod) in enumerate(zip(results, models)):
            samples, pord = hist_samples(res, mod, pnames, thin=thin,
                                         start=sfraction)
            for i, (ax, name) in enumerate(zip(haxes.flatten(), pnames)):
                ax.hist(samples[:,i], bins=kwargs.get("bins", 50),
                        histtype="stepfilled",
                        color=kwargs.get("color", bcolor),
                        alpha=kwargs.get("alpha",0.3))
        
        res, pr, mod = bread.read_pickles(noiseless_run,
                                          model_file=noiseless_run.replace('_mcmc','_model'))
        samples, pord = hist_samples(res, mod, pnames, thin=thin,
                                     start=sfraction)
        truths = [res['obs']['mock_params'][k] for k in pord]
        for i, (ax, name) in enumerate(zip(haxes.flatten(), pnames)):
            ax.hist(samples[:,i], bins=kwargs.get("bins", 50),
                    histtype="stepfilled",
                    color=bcolor,
                    alpha=kwargs.get("alpha",0.3))
            ax.axvline(truths[i], color=kwargs.get('truth_color','k'),
                       label='Mock Truth')
            ax.set_xlabel(name, fontsize=8)
            ax.set_yticklabels([])
            pl.setp(ax.get_xticklabels(), fontsize=8)

        hfig.show()

    else:
        nfig, naxes = pl.subplots(2,2, figsize=(8, 5))
        res, pr, mod = bread.read_pickles(noiseless_run,
                                          model_file=noiseless_run.replace('_mcmc','_model'))
        samples, pord = hist_samples(res, mod, pnames, thin=thin,
                                     start=sfraction)
        ptiles = np.percentile(samples,[16, 50, 84], axis=0)
        ptiles = np.array([pmap.get(k, identity)(ptiles[:,i])
                            for i,k in enumerate(pord)])

        truths = [res['obs']['mock_params'][k] for k in pord]
        truths = np.array([pmap.get(k, identity)(truths[i])
                           for i,k in enumerate(pord)])

        delta = (ptiles - truths).T
        for i, (ax, name) in enumerate(zip(naxes.flatten(), pnames)):
            x = np.arange(nreal+2)-1
            ax.plot(x, x*0 + delta[1,i], color=bcolor)
            ax.fill_between(x, x*0 + delta[0,i], x*0 + delta[2,i],
                            color=bcolor, alpha=0.3)
            ax.axhline(0, linestyle=':', color=kwargs.get('truth_color','k'),
                       label='Mock Truth')
            ax.set_ylabel('$\Delta$' + plabel_map[name], fontsize=8)
            ax.set_xticklabels([])
            pl.setp(ax.get_yticklabels(), fontsize=8)
            ax.set_xlim(-0.5, nreal - 0.5)
            ax.set_ylim(-0.1, 0.1)
            if i > 1:
                ax.set_xlabel('Noise realization #', fontsize=8)
        for j, (res, mod) in enumerate(zip(results, models)):
            samples, pord = hist_samples(res, mod, pnames, thin=thin,
                                         start=sfraction)
            ptiles = np.percentile(samples,[16, 50, 84], axis=0)
            ptiles = np.array([pmap.get(k, identity)(ptiles[:,i])
                            for i,k in enumerate(pord)])
            delta = (ptiles - truths).T
            for i, (ax, name) in enumerate(zip(naxes.flatten(), pnames)):
                yerr = np.array([[delta[2,i] - delta[1,i]],
                                 [ delta[1,i] - delta[0,i]]])
                ax.errorbar([j], [delta[1,i]], yerr, capthick=2, color='black')
                ax.plot([j], [delta[1,i]], 'ok')
        nfig.show()
        nfig.savefig('../tex/figures/noise_realizations.pdf')
