import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from plotting import *

pnames = ['mass','tage','zmet','dust2']
pmap = {'zmet': lambda x: 10**x
        }
plabel_map = {'mass':r'M$_*$',
              'tage': 'Age',
              'zmet': r'$Z/Z_\odot$',
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
    
    truns = ['results/ggc_mock_specphot_linear.u0.t0.3_z0.0_a0.5_5305122_1432333491_mcmc',
             'results/ggc_mock_specphot_linear.u0.t1.1_z0.0_a0.5_5305126_1432333542_mcmc',
             'results/ggc_mock_specphot_linear.u0.t3.0_z0.0_a0.5_5305127_1432333568_mcmc',
             'results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_5280432_1431898211_mcmc']
    zruns = []
    aruns = []
    runs = [truns]#, aruns, zruns]
    colors = pl.rcParams['axes.color_cycle'][:3] + ['red']
    labels = ['0.3', '1.1', '3.0', '9.0']
    pfig, paxes = pl.subplots(len(runs), 1)

    
    for ip, pruns in enumerate(runs):
        ax = np.array(paxes).flat[ip]
        results, models = read_results(pruns)
        pfig, ax = plot_delta_params(results, models, pnames,
                                     sfraction=sfraction, thin=thin,
                                     fax=(pfig, ax),
                                     pmap=pmap, plabel_map=plabel_map,
                                     colors=colors, labels=labels)

        ax.plot(np.arange(4), np.zeros(4), ':k')
        ax.legend(loc=0)
    pfig.show()
