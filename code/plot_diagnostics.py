import sys, pickle
import numpy as np
import matplotlib.pyplot as pl

import bsfh.read_results as bread
from bsfh.gp import ExpSquared as GaussianProcess
from bsfh import sps_basis
sps = sps_basis.StellarPopBasis()
gp = GaussianProcess(None, None)

from plotting import *

if __name__ == '__main__':

    if len(sys.argv) > 1:
        resfile=sys.argv[1]
    else:
        resfile = 'results/ggc_mock.u0.t1.0_z0.0_a0.5_1426268715_mcmc'        
    model_file = resfile.replace('_mcmc','_model')

    showpars_phys = ['mass', 'tage', 'zmet', 'dust2', 'zred', 'sigma_smooth']
    showpars_cal = ['poly_coeffs_1', 'poly_coeffs_2', 'spec_norm',
                    'gp_length', 'gp_amplitude', 'gp_jitter', 'gp_jitter_add',
                    'phot_jitter']

    showpars_pc = ['mass', 'tage', 'dust2', 'zmet',
                   'spec_norm', 'poly_coeffs_1', 'poly_coeffs_2']
    
    result, pr, model = bread.read_pickles(resfile, model_file=model_file)
    ns = result['chain'].shape[0] * result['chain'].shape[1]
    start = int(result['chain'].shape[1]*0.75)
    thin = 2
    
    #sys.exit()
    try:
        efig = bread.param_evol(result, showpars=showpars_phys+showpars_cal)
        efig.savefig(resfile.replace('_mcmc','.pevol.pdf'))
        pl.close(efig)
    except(IndexError):
        print('not enough parameters to show in pevol plot')
    try:
        tfig = bread.subtriangle(result, showpars=showpars_phys,
                                 start=start, thin=thin)
        tfig.savefig(resfile.replace('_mcmc','.pcorner.pdf'))
        pl.close(tfig)
    except(IndexError):
        print('not enough parameters to show in triangle plot')

    try:
        cfig = bread.subtriangle(result, showpars=showpars_cal,
                                 start=start, thin=thin)
        cfig.savefig(resfile.replace('_mcmc','.ccorner.pdf'))
        pl.close(cfig)
    except(IndexError):
        print('not enough parameters to show in triangle plot')

    pcfig = bread.subtriangle(result, showpars=showpars_pc,
                             start=start, thin=thin)
    pcfig.savefig(resfile.replace('_mcmc','.pccorner.pdf'))
    pl.close(pcfig)
