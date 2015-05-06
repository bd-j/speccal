import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from plotting import *

def joint_pdf(res, p1, p2, **kwargs):
    trace, pars = hist_samples(res, res['model'], [p1, p2], **kwargs)
    trace = trace.copy().T
    xbins, ybins, sigma = compute_sigma_level(trace[0], trace[1])
    return xbins, ybins, sigma.T


def compute_sigma_level(trace1, trace2, nbins=20):
    """From a set of traces, bin by number of standard deviations"""
    L, xbins, ybins = np.histogram2d(trace1, trace2, nbins)
    L[L == 0] = 1E-16
    logL = np.log(L)

    shape = L.shape
    L = L.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(L)[::-1]
    i_unsort = np.argsort(i_sort)

    L_cumsum = L[i_sort].cumsum()
    L_cumsum /= L_cumsum[-1]
    
    xbins = 0.5 * (xbins[1:] + xbins[:-1])
    ybins = 0.5 * (ybins[1:] + ybins[:-1])

    return xbins, ybins, L_cumsum[i_unsort].reshape(shape)


if __name__ == "__main__":
    photonly = 'results/ggc_mock_photonly.c0.t9.0_z0.0_a0.5_1430274922_mcmc'
    speconly = 'results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5_1430808300_mcmc'
    specphot = 'results/ggc_mock_specphot.u0.t9.0_z0.0_a0.5_1430380402_mcmc'
    resfiles = [photonly, speconly, specphot]
    clr = ['red','blue', 'magenta']
    results = [bread.read_pickles(rfile, model_file=rfile.replace('mcmc','model'))[0]
               for rfile in resfiles]
    obsdat = results[0]['obs']
    showpars = np.array(['mass', 'tage', 'zmet', 'dust2'])
    npar = len(showpars)

    fig = pl.figure()
    
    gs = gridspec.GridSpec(npar, npar)
    for i, p1 in enumerate(showpars):
        dax = pl.subplot(gs[i,i])
        for n, res in enumerate(results):
            trace, p = hist_samples(res, res['model'], [p1], start=0.5)
            dax.hist(trace, bins = 50, color=clr[n], normed=True,
                     alpha=0.5, histtype='stepfilled')
            if i == 0:
                dax.set_ylabel(p1)
            else:
                dax.set_yticklabels('')
            if i == (npar-1):
                dax.set_xlabel(p2)
            else:
                dax.set_xticklabels('')
            
        for j, p2 in enumerate(showpars[(i+1):]):
            k = j+i+1
            ax = pl.subplot(gs[k, i])

            for n, res in enumerate(results):
                pdf = joint_pdf(res, p2, p1, start=0.5)
                ax.contour(pdf[0], pdf[1], pdf[2], levels = [0.683, 0.955], colors=clr[n])
            if i == 0:
                ax.set_ylabel(p2)
            else:
                ax.set_yticklabels('')
            if k == (npar-1):
                ax.set_xlabel(p1)
                dax.set_xlim(ax.get_xlim())
            else:
                ax.set_xticklabels('')
                
            truths = [obsdat['mock_params'][k] for k in [p1, p2]]
            ax.plot(truths[0], truths[1], 'ok')
    fig.savefig('../tex/figures/combined_post.pdf')
    pl.show()
