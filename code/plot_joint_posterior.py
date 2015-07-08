import sys, pickle, matplotlib
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from plotting import *
nc = newcolors.copy()
#nc.pop('magenta')
matplotlib.colors.cnames.update(nc)

param_name_map = {'tage':'Age (Gyr)',
                  'mass': '$M_*$ $(M_\odot/10^{5})$',
                  'dust2':'$A_V$ (mag)',
                  'zmet': '$\log Z/Z_\odot$',
                  'sigma_smooth': '$\sigma_{{LSF}}$',
                  'zred': '${\it z}$',
                  }
pnmap = param_name_map

def joint_pdf(res, p1, p2, **kwargs):
    trace, pars = hist_samples(res, res['model'], [p1, p2], **kwargs)
    if p1 == 'mass':
        trace[:,1] /= 1e5
    if p2 == 'mass':
        trace[:,0] /= 1e5
    trace = trace.copy().T
    if pars[0] == p1:
        trace = trace[::-1, :]
    xbins, ybins, sigma = compute_sigma_level(trace[0], trace[1])
    return xbins, ybins, sigma.T


def compute_sigma_level(trace1, trace2, nbins=30):
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
    speconly = 'results/ggc_mock_speconly.u0.t9.0_z0.0_a0.5_1431313829_mcmc'
    specphot = 'results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_5280432_1431898211_mcmc'
    resfiles = [photonly, speconly, specphot]
    clr = ['red','blue', 'magenta']
    results = [bread.read_pickles(rfile, model_file=rfile.replace('mcmc','model'))[0]
               for rfile in resfiles]
    obsdat = results[0]['obs']
    showpars = np.array(['mass', 'tage', 'zmet', 'dust2'])
    parlims = np.array([[0.4, 4],
                        [None, None],
                        [-1.0, 0.19],
                        [0, 1.5]])
    npar = len(showpars)

    fig = pl.figure()
    
    gs = gridspec.GridSpec(npar, npar)
    for i, p1 in enumerate(showpars):
        dax = pl.subplot(gs[i,i])
        for n, res in enumerate(results):
            trace, p = hist_samples(res, res['model'], [p1], start=0.5)
            if p1 == 'mass':
                trace /= 1e5
            dax.hist(trace, bins = 30, color=clr[n], normed=True,
                     alpha=0.5, histtype='stepfilled')
            dax.tick_params(axis='both', which='major', labelsize=8)
        # Axis label foo
#        if i == 0:
#            dax.set_ylabel(pnmap.get(p1, p1))
#        else:
        dax.set_yticklabels('')
        if i == (npar-1):
            dax.set_xlabel(pnmap.get(p2, p2), fontsize=12)
        else:
            dax.set_xticklabels('')
        for j, p2 in enumerate(showpars[(i+1):]):
            k = j+i+1
            ax = pl.subplot(gs[k, i])

            for n, res in enumerate(results):
                pdf = joint_pdf(res, p2, p1, start=0.5)
                ax.contour(pdf[0], pdf[1], pdf[2], levels = [0.683, 0.955], colors=clr[n])
                
            # Axis label foo
            if i == 0:
                ax.set_ylabel(pnmap.get(p2, p2), fontsize=12)
            else:
                ax.set_yticklabels('')
            if k == (npar-1):
                ax.set_xlabel(pnmap.get(p1,p1), fontsize=12)
                dax.set_xlim(ax.get_xlim())
            else:
                ax.set_xticklabels('')

            # Axis range foo
            xcur = ax.get_xlim()
            ycur = ax.get_ylim()
            xlims = max([parlims[i,0], xcur[0]]), min([parlims[i,1], xcur[1]])
            ylims = max([parlims[j+i+1,0], ycur[0]]), min([parlims[j+i+1,1], ycur[1]])
            
            ax.set_xlim(*xlims)
            ax.set_ylim(*ylims)
            dax.set_xlim(*xlims)
            dax.tick_params(axis='both', which='major', labelsize=8)
            ax.tick_params(axis='both', which='major', labelsize=8)
            
            truths = np.squeeze(np.copy([obsdat['mock_params'][k] for k in [p1, p2]]))
            if p1 == 'mass':
                truths[0] /= 1e5
            if p2 == 'mass':
                truths[1] /= 1e5
                
            ax.plot(truths[0], truths[1], 'ok')
    fig.savefig('../tex/figures/combined_post.pdf')
    pl.show()
