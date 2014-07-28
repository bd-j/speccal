import matplotlib.pyplot as pl
import numpy as np
import sps_basis
import diagnostics as dg

sps = sps_basis.StellarPopBasis()

showpars = ['mass', 'tage', 'zmet', 'dust2', 'zred', 'sigma_smooth', 'poly_coeffs1', 'spec_norm']

resdir = '/Users/bjohnson/Projects/m31spectra/results/'
fileraw = resdir + 'b192-g242v_020_1400723651.sampler01_mcmc'
filecal = resdir + 'b192-g242s_020_1400736872.sampler01_mcmc'
obj = 'B192.020'
cal = ['raw', 'cal']
filenames = [fileraw, filecal]
outnames = ['{0}.{1}'.format(obj, c) for c in cal]
color = ['red','blue']

start, thin, nspec = 250, 1, 10
ptiles = [16., 50., 84.]
qvalues = np.zeros( [len(filenames), len(showpars), len(ptiles)] )

ndim = len(showpars)
#figure cruft
nx = int(np.floor(np.sqrt(ndim)))
ny = int(np.ceil(ndim*1.0/nx))
figsize = (10, 6)
cfig, axes = pl.subplots(nx, ny, figsize = figsize)

for i, (filen, out) in enumerate(zip(filenames, outnames)):
    # Do the standard diagnostic plots
    _, results, model = dg.diagnostic_plots(filen, sps, showpars =showpars,nspec = nspec,
                                            start = start, thin = thin, outname = out)
    # Now get quantiles and plot distributions
    # This is repeating code in subtriangle and corner
    parnames = np.array(dg.theta_labels( results['model'].theta_desc ))
    ind_show = np.array([p in showpars for p in parnames], dtype = bool)
    if ind_show.sum() != ndim:
        raise KeyError('one of your showpars was not correct')
    
    flatchain = results['chain'][:, start::thin, ind_show]
    dims = flatchain.shape
    flatchain = flatchain.reshape(dims[0] * dims[1], dims[2])
    q = np.array(np.percentile(flatchain, ptiles, axis = 0))
    #print(q.shape, qvalues[i,:,:].shape)
    qvalues[i,:,:] = q.T
    #print(i, axes.shape)#, axes.flatten().shape)
    
    for j, par in enumerate(parnames[ind_show]):
     #   print(j, axes.flatten().shape, axes.shape)
        ax = axes.flatten()[j]
        ax.hist(flatchain[:,j], bins = 30,
                histtype ='stepfilled', alpha = 0.5,
                color = color[i], label = cal[i])
        ax.set_xlabel(par)
        if j == ndim-1:
            ax.legend(loc =0., fontsize ='small')
        ax.set_yticklabels([])
        xt = ax.get_xticks()
        #print(j, xt)
        #ran = xt[-1] - xt[0]
        #s, e =
        if (i == len(cal)-1) and (len(xt) >= 5):
            xt = np.linspace(xt[1], xt[-2], 3)
            ax.set_xticks(xt.copy())

cfig.subplots_adjust(hspace = 0.3, wspace =0.05)
cfig.savefig('{0}.chist.pdf'.format(obj))
cfig.show()
#pl.close(cfig)

pl.figure()
pl.clf()
norm = qvalues[:,:,1].mean(axis =0)
norm = qvalues[1,:,1]
pl.axhline(y = 1.0, linestyle = ':', color = 'k')
for i, out in enumerate(outnames):
    pl.errorbar(np.arange(ndim)+1 + (i-0.5)*0.1, qvalues[i,:,1]/norm,
                yerr = [(qvalues[i,:,1] - qvalues[i,:,0]) / norm,
                        (qvalues[i,:,2] - qvalues[i,:,1]) / norm],
                color = color[i], label = cal[i],
                linestyle = 'None', marker ='o')

pl.xticks(np.arange(ndim)+1, parnames[ind_show], fontsize =12)
pl.ylabel(r'$value / value_{{p50, cal}}$')
pl.legend(loc = 0, fontsize ='small')
pl.ylim(0.5, 1.5)
pl.savefig('{0}.cratio.pdf'.format(obj))
pl.show()
