import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

import diagnostics

import sps_basis
sps = sps_basis.StellarPopBasis(smooth_velocity=False)


def data_figure(results_list, layout = None, shaded = False, **kwargs):
    fig = pl.figure()

    #get the panel geometry if not given
    nobj = len(results_list)
    if layout is None:
        nx = int(np.floor(np.sqrt(nobj)))
        ny = int(np.ceil(nobj*1.0/nx))
        layout = [ny,nx]
    print(ny, nx)
    gs = gridspec.GridSpec(layout[0], layout[1])
    for i, res in enumerate(results_list):
        obs = res['obs']

        if shaded:
            ax = pl.Subplot(fig, gs[i])
            ax = one_data_figure_shaded(obs, ax, color = 'blue', facecolor = 'grey')
            fig.add_subplot(ax)

        else:
            fig, gs = one_data_figure_sep(obs, fig, subplot_spec = gs[i], **kwargs)
    #pl.tight_layout(fig)
    return fig
        
def one_data_figure_shaded(obs, axobj, color='Blue', facecolor='Blue',
                           **kwargs):
    """
    Plot the spectrum on one panel with shading to show the uncertainty.
    """
    
    x, y, e = obs['wavelength'], obs['spectrum'], obs['unc']
    axobj.fill_between(x, y-e, y+e, facecolor='grey', alpha=0.3)
    axobj.plot(x, y, color = color, linewidth = 0.5,**kwargs)

    return axobj

def one_data_figure_sep(obs, fig, subplot_spec=None, **kwargs):
    """
    Plot the spectrum and uncertainties in different panels
    """
    if subplot_spec is None:
        gs = gridspec.GridSpec(2,1,height_ratios = [3,1], hspace=0)
    else:
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, hspace=0,
                                              subplot_spec=subplot_spec,
                                              height_ratios = [3,1])
    
    
    spec = pl.Subplot(fig, gs[0,0])
    spec.plot(obs['wavelength'], obs['spectrum'], **kwargs)
    spec.set_ylabel(r'$f_\lambda \times \, C$')
    pl.setp(spec.get_xticklabels(), visible = False)
    fig.add_subplot(spec)
    unc = pl.Subplot(fig, gs[1,0])
    unc.plot(obs['wavelength'], obs['unc'], **kwargs)
    unc.set_ylabel(r'$\sigma f_\lambda$')
    unc.set_xlabel(r'$\lambda (\AA)$')
    fig.add_subplot(unc)
    return fig, gs
    
def spec_figure(results, alpha=0.3, samples=[-1],
                subplot_spec=None, **kwargs):
    """
    plot stars+dust+neb, then the calibration vector, then the GP
    predicitions, then data and full model, then residuals
    """
    
    obs = results['obs']
    fig = pl.figure(figsize = (10,6))
    gs = gridspec.GridSpec(3, 2, height_ratios=[3,1,1], hspace=0, wspace=0.05)
    
    # Axis variables
    ylabels = [r'$\mu$',r'$f(\alpha)$', r'$\tilde{\Delta}$ (GP)', #r'$\delta[f(\alpha)\mu]$ (GP)',
               'L','o/m', r'(o-m)/$\sigma$']
    color = ['blue', 'green', 'red', 'cyan', 'orange','orange']
    row = [0,1,2,0,1,2]
    col = [0,0,0,1,1,1]
    show_tick = [False, False, True, False, False, True]
    
    # generate axes
    axes = [pl.Subplot(fig, gs[c,r]) for r,c in zip(col, row)]
    # suppress ticks on shared axes
    [pl.setp(ax.get_xticklabels(), visible = sh) for ax,sh in zip(axes, show_tick)]
    # y-labels
    [ax.set_ylabel(label) for ax, label in zip(axes, ylabels)]
    [pl.setp(ax.yaxis, label_position='right') for ax in axes[3:]]
    [ax.yaxis.set_ticks_position('right') for ax in axes[3:]]
    [ax.yaxis.set_ticks_position('both') for ax in axes[3:]]

    # plot the data
    wave, ospec, mask = obs['wavelength'], obs['spectrum'], obs['mask']
    axes[3].plot(wave[mask], ospec[mask], color='magenta', label='Obs', **kwargs)

    # plot posterior draws
    nw, ns, nd = results['chain'].shape
    flatchain = results['chain'].reshape(nw*ns, nd)
    for s in samples:
        theta = flatchain[s,:]
        comps = diagnostics.model_components(theta, results, results['obs'],
                                             sps, photflag=0)
        spec, gp, cal, mask, stars = comps
        full_cal = cal[mask] + gp/stars[mask]
        mwave, ospec, mod = obs['wavelength'][mask],  obs['spectrum'][mask], (spec + gp)
        vecs = [stars[mask], cal[mask], gp,
                spec+gp, mod/ospec, (ospec-mod) / obs['unc'][mask]]
        [ax.plot(mwave, v, color=c, alpha=alpha, **kwargs) for ax,v,c in zip(axes, vecs, color)]

    [a.axhline( int(i==0), linestyle=':', color='black') for i,a in enumerate(axes[-2:])]
    # add axes to the figure
    [fig.add_subplot(ax) for ax in axes]
    return fig

def phot_figure(results, alpha=0.3, samples = [-1], **kwargs):
    """
    Plot the photometry for the model and data (with error bars). Then
    plot residuals
    """
    obs, _, marker = diagnostics.obsdict(results, 1)
    fig = pl.figure()
    gs = gridspec.GridSpec(2,1, height_ratios=[3,1])
    gs.update(hspace=0)
    phot, res = pl.Subplot(fig, gs[0]), pl.Subplot(fig, gs[1])
    res.set_ylabel( r'(obs-model)/$\sigma$')
    phot.set_ylabel('maggies')
    
    # plot posterior draws
    nw, ns, nd = results['chain'].shape
    flatchain = results['chain'].reshape(nw*ns, nd)
    label = 'modeled'
    for s in samples:
        theta = flatchain[s,:]
        comps = diagnostics.model_components(theta, results, obs,
                                             sps, photflag=1)
        spec, gp, cal, mask, stars = comps
        wave, ospec, ounc = obs['wavelength'][mask], obs['spectrum'][mask], obs['unc'][mask]
        phot.plot(wave, spec, label = label,
                  alpha=alpha, color='cyan', marker=marker, **kwargs)
        res.plot(wave, (ospec - spec) / ounc,
                 alpha=alpha, color='cyan', marker=marker, **kwargs)
        label = None
    phot.errorbar(wave, ospec, yerr=ounc,
                  color='magenta')
    phot.plot(wave, ospec, label = 'observed',
              color='magenta', marker=marker, **kwargs)
    phot.legend(loc=0)
    res.axhline(0, linestyle=':', color='black')
    
    fig.add_subplot(phot)
    fig.add_subplot(res)
    
    return fig

def calibration_figure(results_list):
    #plot the calibration vector and the 
    pass

    
def corner_plot(results, showpars=None):
    #just wrap subtriangle
    pass

def comparison_figure(results_list):
    pass

if __name__ == '__main__':
    figext = '.png'
    nsample = 5
    samples = np.random.uniform(0, 1, size=nsample)
    
    results = []
    rdir = '/Users/bjohnson/Projects/cetus/results/'
    res = [rdir+'b192-g242.020.cal_1405648278.sampler01',
           rdir+'b192-g242.020.nocal_1405677518.sampler01']
    name = ['B192 cal.', 'B192 no cal.']

    for i,r in enumerate(res):
        sf, mf = r+'_mcmc', r+'_model'
        result, pr, model = diagnostics.read_pickles(sf, model_file = mf)
        of = result['run_params']['outfile'].split('/')[-1]
        
        ns = result['chain'].shape[0] * result['chain'].shape[1]
        sample = [int(s * ns) for s in samples]
    
        sfig = spec_figure(result, samples=sample,
                        linewidth = 0.5)
        sfig.suptitle(name[i])
        sfig.savefig('sfig.'+ of + figext)
        pfig = phot_figure(result, samples=sample)
        pfig.suptitle(name[i])
        pfig.savefig('pfig.' + of + figext)

        results += [result]
        pl.close(sfig)
        pl.close(pfig)
        
    dfig = data_figure(results,
                       color='b', linewidth=0.5)
    [dfig.axes[i].set_title(name[i]) for i in range(len(name))]

    

    #dfig.savefig(of+'.data.pdf')
