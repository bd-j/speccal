import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec
import triangle

import diagnostics

import sps_basis
sps = sps_basis.StellarPopBasis(smooth_velocity=False)

def data_figure(results_list, layout = None, shaded = False, **kwargs):
    fig = pl.figure(figsize=(10, 5))

    #get the panel geometry if not given
    nobj = len(results_list)
    if layout is None:
        nx = int(np.floor(np.sqrt(nobj)))
        ny = int(np.ceil(nobj*1.0/nx))
        layout = [nx,ny]
    print(ny, nx)
    gs = gridspec.GridSpec(layout[0], layout[1])
    for i, res in enumerate(results_list):
        x, y = i % nx, np.floor(i*1.0 / nx)
        #print(x,y)
        obs = res['obs']

        if shaded:
            ax = pl.Subplot(fig, gs[x,y])
            ax = one_data_figure_shaded(obs, ax, color = 'blue', facecolor = 'grey')
            fig.add_subplot(ax)
        else:
            fig, sub_gs = one_data_figure_sep(obs, fig, subplot_spec = gs[i], **kwargs)
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
                start=0, thin=1,
                subplot_spec=None, xlim=None, **kwargs):
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
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    for s in samples:
        theta = flatchain[s,:]
        comps = diagnostics.model_components(theta, results, results['obs'],
                                             sps, photflag=0)
        spec, gp, cal, mask, stars = comps
        full_cal = cal[mask] + gp/stars[mask]
        mwave, ospec, mod = obs['wavelength'][mask],  obs['spectrum'][mask], (spec + gp)
        vecs = [stars[mask], cal[mask], gp,
                spec+gp, ospec/mod, (ospec-mod) / obs['unc'][mask]]
        [ax.plot(mwave, v, color=c, alpha=alpha, **kwargs) for ax,v,c in zip(axes, vecs, color)]

    [a.axhline( int(i==0), linestyle=':', color='black') for i,a in enumerate(axes[-2:])]
    if xlim is not None:
        [a.set_xlim(xlim) for a in axes]
    # add axes to the figure
    [fig.add_subplot(ax) for ax in axes]
    return fig


def phot_figure(results, alpha=0.3, samples = [-1],
                start=0, thin=1,
                **kwargs):
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
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
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
    

def calibration_figure(cal_result, nocal_result, samples = [-1],
                       start=0, thin=1,
                       alpha = 0.3):
    #plot the calibration vector and the reconstruction of it
    # from the f(alpha) and GP, for cal and uncal.

    fig = pl.figure(figsize=(10, 5))
    
    gs = gridspec.GridSpec(1, 2)
    cal_plot, nocal_plot = pl.Subplot(fig, gs[0]), pl.Subplot(fig, gs[1])

    ra = zip([cal_result, nocal_result], [cal_plot, nocal_plot])
    for i, (result, ax) in enumerate(ra):
        # plot posterior draws
        flatchain = result['chain'][:,start::thin,:]
        flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                      flatchain.shape[2])

        label = 'modeled'
        for s in samples:
            theta = flatchain[s,:]
            comps = diagnostics.model_components(theta, result, result['obs'],
                                             sps, photflag=0)
            spec, gp, cal, mask, stars = comps
            full_cal = cal[mask] + gp/stars[mask]
            mwave = result['obs']['wavelength'][mask]
            ax.plot(mwave, 1/full_cal, label = label,
                          color = 'cyan', alpha = alpha)
            label = None


    cal_plot.axhline(1.0, color='magenta', label='Caldwell')
    cal_plot.legend(loc=0)
    caldwell_cal = cal_result['obs']['spectrum']/nocal_result['obs']['spectrum']
    nocal_plot.plot(cal_result['obs']['wavelength'][mask], caldwell_cal[mask]*0.75, color='magenta', label='Caldwell')
    
    nocal_plot.legend(loc=0)
    
    print(len(cal_result['obs']['spectrum']), len(nocal_result['obs']['spectrum']))

    fig.add_subplot(cal_plot)
    fig.add_subplot(nocal_plot)
    return fig
    
def corner_plot(results, showpars=None, start=0, thin=1):
    #just wrap subtriangle
    """
    Make a triangle plot of the (thinned, latter) samples of the posterior
    parameter space.  Optionally make the plot only for a supplied subset
    of the parameters.
    """

    # pull out the parameter names and flatten the thinned chains
    parnames = np.array(diagnostics.theta_labels(results['model'].theta_desc))
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    truths = results['initial_center']

    # restrict to parameters you want to show
    if showpars is not None:
        ind_show = np.array([p in showpars for p in parnames], dtype= bool)
        flatchain = flatchain[:,ind_show]
        truths = truths[ind_show]
        parnames= parnames[ind_show]
        
    fig = triangle.corner(flatchain, labels = parnames,
                          quantiles=[0.16, 0.5, 0.84], verbose=False,
                          truths = truths)

    return fig

def comparison_figure(results_list):
    pass

if __name__ == '__main__':
    figext = '.png'
    nsample = 5
    samples = np.random.uniform(0, 1, size=nsample)
    showpars_phys = ['mass', 'tage', 'zmet', 'dust2', 'sigma_smooth']
    showpars_cal = ['zmet', 'dust2', 'poly_coeffs1', 'poly_coeffs2', 'gp_length', 'gp_amplitude']

    
    results = []
    rdir = '/Users/bjohnson/Projects/cetus/results/'
    res = [rdir+'b192-g242.020.cal_1405648278.sampler01',
           rdir+'b192-g242.020.nocal_1405677518.sampler01']
    name = ['B192 cal.', 'B192 no cal.']
    
    for i,r in enumerate(res):
        sf, mf = r+'_mcmc', r+'_model'
        result, pr, model = diagnostics.read_pickles(sf, model_file = mf)
        of = result['run_params']['outfile'].split('/')[-1]
        of = of.replace('.','_')
        ns = result['chain'].shape[0] * result['chain'].shape[1]
        sample = [int(s * ns) for s in samples]
    
        sfig = spec_figure(result, samples=sample,
                        linewidth = 0.5, xlim = (3650, 7300))
        sfig.suptitle(name[i])
        sfig.savefig('sfig_'+ of + figext)
        pl.close(sfig)
        
        pfig = phot_figure(result, samples=sample)
        pfig.suptitle(name[i])
        pfig.savefig('pfig_' + of + figext)
        pl.close(pfig)

        tfig = corner_plot(result, showpars = showpars_phys, start=-250)
        tfig.suptitle(name[i])
        tfig.savefig('ptri_' + of + figext)
        pl.close(tfig)

        tfig = corner_plot(result, showpars = showpars_cal, start=-250)
        tfig.suptitle(name[i])
        tfig.savefig('ctri_' + of + figext)
        pl.close(tfig)

               
        results += [result]
        
    dfig = data_figure(results,
                       color='b', linewidth=0.5)
    [dfig.axes[i*2].set_title(name[i]) for i in range(len(name))]
    dfig.savefig('dfig_'+ '_'.join(of.split('_')[0:2]) + figext)

    cfig = calibration_figure(*results, samples =sample)
    [cfig.axes[i].set_title(name[i]) for i in range(len(name))]
    [cfig.axes[i].set_ylabel(r'$f_{{\lambda}}\,/\, f_{{\lambda, obs}}$') for i in range(len(name))]
    cfig.savefig('cfig_'+ '_'.join(of.split('_')[0:2]) + figext)


