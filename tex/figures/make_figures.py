import sys
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec
import triangle

import bsfh.read_results as diagnostics
from bsfh.read_results import model_comp as mcomp

from bsfh import sps_basis
sps = sps_basis.StellarPopBasis(smooth_velocity=False)

pardict = {}
pardict['mass']=r'$m_*$ (M$_\odot$)'
pardict['tage']=r'$t$ (Gyr)'
pardict['zmet']=r'$\log Z/Z_\odot$'
pardict['dust2']=r'$A_V$'
pardict['sigma_smooth']=r'$\sigma \, (\AA)$'
pardict['zred']=r'$z$'
pardict['gp_jitter']=r'$s$'
pardict['gp_amplitude']=r'$a$'
pardict['gp_length']=r'$l \, (\AA)$'
pardict['poly_coeffs1']=r'$c_1$'
pardict['poly_coeffs2']=r'$c_2$'
pardict['spec_norm']=r'$c_0$'
pardict['emission_disp']= r'$\sigma_e \, (\AA)$'
pardict['phot_jitter']=r'$s_{{phot}}$'


def comp_samples(thetas, model, inlog=True, photflag=0):
    specvecs =[]
    obs, _, marker = obsdict(model.obs, photflag)
    wave, ospec, mask = obs['wavelength'], obs['spectrum'], obs['mask']
    mwave, mospec = wave[mask], ospec[mask]
    mounc = obs['unc'][mask]
    if inlog and (photflag == 0):
         mospec = np.exp(mospec)
         mounc *= mospec

    for theta in thetas:
        mu, cal, delta, mask, wave = mcomp(theta, model, sps,
                                           inlog=True, photflag=photflag)
        
        if inlog & (photflag == 0):
            full_cal = np.exp(cal + delta)
            mod = np.exp(mu + cal + delta)
            mu = np.exp(mu)
            cal = np.exp(cal)
            
        elif photflag == 0:
            full_cal = cal + delta/mu
            mod = (mu*cal + delta)

        else:
            mod = mu

        specvecs += [ [mu, cal, delta, mod,
                       mospec/mod, (mospec-mod) / mounc] ]
            
    return wave, mospec, mounc, specvecs


def spec_figure(results, alpha=0.3, samples=[-1],
                start=0, thin=1, inlog=True,
                subplot_spec=None, xlim=None, **kwargs):
    """
    plot stars+dust+neb, then the calibration vector, then the GP
    predicitions, then data and full model, then residuals
    """
    
    fig = pl.figure(figsize = (10,6))
    gs = gridspec.GridSpec(3, 2, height_ratios=[3,1,1], hspace=0, wspace=0.05)
    
    # Axis variables
    ylabels = [r'$\mu$',r'$e^{f(\alpha)}$', r'$\tilde{\Delta}$ (GP)', #r'$\delta[f(\alpha)\mu]$ (GP)',
               r'$\mu \, e^{f(\alpha) + \tilde{\Delta}}$','o/m', r'$\chi$']
    color = ['blue', 'green', 'red', 'magenta', 'orange','orange']
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

    # make posterior draws
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    thetas = [flatchain[s,:] for s in samples]
    thetas += [results['initial_center']]

    mwave, mospec, mounc, specvecs = comp_samples(thetas, results['model'], inlog=inlog)
    
    # plot the data
    axes[3].plot(mwave, mospec, color='grey', label='Obs', **kwargs)
    #plot posterior draws
    for vecs in specvecs[:-1]:
        [ax.plot(mwave, v, color=c, alpha=alpha, **kwargs) for ax,v,c in zip(axes, vecs, color)]
    #plot the minimizer result
    #[ax.plot(mwave, v, color='cyan', alpha=1.0, **kwargs) for ax,v in zip(axes, specvecs[-1])]

    [a.axhline( int(i==0), linestyle=':', color='black') for i,a in enumerate(axes[-2:])]
    if xlim is not None:
        [a.set_xlim(xlim) for a in axes]
    # add axes to the figure
    [fig.add_subplot(ax) for ax in axes]
    return fig

def zoom_spec_figure(results, zoom_region_list,
                     alpha=0.3, samples=[-1],
                     inlog=True,
                     start=0, thin=1, layout=None,
                     subplot_spec=None, xlim=None, **kwargs):

    """
    Plot zoom-ins on selected spectral regions, showing obs, model,
    and residuals.
    """
    fig = pl.figure(figsize=(10, 5))

    # posterior draws
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    thetas = [flatchain[s,:] for s in samples]
    mwave, mospec, mounc, specvecs = comp_samples(thetas, results['model'], inlog=inlog)
    
    #get the panel geometry if not given
    nobj = len(zoom_region_list)
    if layout is None:
        nx = int(np.floor(np.sqrt(nobj)))
        ny = int(np.ceil(nobj*1.0/nx))
        layout = [nx,ny]
        
    #loop over panels, one for each spectral region
    gs = gridspec.GridSpec(layout[0], layout[1])
    for i, reg in enumerate(zoom_region_list):
        x, y = i % nx, np.floor(i*1.0 / nx)
        fig, subgs = one_specregion_figure(mwave, mospec, specvecs,
                                           reg, fig, subplot_spec = gs[i], **kwargs)

    return fig

def one_specregion_figure(wave, obs, specs, reg,
                          fig, subplot_spec=None, **kwargs):

    if subplot_spec is None:
        gs = gridspec.GridSpec(2,1,height_ratios = [3,1], hspace=0)
    else:
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, hspace=0,
                                              subplot_spec=subplot_spec,
                                              height_ratios = [3,1])
    sfig = pl.Subplot(fig, gs[0,0])
    sfig.set_ylabel(r'$f_\lambda \times \, C$')
    
    pl.setp(sfig.get_xticklabels(), visible = False)

    rfig = pl.Subplot(fig, gs[1,0])
    rfig.set_ylabel(r'$\chi$')
    rfig.set_xlabel(r'$\lambda (\AA)$')
    inreg = (wave > reg[0]) & (wave < reg[1])
    #print(inreg.sum())

    for s in specs:
        sfig.plot(wave[inreg], s[3][inreg], alpha=0.3, color = 'magenta')
        rfig.plot(wave[inreg], s[5][inreg], alpha=0.3, **kwargs)
    sfig.plot(wave[inreg], obs[inreg], **kwargs)
    rfig.set_ylim(-3,3)
    fig.add_subplot(sfig)
    fig.add_subplot(rfig)
    return fig, gs


def phot_figure(results, alpha=0.3, samples = [-1],
                start=0, thin=1,
                **kwargs):
    """
    Plot the photometry for the model and data (with error bars). Then
    plot residuals
    """
    fig = pl.figure()
    gs = gridspec.GridSpec(2,1, height_ratios=[3,1])
    gs.update(hspace=0)
    phot, res = pl.Subplot(fig, gs[0]), pl.Subplot(fig, gs[1])
    res.set_ylabel( r'$\chi$')
    phot.set_ylabel('maggies')
    
    # posterior draws
    flatchain = results['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    thetas = [flatchain[s,:] for s in samples]
    mwave, mospec, mounc, specvecs = comp_samples(thetas, results['model'], photflag=1)
    #print(mwave, mospec)
    for vecs in specvecs:
        vv = vecs[0], vecs[-1]
        [ax.plot(mwave, v, color='magenta', alpha=alpha, marker='o', **kwargs)
         for ax, v in zip([phot, res], vv) ]
    
    phot.errorbar(mwave, mospec, yerr=mounc,
                  color='black')
    phot.plot(mwave, mospec, label = 'observed',
              color='black', marker='o', **kwargs)
    phot.legend(loc=0)
    res.axhline(0, linestyle=':', color='grey')
    
    fig.add_subplot(phot)
    fig.add_subplot(res)
    
    return fig
    

def calibration_figure(cal_result, nocal_result = None, samples = [-1],
                       start=0, thin=1, inlog=True,
                       alpha = 0.3, **kwargs):
    #plot the calibration vector and the reconstruction of it
    # from the f(alpha) and GP, for cal and uncal.

    fig = pl.figure()
    gs = gridspec.GridSpec(1,1)
    ax = pl.Subplot(fig, gs[0])
    ax.set_ylabel(r'$\mu/e^{f(\alpha) + \tilde{\Delta}}$')

    if nocal_result is None:
        result = cal_result
    else:
        result = nocal_result
        
    # plot posterior draws
    flatchain = result['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    thetas = [flatchain[s,:] for s in samples]
    mwave, mospec, mounc, specvecs = comp_samples(thetas, result['model'], photflag=0)
    [ax.plot(mwave, vec[0]/vec[3], color='magenta', alpha=alpha) for vec in specvecs]
    
    if nocal_result is None:
        ax.axhline(1.0, color='grey', label='Caldwell')
    else:
        mwave = cal_result['obs']['wavelength'][mask]
        mospec_cal = cal_result['obs']['spectrum'][mask]
        mospec_nocal = nocal_result['obs']['spectrum'][mask]
        if inlog:
            caldwell = np.exp(mospec_cal - mospec_nocal)
        else:
            caldwell = mospec_cal/mospec_nocal
        ax.plot(mwave, caldwell, color='grey', label='Caldwell')
    
    ax.legend(loc=0)
    fig.add_subplot(ax)
        
    return fig
    
def corner_plot(results, showpars=None, start=0, thin=1):
    #just wrap subtriangle
    """
    Make a triangle plot of the (thinned, latter) samples of the posterior
    parameter space.  Optionally make the plot only for a supplied subset
    of the parameters.
    """

    # pull out the parameter names and flatten the thinned chains
    parnames = np.array(results['model'].theta_labels())
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
        parlabels = [pardict[p] for p in parnames]
    fig = triangle.corner(flatchain, labels = parlabels,
                          quantiles=[0.16, 0.5, 0.84], verbose=False,
                          truths = truths)

    return fig

def obsdict(obs, photflag):
    """
    Return a dictionary of observational data, generated depending on
    whether you're matching photometry or spectroscopy.
    """
    if photflag == 0:
        outn = 'spectrum'
        marker = None
    elif photflag == 1:
        outn = 'sed'
        marker = 'o'
        obs = obs.copy()
        obs['wavelength'] = np.array([f.wave_effective for f in obs['filters']])
        obs['spectrum'] = 10**(0-0.4 * obs['mags'])
        obs['unc'] = obs['mags_unc'] * obs['spectrum']
        obs['mask'] = obs['mags_unc'] > 0
        
    return obs, outn, marker

def data_figure(results_list, layout = None, shaded = False, **kwargs):
    fig = pl.figure(figsize=(10, 5))

    #get the panel geometry if not given
    nobj = len(results_list)
    if layout is None:
        nx = int(np.floor(np.sqrt(nobj)))
        ny = int(np.ceil(nobj*1.0/nx))
        layout = [nx,ny]
        
    gs = gridspec.GridSpec(layout[0], layout[1])
    for i, res in enumerate(results_list):
        x, y = i % nx, np.floor(i*1.0 / nx)
        obs = res['obs']

        if shaded:
            ax = pl.Subplot(fig, gs[x,y])
            ax = one_data_figure_shaded(obs, ax, color = 'blue', facecolor = 'grey')
            fig.add_subplot(ax)
        else:
            fig, sub_gs = one_data_figure_sep(obs, fig, subplot_spec = gs[i], **kwargs)
            
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


def gp_figure(result, start=0, thin=1,
              inlog=True, samples=[-1]):
    fig = pl.figure(figsize = (10,6))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1,1,1], hspace=0, wspace=0.05)
    
    # Axis variables
    ylabels = [r'$(\sigma/d)^2$', r'$\ln d$', r'$\ln \mu$', r'$f(\alpha)$', r'$\Delta$ (residual)',
                r'$\chi$ (linear)']
    color = ['magenta',  'grey', 'blue', 'green', 'red', 'orange']
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

    # make posterior draws
    flatchain = result['chain'][:,start::thin,:]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])
    thetas = [flatchain[s,:] for s in samples]
    thetas += [result['initial_center']]

    mwave, mospec, mounc, specvecs = comp_samples(thetas, result['model'], inlog=inlog)
    for vec in specvecs:
        vv = [(mounc/mospec)**2, np.log(mospec), np.log(vec[0]), np.log(vec[1]), vec[2], vec[-1]]
        [ax.plot(mwave,v, color=c, alpha=0.3) for ax, v, c in zip(axes[2:], vv[2:], color[2:])]
    [ax.plot(mwave,v, color=c, alpha=1.0) for ax, v, c in zip(axes[:2], vv[:2], color[:2])]
    [fig.add_subplot(ax) for ax in axes]
    return fig

def comparison_figure(results_list):
    pass

if __name__ == '__main__':
    figext = '.png'
    nsample = 5
    samples = np.random.uniform(0, 1, size=nsample)
    showpars_phys = ['mass', 'tage', 'zmet', 'dust2', 'sigma_smooth']
    showpars_cal = ['dust2', 'spec_norm','poly_coeffs1', 'poly_coeffs2', 'gp_length', 'gp_amplitude', 'gp_jitter']
    zoom_regions = [[3750,4100.], [6500, 6600.], [5850, 5950], [5000, 5400]]
    
    results = []
    rdir = '/Users/bjohnson/Projects/cetus/results/'
    #res = [rdir+'b192-g242.020.cal_1405648278.sampler01',
    #       rdir+'b192-g242.020.nocal_1405677518.sampler01']
    #res = [rdir+'b192-g242.225.cal_1407376313.sampler01']
    #res = [rdir+'b192-g242.225.cal_1407608570.sampler01']
    #res = [rdir+'b192-g242.225.cal_1409443437.sampler01']
    res = [#rdir+'b192-g242.225.cal_1407608570.sampler01']
           #rdir+'b192-g242.225.cal_1409443437.sampler01']
           #rdir+'b192-g242.225.cal_1409477803.sampler01',
           #rdir+'b192-g242.225.cal_1409534549.sampler01',
           rdir+'b192-g242.225.cal_1412345250.sampler01']
    inlog = True

    name = ['B192 cal.', 'B192 no cal.']
    
    for i,r in enumerate(res):
        sf, mf = r+'_mcmc', r+'_model'
        result, pr, model = diagnostics.read_pickles(sf, model_file = mf)
        result['model'] = model
        #best = np.argmin([p.fun for p in powell_guesses])
        #result['optimizer_results'] = pr[best]
        of = result['run_params']['outfile'].split('/')[-1]
        of = of.replace('.','_')
        ns = result['chain'].shape[0] * result['chain'].shape[1]
        sample = [int(s * ns) for s in samples]

        gfig = gp_figure(result, samples=sample, inlog=inlog)
        gfig.axes[0].set_yscale('log')
        gfig.savefig('gp_figure.png')
        
        
        sfig = spec_figure(result, samples=sample,
                        linewidth = 0.5, xlim = (3650, 7300),
                        inlog=inlog)
        
        sfig.suptitle(name[i])
        sfig.savefig('sfig_'+ of + figext)
        pl.close(sfig)
        
        
        pfig = phot_figure(result, samples=sample)
        pfig.suptitle(name[i])
        pfig.savefig('pfig_' + of + figext)
        pl.close(pfig)

        tfig = corner_plot(result, showpars = showpars_phys, start=-1250)
        tfig.suptitle(name[i])
        tfig.savefig('ptri_' + of + figext)
        pl.close(tfig)

        tfig = corner_plot(result, showpars = showpars_cal, start=-1250)
        tfig.suptitle(name[i])
        tfig.savefig('ctri_' + of + figext)
        pl.close(tfig)

        zfig = zoom_spec_figure(result, zoom_regions,
                                inlog=inlog,
                                samples=sample, color='grey')
        zfig.savefig('zfig_'+ of + figext)
        pl.close(zfig)
               
        results += [result]

    results += [None]
    #sys.exit()
    #dfig = data_figure(results,
    #                   color='b', linewidth=0.5)
    #[dfig.axes[i*2].set_title(name[i]) for i in range(len(name))]
    #dfig.savefig('dfig_'+ '_'.join(of.split('_')[0:2]) + figext)

    cfig = calibration_figure(results[0], nocal_result = results[1],
                              samples=sample, inlog=inlog)
    cfig.savefig('cfig_'+ '_'.join(of.split('_')[0:2]) + figext)



    print(of)
