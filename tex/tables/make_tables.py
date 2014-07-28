import diagnostics

pardict = {}
pardict[r'$m_*$']='mass'
pardict[r'$t$']='tage'
pardict[r'$Z$']='zmet'
pardict[r'$A_V$']='dust2'
pardict[r'$\sigma$']='sigma_smooth'
pardict[r'$z$']='zred'
pardict[r'$s$']='gp_jitter'
pardict[r'$a$']='gp_amplitude'
pardict[r'$l$']='gp_length'
#pardict[r'$c_2$']='poly_coeffs'
#pardict[r'$c_1$']='poly_coeffs'
             
def write_prior_table(plist, pardict, dir = './'):
    
    pnames = [p['name'] for p in plist]
    
    freelines =[]
    fixedlines =[]
    fmt_free = "{:10s} & {:10s} & {:7.3f} & {:7.3f} \\\ \n"
    fmt_fixed = "{:10s} & {10s} & {:7.3f} \\\ \n"
    
    for k,v in pardict.iteritems():
        ind = pnames.index(v)
        unit = plist[ind].get('units','\\nodata')
        if plist[ind]['isfree']:
            tmp = plist[ind]['prior_args']
            mini, maxi = tmp['mini'], tmp['maxi']
            freelines += [fmt_free.format(k, unit, mini, maxi)]
        else:
            val = plist[ind]['value']
            fixedlines += [fmt_fixed.format(k, unit, val)]
            

    prior_table = open(dir+'priors.tbl','w')
    prior_table.write("""
\\begin{deluxetable}{lccc}
\\tablecaption{Prior Parameter Ranges \\label{tbl:priors}}
\\tablehead{
\\colhead{Parameter} & \\colhead{Units} & \\colhead{Min.} & \\colhead{Max.}}
\\startdata
""")
    for l in freelines:
        prior_table.write(l)
    prior_table.write("""
\\enddata
\\end{deluxetable}
"""
    )
    prior_table.close()
    return freelines, fixedlines
    
if __name__ == '__main__':

    results = []
    rdir = '/Users/bjohnson/Projects/cetus/results/'
    res = [rdir+'b192-g242.020.cal_1405648278.sampler01',
           rdir+'b192-g242.020.nocal_1405677518.sampler01']
    name = ['B192 cal.', 'B192 no cal.']
    
    sf, mf = res[0]+'_mcmc', res[0]+'_model'
    result, pr, model = diagnostics.read_pickles(sf, model_file = mf)
    plist = result['plist']
    free, fixed = write_prior_table(plist, pardict)
