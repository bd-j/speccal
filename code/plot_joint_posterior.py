import sys, pickle
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec

from bsfh import read_results as bread
from bsfh import sps_basis
from bsfh.gp import GaussianProcess



if __name__ == "__main__":
    photonly = ''
    speconly = ''
    specphot = ''

    results = [bread.read_pickles(rfile, model_file=rfile.replace('mcmc','model'))[0]
               for rfile in resfiles]

    showpars = np.array(['mass', 'tage', 'zmet', 'dust2'])

    gs2 = gridspec.GridSpec(4, 2)
    for i, p1 in enumerate(showpars):
        for j, p2 in enumerate(showpars[i:])
            k = j+i
            ax = axes[i, k]
            pdf = joint_pdf(res, p1, p2)
            
    
