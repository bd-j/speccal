import os
import numpy as np

sdir = os.path.join(os.environ.get('PROJECTS'), 'speccal')
datadir= "$PROJECTS/speccal/data/ggclib/mocks/miles/"
jobscriptdir = os.path.join(sdir, 'code','jobscripts')
accounts = {'PHAT': 'TG-AST130057'}

def make_stampede_mock_job(ncpu=16, niter=1024, do_powell=True,
                           paramfile='mock_specphot_linear',
                           params='u0.t12.z0.0.a0.5',
                           partition='normal',
                           walltime=5.0,
                           jobscriptdir=jobscriptdir,
                           account='TG-AST130057',
                           noisefactor=10.0):

    nwalkers = (ncpu - 1) * 2
    wt = '{0:02.0f}:{1:02.0f}:00'.format(walltime, np.mod(walltime, 1) * 60)
    jobname = '{0}.{1}'.format(paramfile, params)
    scriptname = os.path.join(jobscriptdir, '{0}.sh'.format(jobname))
    out = open(scriptname, 'w')
    out.write("#!/bin/bash\n\n")
    out.write("###partition\n"
              "#SBATCH -p {}\n\n".format(partition))
    out.write("### Requested number of nodes\n"
              "#SBATCH -n {0}\n\n".format(ncpu))
    out.write("### Requested computing time\n"
              "#SBATCH -t {0}\n\n".format(wt))
    out.write("### Account\n"
              "#{0}\n\n".format(account))
    out.write("### Job name\n"
              "#SBATCH -J '{0}'\n\n".format(jobname))
    out.write("### output and error logs\n"
              "#SBATCH -o {0}_%j.out\n"
              "#SBATCH -e {0}_%j.err\n\n".format(jobname))
    out.write("\n ibrun python-mpi $PROJECTS/speccal/code/prospectr.py \\\n")
    out.write(" --param_file=$PROJECTS/speccal/code/paramfiles/{0}.py \\\n".format(paramfile))
    out.write(" --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.{0}.pkl \\\n".format(params))
    out.write(" --outfile=$PROJECTS/speccal/code/results/{0}_$SLURM_JOB_ID \\\n".format(jobname))
    out.write(" --nwalkers={0}"
              " --niter={1}"
              " --do_powell={2}"
              " --noisefactor={3}".format(nwalkers, niter, str(do_powell), noisefactor))
    out.close()
    return jobname, scriptname
