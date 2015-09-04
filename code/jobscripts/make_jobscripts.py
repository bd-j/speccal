import os
import numpy as np

sdir = os.path.join(os.environ.get('PROJECTS'), 'speccal')
datadir= "$PROJECTS/speccal/data/ggclib/mocks/miles/"
jobscriptdir = os.path.join(sdir, 'code','jobscripts')
accounts = {'PHAT': 'TG-AST130057'}

def make_mock_job(ncpu=8, niter=1024, nwalkers=32, do_powell=True,
                  paramfile='mock_specphot_linear',
                  params='u0.t12.z0.0.a0.5', noisefactor=10.0,
                  machine='stampede', account='TG-AST130057',
                  partition='normal', walltime=5.0,
                  jobscriptdir=jobscriptdir):

    jobname = '{0}.{1}'.format(paramfile, params)
    scriptname = os.path.join(jobscriptdir, '{0}.sh'.format(jobname))
    out = open(scriptname, 'w')

    if machine =='stampede':
        hdr, outfile, nwalkers = stampede_header(ncpu=ncpu, account=account,
                                                 partition=partition,
                                                 walltime=walltime,
                                                 jobname=jobname)
    elif machine == 'mac':
        assert ncpu <= 8
        nwalkers -= np.mod(nwalkers, (ncpu-1) * 2)
        hdr = ("#!/bin/bash\n\n"
               "mpirun -np {0} python $PROJECTS/speccal/code/prospectr.py \\\n".format(ncpu))
        outfile = "$PROJECTS/speccal/code/results/{}".format(jobname)

    
    out.write(hdr)
    out.write(" --param_file=$PROJECTS/speccal/code/paramfiles/{0}.py \\\n".format(paramfile))
    out.write(" --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.{0}.pkl \\\n".format(params))
    out.write(" --outfile={0} \\\n".format(outfile))
    out.write(" --nwalkers={0}"
              " --niter={1}"
              " --do_powell={2}"
              " --noisefactor={3}".format(nwalkers, niter, str(do_powell), noisefactor))
    out.close()
    return jobname, scriptname


def make_real_job(ncpu=8, niter=1024, nwalkers=32, do_powell=False,
                  paramfile='real_specphot_otherjitter',
                  objname='NGC1851', calibrated=True,
                  noisefactor=5.0,
                  machine='stampede', account='TG-AST130057',
                  partition='normal', walltime=5.0,
                  jobscriptdir=jobscriptdir):
    """Make a jobscript for real spectra.
    """
    
    jobname = '{0}.{1}.cal{2}'.format(paramfile, objname, calibrated)
    scriptname = os.path.join(jobscriptdir, '{0}.sh'.format(jobname))
    out = open(scriptname, 'w')
    if machine =='stampede':
        hdr, outfile, nwalkers = stampede_header(ncpu=ncpu, account=account,
                                                 partition=partition,
                                                 walltime=walltime,
                                                 jobname=jobname)
    elif machine == 'mac':
        assert ncpu <= 8
        nwalkers -= np.mod(nwalkers, (ncpu-1) * 2)
        hdr = ("#!/bin/bash\n\n"
               "mpirun -np {0} python $PROJECTS/speccal/code/prospectr.py \\\n".format(ncpu))
        outfile = "$PROJECTS/speccal/code/results/{}".format(jobname)

    out.write(hdr)
    out.write(" --param_file=$PROJECTS/speccal/code/paramfiles/{0}.py \\\n".format(paramfile))
    out.write(" --objname={0} --calibrated={1} \\\n".format(objname, str(calibrated)))
    out.write(" --outfile={0} \\\n".format(outfile))
    out.write(" --nwalkers={0}"
              " --niter={1}"
              " --do_powell={2}"
              " --noisefactor={3}".format(nwalkers, niter, str(do_powell), noisefactor))
    out.close()
    return jobname, scriptname


def stampede_header(ncpu=16, account='TG-AST130057',
                    partition='normal', walltime=5.0,
                    jobname=None):

    nwalkers = (ncpu - 1) * 2
    wt = '{0:02.0f}:{1:02.0f}:00'.format(walltime, np.mod(walltime, 1) * 60)
    outfile = "$PROJECTS/speccal/code/results/{0}_$SLURM_JOB_ID".format(jobname)

    hdr = ''
    hdr += "#!/bin/bash\n\n"
    hdr += ("###partition\n"
            "#SBATCH -p {}\n\n").format(partition)
    hdr += ("### Requested number of nodes\n"
            "#SBATCH -n {0}\n\n").format(ncpu)
    hdr += ("### Requested computing time\n"
            "#SBATCH -t {0}\n\n").format(wt)
    hdr += ("### Account\n"
            "#SBATCH -A {0}\n\n").format(account)
    hdr += ("### Job name\n"
            "#SBATCH -J '{0}'\n\n").format(jobname)
    hdr += ("### output and error logs\n"
            "#SBATCH -o {0}_%j.out\n"
            "#SBATCH -e {0}_%j.err\n\n").format(jobname)
    hdr += "\n ibrun python-mpi $PROJECTS/speccal/code/prospectr.py \\\n"

    return hdr, outfile, nwalkers
