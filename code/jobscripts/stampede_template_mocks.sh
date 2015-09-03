#!/bin/bash

### Name of the job

### Requested number of nodes
#SBATCH -n ${ncpu}
### Requested computing time in minutes
#SBATCH -t 4:00:00
###partition
#SBATCH -p normal
### memory per cpu, in MB

### Account
### startup
###SBATCH -A TG-AST140054
### PHAT
###SBATCH -A TG-AST130057
### advance
#SBATCH -A TG-AST150015

### Job name
#SBATCH -J 'speccal_mock'
### output and error logs
#SBATCH -o specphot_mock_lin_%j.out
#SBATCH -e specphot_mock_lin_%j.err

###python make_ggc_mock --tage=1.0 --zmet=0.0 --dust2=0.5 --add_noise=False --apply_cal=True

### params = u0.t9.0_z0.0_a0.5
###specphot
ibrun python-mpi $PROJECTS/speccal/code/prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.${params}.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.${params}_$SLURM_JOB_ID \
  --nwalkers=${nwalker} --niter=${niter} --do_powell=True 

