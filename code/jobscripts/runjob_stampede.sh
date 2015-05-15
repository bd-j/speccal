#!/bin/bash

### Name of the job

### Requested number of nodes
#SBATCH -n 256
### Requested computing time in minutes
#SBATCH -t 4:00:00
###partition
#SBATCH -p normal
### memory per cpu, in MB

### Account
### PHAT
#SBATCH -A TG-AST140054

### Job name
#SBATCH -J 'speccal_mock'
### output and error logs
#SBATCH -o speconly_uncal_mock_%j.out
#SBATCH -e speconly_uncal_mock_%j.err


#python make_ggc_mock --tage=1.0 --zmet=0.0 --dust2=0.5 --add_noise=False --apply_cal=True


#specphot
ibrun python-mpi $PROJECTS/speccal/code/prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_$SLURM_JOB_ID \
  --nwalkers=510 --niter=2048 --do_powell=True 

