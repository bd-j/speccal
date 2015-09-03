#!/bin/bash

# speconly_ideal
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_ideal.py \
       --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
       --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5 \
       --nwalkers=126 --niter=1024 --do_powell=False \
       --noisefactor=10.0

# speconly_calibrated, lo S/N
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=False \
  --noisefactor=10.0

# speconly_uncalibrated, lo S/N
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_linear.py \
       --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
       --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.u0.t9.0_z0.0_a0.5 \
       --nwalkers=254 --niter=6144 --do_powell=False \
       --noisefactor=10.0

# photonly
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_photonly.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_photonly.c0.t9.0_z0.0_a0.5 \
  --nwalkers=64 --niter=2048 --do_powell=False 

# somephot linear
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/mock_somephot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_griz_linear.u0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True \
  --noisefactor=10.0

# specphot linear
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True  \
  --noisefactor=10.0

# specphot linear noisy
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u1.t12.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.u1.t12.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True \
  --noisefactor=10.0

# specphot linear vary params

#real data
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/real_specphot.py \
  --objname=NGC1851 --noisefactor=5.0 --calibrated=True \
  --outfile=$PROJECTS/speccal/code/results/ggc_ngc1851 \
  --nwalkers=126 --niter=1024 --do_powell=True \
  --noisefactor=5
