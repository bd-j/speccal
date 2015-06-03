#!/bin/bash


#python make_ggc_mock --tage=1.0 --zmet=0.0 --dust2=0.5 --add_noise=False --apply_cal=True


#speconly_ideal
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_ideal.py \
       --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
       --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5 \
       --nwalkers=126 --niter=1024 --do_powell=False 

# speconly_calibrated, lo S/N
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_losn.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=False


#speconly_uncalibrated, lo S/N
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_speconly_losn.py \
       --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
       --outfile=$PROJECTS/speccal/code/results/ggc_mock_speconly.u0.t9.0_z0.0_a0.5 \
       --nwalkers=254 --niter=6144 --do_powell=False

# photonly
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_photonly.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.c0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_photonly.c0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=2048 --do_powell=False 

#specphot
mpirun -np 5 python prospectr.py --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot.u0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True --gptype=george


#specphot linear
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u0.t9.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True

#specphot linear noisy
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/mock_specphot_linear.py \
  --filename=$PROJECTS/speccal/data/ggclib/mocks/miles/ggc_mock.u1.t12.0_z0.0_a0.5.pkl \
  --outfile=$PROJECTS/speccal/code/results/ggc_mock_specphot_linear.u1.t12.0_z0.0_a0.5 \
  --nwalkers=126 --niter=1024 --do_powell=True


#real data
mpirun -np 5 python prospectr.py \
  --param_file=$PROJECTS/speccal/code/paramfiles/real_specphot.py \
  --objname=NGC1851 --noisefactor=5.0 --calibrated=True \
  --outfile=$PROJECTS/speccal/code/results/ggc_ngc1851 \
  --nwalkers=126 --niter=1024 --do_powell=True
