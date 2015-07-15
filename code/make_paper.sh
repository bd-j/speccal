#!bin/sh

# define a version
version=""

# First we make the mocks with fiducial parameters
fiducial_params = {'tage':12.0, 'zmet':0.0, 'dust2':0.5}
fall = "t{tage:3.1f}_z{zmet:3.1f}_a{dust2:3.1f}".format(**fiducial_params)

ft, fa, fz = make_ggc_mocks.main(fiducial_params=fiducial+params,
				 do_noise=True, do_vary=True,
				 version=version)

### Then we build jobsripts
import make_jobscripts
make_job = make_jobscripts.make_stampede_mock_job

# speconly_ideal
jobn, ideal_script = make_job(ncpu=128, niter=1024,
			      walltime=3.0,
  			      paramfile='mock_speconly_ideal',
			      params='c0.{0}'.format(fall),
			      do_powell=False, noisefactor=10)

# speconly_calibrated, lo S/N
jobn, sonlycal_script = make_job(ncpu=128, niter=6144,
				 walltime=5.0,
  			         paramfile='mock_speconly_linear',
			         params='c0.{0}'.format(fall),
			         do_powell=False, noisefactor=10)

# speconly_uncalibrated, lo S/N
jobn, sonlycal_script = make_job(ncpu=128, niter=6144,
  			         paramfile='mock_speconly_linear',
			         params='u0.{0}'.format(fall),
			         do_powell=False, noisefactor=10)

# photonly
jobn, photonly_script = make_job(ncpu=64, niter=1024,
				 walltime=2,
  			         paramfile='mock_photonly',
			         params='c0.{0}'.format(fall),
			         do_powell=False, noisefactor=10)

# specphot uncalibrated
jobn, specphot_script = make_job(ncpu=128, niter=2048,
  			         param_file='mock_specphot_linear',
			         params='c0.{0}'.format(fall),
			         do_powell=False, noisefactor=10)
# noise_runs

# physical parameter variation runs

# nphot runs

# And run them
#./runjobs.sh

# and assign output names to variables
specphot_ideal="results/ggc_mock_ideal.c0.t9.0_z0.0_a0.5_1430261146_mcmc"
speconly_cal="results/ggc_mock_speconly.c0.t9.0_z0.0_a0.5_1430808300_mcmc"
speconly_uncal="results/ggc_mock_speconly.u0.t9.0_z0.0_a0.5_1431313829_mcmc"
specphot_uncal="results/ggc_mock_specphot_linear.u0.t9.0_z0.0_a0.5_5280432_1431898211_mcmc"
photonly="results/ggc_mock_photonly.c0.t9.0_z0.0_a0.5_1430274922_mcmc"

noise_runs="ls results/ggc_mock_specphot_linear.u[1-9]*.${fall}_*_mcmc"
age_runs="ls results/ggc_mock_specphot_linear.u0.t*_z${fz}_a${fa}_*_mcmc"
z_runs="ls results/ggc_mock_specphot_linear.u0.t${ft}_z*_a${fa}_*_mcmc"
dust_runs="ls results/ggc_mock_specphot_linear.u0.t${ft}_z${fz}_a*_*_mcmc"

real_cal="results/ggc_ngc1851_1432787257_mcmc"
real_uncal="results/ggc_ngc1851_uncal_tightprior_1433448554_mcmc"


### Then we make figures for mocks ###
#fig 1
python plot_data.py
#fig 2
python plot_ideal_dashboard.py $specphot_ideal green "Ideal, Perfectly Known calibration"
cp "${specphot_ideal%_mcmc}.dashboard.pdf" ../tex/figures/ideal.pdf
#fig 3
python plot_dashboard.py $speconly_cal blue "Spectroscopy Only, Calibrated data"
cp "${speconly_cal%_mcmc}.dashboard.pdf" ../tex/figures/speconly_calibrated.pdf
#fig 4
python plot_dashboard.py $speconly_uncal blue "Spectroscopy Only, Unalibrated data"
cp "${speconly_uncal%_mcmc}.dashboard.pdf" ../tex/figures/speconly_uncalibrated.pdf
#fig 5
python plot_photonly_dashboard.py $photonly red "Photometry Only"
cp "${photonly%_mcmc}.dashboard.pdf" ../tex/figures/photonly.pdf
#fig 6
python plot_dashboard.py $specphot_uncal maroon "Spectroscopy and Photometry, Uncalibrated data"
cp "${specphot_uncal%_mcmc}_dashboard.pdf" ../tex/figures/specphot_uncalibrated.pdf
# Fig components
python plot_components.py $specphot_uncal maroon
cp "${specphot_uncal%_mcmc}.components.pdf" ../tex/figures/components.pdf
#fig 7
python plot_joint_posterior.py $speconly_uncal $photonly $specphot_uncal
#fig 8
python plot_noise_realizations.py $noise_runs
#fig 9
python plot_phys_summary $age_runs $z_runs $dust_runs


#fig 10
python plot_real_dashboard.py $real_cal
#fig 11
python plot_real_dashboard.py $real_uncal

#Then we compile the paper
cd ../tex
make -B
