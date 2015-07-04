#!bin/sh

# First we make the mocks with fiducial parameters
ft="9.0"
fz="0.0"
fa="0.5"
fall="t${ft}_z${fz}_a${fa}"
#python make_ggc_mocks #$ft $fz $fa

### Then we fit the mocks and assign output names to variables
#./runjobs.sh
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
