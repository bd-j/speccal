A Library of Integrated Spectra of Galactic Globular Clusters

Ricardo P. Schiavon, James A. Rose, Stephane Courteau & Lauren MacArthur

Abstract: We present a new library of integrated spectra of 40 Galactic
globular clusters, obtained with the Blanco 4-m telescope and the R-C
spectrograph at the Cerro Tololo Interamerican Observatory.  The spectra
cover the range ~ 3350 - 6430 Angstrom with ~ 3.1 Angstrom (FWHM)
resolution. The spectroscopic observations and data reduction were
designed to integrate the full projected area within the cluster core
radii in order to properly sample the light from stars in all relevant
evolutionary stages. The S/N values of the flux-calibrated spectra range
from 50 to 240/Angstrom at 4000 Angstrom and from 125 to 500/Angstrom
at 5000 Anstrom. The selected targets span a wide range of cluster
parameters, including metallicity, horizontal-branch morphology, Galactic
coordinates, Galactocentric distance, and concentration. The total sample
is thus fairly representative of the entire Galactic globular cluster
population and should be valuable for comparison with similar integrated
spectra of unresolved stellar populations in remote systems. For most of
the library clusters, our spectra can be coupled with deep color-magnitude
diagrams and reliable metal abundances from the literature to enable
the calibration of stellar population synthesis models.


This README file contains basic information for users of the spectral
library. For further details on the observations and data reduction,
please refer to Schiavon et al. 2005, ApJS, XXX, XXX.


1) For each cluster the final flux-calibrated spectrum is presented in
a single FITS file. These files are named NGC****_x_n.fits, where:
     x = a,b,c refers to different exposures of the same cluster.
     n = 1,2 refers to different aperture extractions
Please, refer to the paper for further details.

2) For each of the above files, there is a corresponding auxiliary
FITS file in multi-spectrum format, where the spectra are NOT
flux-calibrated. The latter files are named NGC****_x_n.aux.fits. The
content of these multi-spectrum FITS files is as follows:

  -- The first band contains the variance-weighted, CR-cleaned,
sky-subtracted, spectrum.

  -- The second band contains the unweighted, uncleaned, 
sky-subtracted, spectrum.

  -- The third band contains the sky spectrum obtained from the separate
sky frame, except for few the cases where sky subtraction was performed
using the ends of the slit. This is the case of NGC 6352, 6553, and 7089.
For those clusters, the sky spectrum recorded come from the ends of the
slit.

  -- The fourth band contains the S/N ratio per spectral pixel, obtained
by dividing the sky-subtracted spectrum by the quadrature sum of the
variance spectrum in the cluster and in the sky.


