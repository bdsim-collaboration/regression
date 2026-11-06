lb1: laserflux, wavelength=532.0*nm,
	m2=1.2,
	pulseDuration=30.0*ps,
	pulseEnergy=3*mJ,
	w0=2000*um,
	laserPolarization1=0,
	laserPolarization2=0,
	laserPolarization3=1.0;

l1: laserwire, l=2.5*m, wireLength=9*mm,laserOffsetTheta=0,
	laserOffsetPhi=0, laserOffsetX=0*m,laserOffsetY=0*m,laserOffsetZ=0*m,
	laserBeam="lb1";

dipole: rbend, l=20*m, angle=0.00226709, apertureType="elliptical", aper1=200*mm, aper2=35*mm;

sep: drift, l=80*m, apertureType="rectangular", aper1=300*mm, aper2=35*mm;

lat : line=(l1, dipole, sep);
use, lat;

sample, all;

option, physicsList="laser_cumulative_compton_scattering";

beam, particle = "e-",
      energy = 45.6*GeV,
      distrType = "reference";