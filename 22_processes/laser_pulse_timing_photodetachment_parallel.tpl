! configuration file based on documentation from S. Alden in "BDSIM Laser Documentation" (2024)
! defining the laserflux 
lb1: laserflux, wavelength=1064.0*nm,
m2=1.8, 
pulseDuration=100.0*ps,
pulseEnergy={{LASER_ENERGY}}*mJ,
w0=1.0*mm;

! defining the parallel laserwire
laser1: laserwire, l=10.0*m, wireLength=10.0*m, laserOffsetTheta={{LASER_OFFSET_THETA}}, laserOffsetPhi=0, laserOffsetY=0.0*m, laserBeam="lb1";

! defining a drift (not sure it's needed)
drift1: drift, l=3.280275*m;

! defining the lattice
lattice: line=(laser1);
use, lattice;

! defining the beam
beam, particle="ion 1 1 -1",
energy={{BEAM_ENERGY}}*GeV,
distrType="gauss",
sigmaT=2.0e-10*s;

! defining options
option, physicsList="laser_photo_detachment",
storeSamplerCharge=1,
storeTrajectory=1,
trajConnect=1, 
scaleFactorLaser=100, 
storeTrajectoryDepth=1, 
storeTrajectoryParticleID="11",
storeTrajectoryLocal=1,
storeTrajectorySecondaryParticles=1,
storeTrajectoryTime=1,
storeTrajectoryLinks=1,
storeTrajectoryAllVariables=1;

! defining samplers
sample, range=laser1;