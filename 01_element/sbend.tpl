d1: drift, l=1*m;
t1: sbend, l={{ LENGTH }}*m, angle={{ ANGLE }}, horizontalWidth={{ HORIZONTAL_WIDTH }}*m , tilt={{ TILT1 }}, e1={{ E1 }}, e2={{ E2 }} , h1={{ H1 }}, h2= {{ H2 }}, fint= {{ FINT }}, fintx= {{ FINTX}}, fintK2= {{ FINTK2 }}, fintxK2={{ FINTXK2}}, hgap= {{ HGAP }};
d2: drift, l=1*m;

l0 : line = (d1, t1, d2);

use, period=l0;

sample, all;

option, sampleElementsWithPoleface=1;

beam, particle="{{ PARTICLE_TYPE }}",
      kineticEnergy={{ BEAM_ENERGY }}*MeV,
      X0=0.0*m,
      Xp0=0.0,
      Y0=0.0*m,
      Yp0=0.0, 
      alfx=0,
      alfy=0,
      betx=4*m,
      bety=4*m,
      dispx=0.0*m,
      dispxp=0.0,
      dispy=0.0*m,
      dispyp=0.0,
      distrType="gausstwiss",
      emitx=5e-7*m,
      emity=5e-7*m,
      sigmaE=0.0,
      sigmaT=1e-11;

option, integratorSet= "{{ INTEGRATOR }}"; !"geant4"; !"bdsimmatrix",
