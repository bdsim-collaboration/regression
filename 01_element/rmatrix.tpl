d1: drift, l=1*m;
rm1: rmatrix, l={{ LENGTH }}*m,
     rmat11=1, rmat12=0, rmat13=0, rmat14=0,
     rmat21={{ RMAT21 }}, rmat22=1, rmat23=0, rmat24=0,
     rmat31=0, rmat32=0, rmat33=1, rmat34=0,
     rmat41=0, rmat42=0, rmat43={{ RMAT43 }}, rmat44=1;
d2: drift, l=1*m;

l0 : line = (d1,rm1,d2);

use, period=l0;

sample, all;

beam, particle="e-",
      energy={{ BEAM_ENERGY }}*GeV,
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
