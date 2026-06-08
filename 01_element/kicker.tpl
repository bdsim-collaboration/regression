d1: drift, l=1*m;
k1: kicker, hkick={{ HKICK }}, vkick={{ VKICK }};
d2: drift, l=1*m;

l0 : line = (d1,k1,d2);

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
      emitx=1e-9*m,
      emity=1e-9*m,
      sigmaE=0.0,
      sigmaT=1e-11;
