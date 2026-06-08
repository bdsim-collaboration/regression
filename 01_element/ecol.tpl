d1: drift, l=1*m;
ecol1: ecol, l={{ LENGTH }}*m, material="graphite", xsize={{ XSIZE }}*m, ysize={{ YSIZE }}*m;
d2: drift, l=1*m;

l0 : line = (d1,ecol1,d2);

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
      betx={{ BETX }}*m,
      bety={{ BETY }}*m,
      dispx=0.0*m,
      dispxp=0.0,
      dispy=0.0*m,
      dispyp=0.0,
      distrType="gausstwiss",
      emitx={{ EMITX }}*m,
      emity={{ EMITY }}*m,
      sigmaE=0.0,
      sigmaT=1e-11;
