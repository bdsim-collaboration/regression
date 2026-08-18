
d1: drift, l={{ LENGTH }}*m;


l0 : line = (d1);

use, period=l0;

sample, all;

beam, particle="{{P_TYPE}}",
      kineticEnergy={{ BEAM_ENERGY }}*GeV,
      Xp0=0,
      Yp0=0,
      distrType="reference"; 
      
option,  uprootCompatible=1, storeSamplerPolarCoords=1, storeSamplerAll=1; 
