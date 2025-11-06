d1: drift, l=1*m, aper1=10*cm;
AlBlock: rcol, l=2*m, material="G4_Al";

l0 : line = (d1,AlBlock);

use, period=l0;

sample, all;

option, physicsList = "g4FTFP_BERT";

option, muonSplittingFactor={{ SPLITTING_FACTOR }};

beam, particle="pi+",
      kineticEnergy=10*GeV;

sampler: samplerplacement, z=3*m, aper1=5*m, aper2=5*m, shape="rectangular", partID={13,-13};