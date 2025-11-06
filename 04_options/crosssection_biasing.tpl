d1: drift, l=1*m, aper1=10*cm;
AlBlock: rcol, l=2*m, material="G4_Al";

l0 : line = (d1,AlBlock);

use, period=l0;

sample, all;

option, physicsList = "g4FTFP_BERT";

piPlusBias:   xsecBias, particle="pi+", proc="Decay", xsecfact={{ BIAS_FACTOR }}, flag=1;
piMinusBias:  xsecBias, particle="pi-", proc="Decay", xsecfact={{ BIAS_FACTOR }}, flag=1;

option, biasForWorldVolume="piPlusBias piMinusBias";
option, biasForWorldContents="piPlusBias piMinusBias";
option, biasForWorldVacuum = "piPlusBias piMinusBias";

beam, particle="pi+",
      kineticEnergy=10*GeV;

sampler: samplerplacement, z=3*m, aper1=5*m, aper2=5*m, shape="rectangular", partID={13,-13};