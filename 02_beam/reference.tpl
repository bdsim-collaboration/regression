d1: drift, l=1*m;
t1: drift, l={{ LENGTH }}*m;
d2: drift, l=1*m;

l0 : line = (d1,t1,d2);

use, period=l0;

beam,  particle="e-",
       energy={{ ENERGY }} * GeV,
       distrType="reference";


sample, all;