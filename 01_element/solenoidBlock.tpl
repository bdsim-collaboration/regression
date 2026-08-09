{%- block gmad -%}
cooldef1: coolingchannel, 
    nCoils=1,
    coilInnerRadius = {0.25*m},
    coilRadialThickness = {0.1693*m},
    coilLengthZ = {0.14*m},
    coilCurrent = {500.0},
    coilOffsetX = {0*m},
    coilOffsetY = {0.0*cm},
    coilOffsetZ = {0.0*m},
    coilTiltX ={0.0*degrees},
    coilTiltY = {0.0*degrees},
    coilTiltZ = {0.0*degrees},
    coilMaterial = {"G4_Cu"},
    nSheets = 20, 
    onAxisTolerance = 1e-5,
    nDipoles = 0,
    dipoleAperture= {0.0},
    dipoleLengthZ = {0.0},
    dipoleFieldStrength = {0.0},
    dipoleEngeCoefficient = {0},
    dipoleOffsetZ = {0.0},
    
    nAbsorbers = 0,
    absorberType = {"cylinder"},
    absorberMaterial = {"G4_LITHIUM_HYDRIDE"},
    absorberOffsetZ = {0.0},
    absorberCylinderLength = {6.537*cm},
    absorberCylinderRadius = {15*cm},
    absorberWedgeOpeningAngle = {0},
    absorberWedgeHeight = {0},
    absorberWedgeRotationAngle = {0},
    absorberWedgeOffsetX = {0},
    absorberWedgeOffsetY = {0},
    absorberWedgeApexToBase = {0},
    
    nRFCavities = 0,
    rfOffsetZ = {0},
    rfTimeOffset={0},
    rfLength={0},
    rfVoltage={0},
    rfPhase={0.0},
    rfFrequency={0.0},
    rfWindowThickness = {0*um},
    rfWindowMaterial = {"G4_Be"},
    rfWindowRadius = {15*cm},
    rfCavityMaterial = {"G4_Cu"},
    rfCavityVacuumMaterial = {"vacuum"},
    rfCavityRadius = {163*mm},
    rfCavityThickness = {3*mm},
    
    magneticFieldModel="solenoidblock",
    electricFieldModel="rfpillbox",
    dipoleFieldModel="dipoleenge";

mc1: muoncooler, l=1.0*m, horizontalWidth=2.*m, coolingDefinition="cooldef1";

c1: rcol, l=1*mm, xsize=0*mm, ysize=0*mm, material="graphite", outerDiameter=200*cm;

d1: drift, l=1*mm;

lat: line=(mc1,d1,c1);

use, period=lat;

sample, all;

option, checkOverlaps=1,
        stopSecondaries=1,
        collimatorsAreInfiniteAbsorbers=1,
        maximumStepLength=20*mm,
        integratorSet="geant4",
        physicsList="qgsp_bic em muon";

beam, particle="mu+",
      momentum=200*MeV,
      distrType="userfile",
      distrFile="solenoidBlockBeam.dat",
      distrFileFormat = "x[mm]:xp[rad]:y[mm]:yp[rad]:z[mm]:Ek[MeV]",
      nlinesIgnore=1;
{%- endblock -%}

{%- block sampler -%}
X Y
0.0 0.0
0.010566326789557934 0.0034909385722130537
0.021383726969361305 0.006933522876352072
0.032705288380384445 0.010242943651974201
0.04478023946285248 0.0132580092176795
0.05782969295978546 0.015696242451667786
0.07198864966630936 0.017108118161559105
0.08719359338283539 0.016847513616085052
0.10299701988697052 0.01410119142383337
0.11831708252429962 0.008057610131800175
0.13121087849140167 -0.0016753420932218432
0.13890157639980316 -0.01433532778173685
0.1384236216545105 -0.027060799300670624
0.12826643884181976 -0.03410504758358002
0.11147365719079971 -0.026813294738531113
0.10156432539224625 0.0038850316777825356
0.13186359405517578 0.05205060914158821
0.21759480237960815 0.0464925542473793
0.2112365961074829 -0.03917374461889267
0.13047851622104645 0.03701428696513176
0.18474081158638 -0.04443580284714699
{%- endblock -%}

{%- block beam -%}
"x [mm]" "px [rad]" "y [mm]" "py [rad]" "z [mm]" "Ek [MeV]"
0.0 0.0 0.0 0.0 0.0 120.53555547653541 
10.0 0.0 0.0 0.0 0.0 120.53555547653541 
20.0 0.0 0.0 0.0 0.0 120.53555547653541
30.0 0.0 0.0 0.0 0.0 120.53555547653541
40.0 0.0 0.0 0.0 0.0 120.53555547653541
50.0 0.0 0.0 0.0 0.0 120.53555547653541
60.0 0.0 0.0 0.0 0.0 120.53555547653541
70.0 0.0 0.0 0.0 0.0 120.53555547653541
80.0 0.0 0.0 0.0 0.0 120.53555547653541
90.0 0.0 0.0 0.0 0.0 120.53555547653541
100.0 0.0 0.0 0.0 0.0 120.53555547653541
110.0 0.0 0.0 0.0 0.0 120.53555547653541
120.0 0.0 0.0 0.0 0.0 120.53555547653541
130.0 0.0 0.0 0.0 0.0 120.53555547653541
140.0 0.0 0.0 0.0 0.0 120.53555547653541
150.0 0.0 0.0 0.0 0.0 120.53555547653541
160.0 0.0 0.0 0.0 0.0 120.53555547653541
170.0 0.0 0.0 0.0 0.0 120.53555547653541
180.0 0.0 0.0 0.0 0.0 120.53555547653541
190.0 0.0 0.0 0.0 0.0 120.53555547653541
200.0 0.0 0.0 0.0 0.0 120.53555547653541
{%- endblock -%}
