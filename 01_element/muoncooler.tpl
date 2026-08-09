{%- block gmad -%}
cooldef1: coolingchannel, 
    nCoils=62,
    coilInnerRadius = {250*mm},
    coilRadialThickness = {0.001*mm},
    coilLengthZ = {0.177486988*m},
    coilCurrent = {6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01, 6072100.01, -6072100.01},
    coilOffsetX = {0.0},
    coilOffsetY = {0.0},
    coilOffsetZ = {-12.267157421, -11.732842579, -11.467157421000001, -10.932842579, -10.667157421, -10.132842579, -9.867157421000002, -9.332842579000001, -9.067157421000001, -8.532842579, -8.267157421, -7.732842579, -7.4671574210000005, -6.932842579, -6.667157421000001, -6.132842579, -5.867157421000001, -5.332842579, -5.067157421000001, -4.532842579, -4.267157421, -3.732842579, -3.467157421, -2.9328425790000003, -2.667157421, -2.1328425790000005, -1.867157421, -1.3328425790000002, -1.067157421, -0.532842579, -0.26715742100000006, 0.26715742100000006, 0.532842579, 1.067157421, 1.3328425790000002, 1.867157421, 2.1328425790000005, 2.667157421, 2.9328425790000003, 3.467157421, 3.732842579, 4.267157421, 4.532842579, 5.067157421000001, 5.332842579, 5.867157421000001, 6.132842579, 6.667157421000001, 6.932842579, 7.4671574210000005, 7.732842579, 8.267157421, 8.532842579, 9.067157421000001, 9.332842579000001, 9.867157421000002, 10.132842579, 10.667157421, 10.932842579, 11.467157421000001, 11.732842579, 12.267157421},
    coilTiltX = {0.0},
    coilTiltY = {0.0},
    coilTiltZ = {0.0},
    coilMaterial = {"G4_Cu"},
    onAxisTolerance = 1e-2,
    
    nDipoles= 62,
    dipoleAperture= {0.2}, 
    dipoleLengthZ = {0.1},
    dipoleOffsetZ = {-12.24, -11.76, -11.440000000000001, -10.96, -10.64, -10.16, -9.840000000000002, -9.360000000000001, -9.040000000000001, -8.56, -8.24, -7.76, -7.44, -6.96, -6.640000000000001, -6.16, -5.840000000000001, -5.36, -5.040000000000001, -4.5600000000000005, -4.24, -3.7600000000000002, -3.44, -2.9600000000000004, -2.64, -2.1600000000000006, -1.84, -1.36, -1.0400000000000003, -0.56, -0.24000000000000002, 0.24000000000000002, 0.56, 1.0400000000000003, 1.36, 1.84, 2.1600000000000006, 2.64, 2.9600000000000004, 3.44, 3.7600000000000002, 4.24, 4.5600000000000005, 5.040000000000001, 5.36, 5.840000000000001, 6.16, 6.640000000000001, 6.96, 7.44, 7.76, 8.24, 8.56, 9.040000000000001, 9.360000000000001, 9.840000000000002, 10.16, 10.64, 10.96, 11.440000000000001, 11.76, 12.24},
    dipoleFieldStrength = {0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2, -0.2, 0.2, 0.2, -0.2},
    dipoleEngeCoefficient = {5.5},
    
    nAbsorbers = 0,
    absorberType = {"cylinder"},
    absorberMaterial = {"G4_LITHIUM_HYDRIDE"},
    absorberOffsetZ = {0.0},
    absorberCylinderLength = {0.0},
    absorberCylinderRadius = {0.0},
    absorberWedgeOpeningAngle = {0.0},
    absorberWedgeHeight = {0.0},
    absorberWedgeRotationAngle = {0.0},
    absorberWedgeOffsetX = {0},
    absorberWedgeOffsetY = {0},
    absorberWedgeApexToBase = {0.0},
    
    nRFCavities = 0,
    rfOffsetZ = {0.0},
    rfTimeOffset = {0.0},
    rfLength={0.0},
    rfVoltage={0.0},
    rfPhase={0.0},
    rfFrequency={0.0},
    rfWindowThickness = {0*um},
    rfWindowMaterial = {"G4_Be"},
    rfWindowRadius = {0.0},
    rfCavityMaterial = {"G4_Cu"},
    rfCavityVacuumMaterial = {"vacuum"},
    rfCavityRadius = {0.0},
    rfCavityThickness = {0.0},
    
    magneticFieldModel="solenoidsheet",
    dipoleFieldModel="dipole",
    electricFieldModel="rfpillbox";

mc1: muoncooler, l=28*m, horizontalWidth=0.8*m, coolingDefinition="cooldef1";

s1: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-3200.0*mm, shape="rectangular", aper1=5*m, aper2=5*m;
s2: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-1600.0*mm, shape="rectangular", aper1=5*m, aper2=5*m;
s3: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=0.0*mm, shape="rectangular", aper1=5*m, aper2=5*m;
s4: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=1600.0*mm, shape="rectangular", aper1=5*m, aper2=5*m;
s5: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=3200.0*mm, shape="rectangular", aper1=5*m, aper2=5*m;

d1: drift, l=1*mm;
c1: rcol, l=1*mm, xsize=0*mm, ysize=0*mm, material="graphite", outerDiameter=200*cm;
lat: line=(mc1,d1,c1);

use, period=lat;

option, checkOverlaps=1,
        stopSecondaries=1,
        storeTrajectories=1,
        collimatorsAreInfiniteAbsorbers=1,
        maximumStepLength=20*mm,
        integratorSet="geant4",
        physicsList="qgsp_bic em muon";

beam, particle="mu+",
      momentum=200*MeV,
      distrType="userfile",
      distrFile="muoncoolerBeam.dat",
      distrFileFormat = "x[m]:xp[rad]:y[m]:yp[rad]:z[m]:P[GeV]",
      nlinesIgnore=1;

{%- endblock -%}

{%- block beam -%}
"x [m]" "px [rad]" "y [m]" "py [rad]" "z [m]" "P [GeV]"
0.010894217528402805 0.0 0.008229246363043785 0.0 8.8 0.210 
0.006983276456594467 0.0 0.006129929795861244 0.0 8.8 0.190
0.00862 0.0 0.00702 0.0 8.8 0.200

{%- endblock -%}
