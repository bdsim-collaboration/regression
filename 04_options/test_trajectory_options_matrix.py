import os

import pybdsim
import pytest


TEST_CASES = [
    ("storeTrajectory", "storeTrajectory=1"),
    ("storeTrajectories", "storeTrajectories=1"),
    ("storeTrajectoryDepth", "storeTrajectory=1, storeTrajectoryDepth=1"),
    ("storeTrajectoryELossSRange", 'storeTrajectory=1, storeTrajectoryELossSRange="0.0:10.0"'),
    ("storeTrajectoryEnergyThreshold", "storeTrajectory=1, storeTrajectoryEnergyThreshold=0.01*GeV"),
    ("storeTrajectoryParticle", 'storeTrajectory=1, storeTrajectoryParticle="e-"'),
    ("storeTrajectoryParticleID", 'storeTrajectory=1, storeTrajectoryParticleID="11"'),
    ("storeTrajectorySecondaryParticles", "storeTrajectory=1, storeTrajectorySecondaryParticles=1"),
    ("storeTrajectorySamplerID", 'storeTrajectory=1, storeTrajectorySamplerID="t1"'),
    ("storeTrajectoryTransportationSteps", "storeTrajectory=1, storeTrajectoryTransportationSteps=0"),
    ("trajNoTransportation", "storeTrajectory=1, trajNoTransportation=1"),
    ("trajectoryConnect", "storeTrajectory=1, trajectoryConnect=1"),
    (
        "trajectoryFilterLogicAND",
        'storeTrajectory=1, storeTrajectoryParticleID="11", storeTrajectorySamplerID="t1", trajectoryFilterLogicAND=1',
    ),
    ("trajCutGTZ", "storeTrajectory=1, trajCutGTZ=0.0*m"),
    ("trajCutLTR", "storeTrajectory=1, trajCutLTR=10.0*m"),
    ("storeTrajectoryAllVariables", "storeTrajectory=1, storeTrajectoryAllVariables=1"),
    ("storeTrajectoryIon", "storeTrajectory=1, storeTrajectoryIon=1"),
    ("storeTrajectoryKineticEnergy", "storeTrajectory=1, storeTrajectoryKineticEnergy=0"),
    ("storeTrajectoryLocal", "storeTrajectory=1, storeTrajectoryLocal=1"),
    ("storeTrajectoryLinks", "storeTrajectory=1, storeTrajectoryLinks=1"),
    ("storeTrajectoryMaterial", "storeTrajectory=1, storeTrajectoryMaterial=1"),
    ("storeTrajectoryMomentumVector", "storeTrajectory=1, storeTrajectoryMomentumVector=1"),
    ("storeTrajectoryProcesses", "storeTrajectory=1, storeTrajectoryProcesses=1"),
    ("storeTrajectoryStepPoints", "storeTrajectory=1, storeTrajectoryStepPoints=2"),
    (
        "storeTrajectoryStepPointLast",
        "storeTrajectory=1, storeTrajectoryStepPoints=2, storeTrajectoryStepPointLast=1",
    ),
    ("storeTrajectoryTime", "storeTrajectory=1, storeTrajectoryTime=1"),
]


@pytest.mark.parametrize("label,option_lines", TEST_CASES)
def test_trajectory_option(label, option_lines):
    os.chdir(os.path.dirname(__file__))

    base_name = f"trajectory_options_matrix_{label}"
    template_name = "trajectory_options_matrix.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"
    optics_name = base_name + "_optics.root"

    data = {
        "LENGTH": "1.0",
        "BEAM_ENERGY": "1",
        "OPTION_LINES": option_lines,
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 100, 1)

    assert os.path.exists(root_name)

    pybdsim.Run.RebdsimOptics(root_name, optics_name)

    assert os.path.exists(optics_name)
