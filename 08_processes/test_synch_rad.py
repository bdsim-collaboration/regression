import pytest
import pybdsim
import os

import synch_rad_analysis

def test(geant4_version, bdsim_version,
         test_length, testlength_primaries, testdata_store) :

    os.chdir(os.path.dirname(__file__))

    base_name     = "synch_rad"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"

    nprimary = testlength_primaries.get_nprimary(__file__,test_length)
    pybdsim.Run.Bdsim(gmad_name, base_name, ngenerate=nprimary, seed=1)

    ac = synch_rad_analysis.analysis(root_name)
    ac.write_persistent_data("synch_rad.json")
