import pytest
import pybdsim
import numpy as np
import os


def test():
    os.chdir(os.path.dirname(__file__))

    base_name     = "hkicker"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    hkick = 0.001   # horizontal kick [rad], positive = +x direction
    data = {
        'HKICK'       : hkick,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("h1.")

    mean_xp = np.mean(s['xp'])
    print("Expected xp deflection:", hkick)
    print("Measured mean xp:", mean_xp)

    assert mean_xp == pytest.approx(hkick, rel=0.05)
