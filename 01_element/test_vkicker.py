import pytest
import pybdsim
import numpy as np
import os


def test():
    os.chdir(os.path.dirname(__file__))

    base_name     = "vkicker"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    vkick = 0.001   # vertical kick [rad], positive = +y direction
    data = {
        'VKICK'       : vkick,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("v1.")

    mean_yp = np.mean(s['yp'])
    print("Expected yp deflection:", vkick)
    print("Measured mean yp:", mean_yp)

    assert mean_yp == pytest.approx(vkick, rel=0.05)
