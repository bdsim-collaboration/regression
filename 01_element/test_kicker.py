import pytest
import pybdsim
import numpy as np
import os


def test():
    os.chdir(os.path.dirname(__file__))

    base_name     = "kicker"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    hkick = 0.001   # horizontal kick [rad]
    vkick = 0.0005  # vertical kick [rad]
    data = {
        'HKICK'       : hkick,
        'VKICK'       : vkick,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("k1.")

    mean_xp = np.mean(s['xp'])
    mean_yp = np.mean(s['yp'])
    print("Expected xp:", hkick, "Measured:", mean_xp)
    print("Expected yp:", vkick, "Measured:", mean_yp)

    assert mean_xp == pytest.approx(hkick, rel=0.05)
    assert mean_yp == pytest.approx(vkick, rel=0.05)
