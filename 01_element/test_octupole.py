import pytest
import pybdsim
import numpy as np
import os


def test():
    os.chdir(os.path.dirname(__file__))

    base_name     = "octupole"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"
    optics_name   = base_name + "_optics.root"

    l  = 0.5
    k3 = 10.0   # weak octupole: nonlinear kicks are negligible at test emittances
    data = {
        'LENGTH'      : l,
        'K3'          : k3,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)
    pybdsim.Run.RebdsimOptics(root_name, optics_name)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1.", "oct1.", size=6, average=True)

    # An octupole has no linear kicks, so to first order it acts as a drift of length l
    ref_rmat = [[1, l, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, l, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0]]

    print(pybdsim.Testing.round_matrix(rmat, 3))
    print(pybdsim.Testing.round_matrix(ref_rmat, 3))
    print(pybdsim.Testing.max_matrix_diff(rmat, ref_rmat))

    assert pybdsim.Testing.compare_matrix(rmat, ref_rmat)
