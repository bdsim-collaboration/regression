import pytest
import pybdsim
import numpy as np
import os


def test():
    """Verify that an rmatrix element applies the specified 4x4 transfer
    matrix at its midpoint, producing the expected compound transfer map
    (drift L/2) * R_custom * (drift L/2)."""
    os.chdir(os.path.dirname(__file__))

    base_name     = "rmatrix"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"
    optics_name   = base_name + "_optics.root"

    l      = 2.0
    rmat21 = -0.5   # thin horizontal defocusing: x' += rmat21 * x
    rmat43 =  0.5   # thin vertical defocusing:   y' += rmat43 * y
    data = {
        'LENGTH'      : l,
        'RMAT21'      : rmat21,
        'RMAT43'      : rmat43,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)
    pybdsim.Run.RebdsimOptics(root_name, optics_name)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1.", "rm1.", size=6, average=True)

    # The rmatrix element applies R_custom at its midpoint flanked by two drift(l/2) sections.
    # Compound transfer: T = drift(l/2) * R_custom * drift(l/2)
    half_l = l / 2.0
    D = np.array([[1, half_l], [0, 1]])          # 2x2 drift
    Rx = np.array([[1,      0], [rmat21, 1]])    # thin H kick
    Ry = np.array([[1,      0], [rmat43, 1]])    # thin V kick

    Tx = D @ Rx @ D   # horizontal block
    Ty = D @ Ry @ D   # vertical block

    ref_rmat = [[Tx[0, 0], Tx[0, 1], 0,        0,        0, 0],
                [Tx[1, 0], Tx[1, 1], 0,        0,        0, 0],
                [0,        0,        Ty[0, 0], Ty[0, 1], 0, 0],
                [0,        0,        Ty[1, 0], Ty[1, 1], 0, 0],
                [0,        0,        0,        0,        1, 0],
                [0,        0,        0,        0,        0, 0]]

    print(pybdsim.Testing.round_matrix(rmat, 3))
    print(pybdsim.Testing.round_matrix(ref_rmat, 3))
    print(pybdsim.Testing.max_matrix_diff(rmat, ref_rmat))

    assert pybdsim.Testing.compare_matrix(rmat, ref_rmat, 5e-2)
