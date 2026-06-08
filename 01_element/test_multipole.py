import pytest
import pybdsim
import numpy as np
from numpy import sin, cos, sinh, cosh, sqrt
import os


def test():
    """Verify that a multipole with a pure quadrupole component matches
    the analytical thick-lens quadrupole transfer matrix."""
    os.chdir(os.path.dirname(__file__))

    base_name     = "multipole"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"
    optics_name   = base_name + "_optics.root"

    l   = 1.0
    k1  = -1.0   # k1 < 0 => horizontally focusing in BDSIM convention
    # knl[0] is the integrated quadrupole strength (k1 * L)
    k1L = k1 * l
    data = {
        'LENGTH'      : l,
        'KNL'         : '{' + str(k1L) + '}',
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1)
    pybdsim.Run.RebdsimOptics(root_name, optics_name)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1.", "m1.", size=6, average=True)

    # Analytical thick-lens quadrupole matrix (k1 < 0 = horizontal focusing)
    a = sqrt(abs(k1)) * l
    ref_rmat = [[cos(a),              1/sqrt(abs(k1))*sin(a),  0,                       0,                      0, 0],
                [-sqrt(abs(k1))*sin(a), cos(a),                0,                       0,                      0, 0],
                [0,                    0,                       cosh(a),                 1/sqrt(abs(k1))*sinh(a), 0, 0],
                [0,                    0,                       sqrt(abs(k1))*sinh(a),   cosh(a),                 0, 0],
                [0,                    0,                       0,                       0,                       1, 0],
                [0,                    0,                       0,                       0,                       0, 0]]

    print(pybdsim.Testing.round_matrix(rmat, 3))
    print(pybdsim.Testing.round_matrix(ref_rmat, 3))
    print(pybdsim.Testing.max_matrix_diff(rmat, ref_rmat))

    assert pybdsim.Testing.compare_matrix(rmat, ref_rmat, 5e-2)
