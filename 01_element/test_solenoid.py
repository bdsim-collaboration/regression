import pytest
import pybdsim
import numpy as np
from numpy import cos, sin
import os


def test():
    os.chdir(os.path.dirname(__file__))

    base_name     = "solenoid"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"
    optics_name   = base_name + "_optics.root"

    l  = 2.0
    ks = 0.5   # solenoid strength B0/Brho [1/m]
    data = {
        'LENGTH'      : l,
        'KS'          : ks,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, 5000, 1)
    pybdsim.Run.RebdsimOptics(root_name, optics_name)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1.", "sol1.", size=6, average=True)

    # Analytical solenoid transfer matrix
    # K = ks/2, phi = K*L
    # C = cos(phi), S = sin(phi)
    K   = ks / 2.0
    phi = K * l
    C   = cos(phi)
    S   = sin(phi)

    ref_rmat = [[C**2,        S*C/K,    S*C,      S**2/K, 0, 0],
                [-K*S*C,      C**2,    -K*S**2,   S*C,    0, 0],
                [-S*C,       -S**2/K,  C**2,      S*C/K,  0, 0],
                [K*S**2,     -S*C,    -K*S*C,     C**2,   0, 0],
                [0,           0,       0,          0,      1, 0],
                [0,           0,       0,          0,      0, 0]]

    print(pybdsim.Testing.round_matrix(rmat, 3))
    print(pybdsim.Testing.round_matrix(ref_rmat, 3))
    print(pybdsim.Testing.max_matrix_diff(rmat, ref_rmat))

    assert pybdsim.Testing.compare_matrix(rmat, ref_rmat, 5e-2)
