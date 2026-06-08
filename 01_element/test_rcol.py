import pytest
import pybdsim
import numpy as np
import os


def test_absorption():
    """Verify that a narrow rcol absorbs most of a wide beam."""
    os.chdir(os.path.dirname(__file__))

    base_name     = "rcol"
    template_name = base_name + ".tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    n_generate = 3000
    # Beam with sigma_x ~ 7 mm (betx=100m, emitx=5e-7m -> sigma=sqrt(5e-5)~7mm)
    # Collimator xsize = 1 mm << sigma_x => significant absorption
    data = {
        'LENGTH'      : 0.5,
        'XSIZE'       : 0.001,   # 1 mm half-aperture
        'YSIZE'       : 0.05,    # 50 mm half-aperture (wide, passes all y)
        'BETX'        : 100,
        'BETY'        : 4,
        'EMITX'       : 5e-7,
        'EMITY'       : 5e-7,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, n_generate, 1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("rcol1.")
    n_survived = len(s)

    print("Generated:", n_generate)
    print("Survived:", n_survived)
    print("Transmission fraction:", n_survived / n_generate)

    # With sigma_x ~ 7 mm and xsize = 1 mm, transmission << 50 %
    assert n_survived < 0.5 * n_generate


def test_transmission():
    """Verify that a wide rcol passes a narrow beam with near-100% efficiency."""
    os.chdir(os.path.dirname(__file__))

    base_name     = "rcol_wide"
    template_name = "rcol.tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    n_generate = 3000
    # Beam with sigma_x = sigma_y ~ 1.4 mm (betx=bety=4m, emitx=emity=5e-7m)
    # Collimator xsize = ysize = 5 cm >> sigma => nearly all particles pass
    data = {
        'LENGTH'      : 0.5,
        'XSIZE'       : 0.05,    # 50 mm half-aperture >> sigma_x ~ 1.4 mm
        'YSIZE'       : 0.05,    # 50 mm half-aperture >> sigma_y ~ 1.4 mm
        'BETX'        : 4,
        'BETY'        : 4,
        'EMITX'       : 5e-7,
        'EMITY'       : 5e-7,
        'BEAM_ENERGY' : '1'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, n_generate, 2)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("rcol1.")
    n_survived = len(s)

    print("Generated:", n_generate)
    print("Survived:", n_survived)
    print("Transmission fraction:", n_survived / n_generate)

    # With sigma << xsize, expect near-100% transmission
    assert n_survived > 0.90 * n_generate
