import pybdsim
import os
import pytest

def test(test_length, testlength_primaries, testdata_store):

    os.chdir(os.path.dirname(__file__))

    base_name     = "field_2d_nearest"
    template_name = "field_map_interpolator.tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    params = {
        'FIELD_TYPE': 'bmap2d',
        'FIELD_FORMAT': 'bdsim2d',
        'FIELD_FILE': '2dexample.dat',
        'FIELD_LENGTH': '1.0',
        'BEAM_ENERGY': '1.0',
        'FIELD_INTERPOLATOR': 'nearest',
    }

    nprimary = testlength_primaries.get_nprimary(__file__, test_length)

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, params)
    pybdsim.Run.Bdsim(gmad_name, base_name, nprimary, 1)

    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    sampler_data = data.get_sampler("d2.")

    sampler_Sigma_x  = sampler_data['x'].std()
    sampler_Sigma_xp = sampler_data['xp'].std()
    sampler_mean_x   = sampler_data['x'].mean()
    sampler_Sigma_y  = sampler_data['y'].std()
    sampler_Sigma_yp = sampler_data['yp'].std()
    sampler_mean_y   = sampler_data['y'].mean()

    assert pytest.approx(sampler_Sigma_x, rel=1e-3)  == 0.003369653184855525
    assert pytest.approx(sampler_Sigma_xp, rel=1e-3) == 0.0017251188722198771
    assert pytest.approx(sampler_mean_x, rel=1e-3)   == 3.561127408404445e-05
    assert pytest.approx(sampler_Sigma_y, rel=1e-3)  == 0.021222188911178278
    assert pytest.approx(sampler_Sigma_yp, rel=1e-3) == 0.021194220357673686
    assert pytest.approx(sampler_mean_y, rel=1e-3)   == -0.00018837782195287218

    testdata_store.add_test_output(__file__, root_name, "root", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_x, "Sigma_x", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_xp, "Sigma_xp", nprimary)
    testdata_store.add_test_object(__file__, sampler_mean_x, "mean_x", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_y, "Sigma_y", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_yp, "Sigma_yp", nprimary)
    testdata_store.add_test_object(__file__, sampler_mean_y, "mean_y", nprimary)
