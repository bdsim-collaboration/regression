import pybdsim
import os
import pytest

def test(test_length, testlength_primaries, testdata_store):

    os.chdir(os.path.dirname(__file__))

    base_name     = "field_3d_linearmag"
    template_name = "field_map_interpolator.tpl"
    gmad_name     = base_name + ".gmad"
    root_name     = base_name + ".root"

    params = {
        'FIELD_TYPE': 'bmap3d',
        'FIELD_FORMAT': 'bdsim3d',
        'FIELD_FILE': '3dexample.dat',
        'FIELD_LENGTH': '1.0',
        'BEAM_ENERGY': '10.0',
        'FIELD_INTERPOLATOR': 'linearmag',
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

    assert pytest.approx(sampler_Sigma_x, rel=1e-3)  == 0.014327077885869297
    assert pytest.approx(sampler_Sigma_xp, rel=1e-3) == 0.007516079445364859
    assert pytest.approx(sampler_mean_x, rel=1e-3)   == 0.5661105772137642
    assert pytest.approx(sampler_Sigma_y, rel=1e-3)  == 0.004486432469495695
    assert pytest.approx(sampler_Sigma_yp, rel=1e-3) == 0.002270183744965747
    assert pytest.approx(sampler_mean_y, rel=1e-3)   == 0.22307491642832755

    testdata_store.add_test_output(__file__, root_name, "root", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_x, "Sigma_x", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_xp, "Sigma_xp", nprimary)
    testdata_store.add_test_object(__file__, sampler_mean_x, "mean_x", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_y, "Sigma_y", nprimary)
    testdata_store.add_test_object(__file__, sampler_Sigma_yp, "Sigma_yp", nprimary)
    testdata_store.add_test_object(__file__, sampler_mean_y, "mean_y", nprimary)
