import pybdsim
import os
import pytest

def test():

    os.chdir(os.path.dirname(__file__))

    # === Setup filenames ===
    base_name = "field_2d"
    template_name = "field_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    # === Test parameters ===
    params = {
        'FIELD_TYPE': 'bmap2d',
        'FIELD_FORMAT': 'bdsim2d',
        'FIELD_FILE': '2dexample.dat',
        'FIELD_LENGTH': '1.0',
        'BEAM_ENERGY': '1.0'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, params)
    pybdsim.Run.Bdsim(gmad_name, base_name, ngenerate=5000, seed=1)

    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    sampler_data = data.get_sampler("d2.")

    sampler_number = len(sampler_data)
    sampler_x = sampler_data['x']
    sampler_xp = sampler_data['xp']
    sampler_y = sampler_data['y']
    sampler_yp = sampler_data['yp']
    sampler_energy = sampler_data['energy']

    sampler_x_sigma = sampler_x.std()
    sampler_xp_sigma = sampler_xp.std()
    sampler_x_mean = sampler_x.mean()
    sampler_y_sigma = sampler_y.std()
    sampler_yp_sigma = sampler_yp.std()
    sampler_y_mean = sampler_y.mean()
    sampler_energy_sigma = sampler_energy.std()
    sampler_energy_mean = sampler_energy.mean()

    assert sampler_number == 5000
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 0.018419807674184555
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.012275748730702854
    assert pytest.approx(sampler_x_mean, rel=1e-4) == 0.00010576244183662311
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.0028903558879167703
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.0032919311481443892
    assert pytest.approx(sampler_y_mean, rel=1e-4) == 1.0439962574434957e-05
    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.019783911947098665
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.9995454897522926


