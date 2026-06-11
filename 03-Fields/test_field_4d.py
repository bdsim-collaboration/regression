import pybdsim
import os
import pytest

def test():

    os.chdir(os.path.dirname(__file__))

    # === Setup filenames ===
    base_name = "field_4d"
    template_name = "field_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    # === Test parameters ===
    params = {
        'FIELD_TYPE': 'bmap4d',
        'FIELD_FORMAT': 'bdsim4d',
        'FIELD_FILE': '4dexample.dat',
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
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 0.0017761945113131978
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.00035136689967708905
    assert pytest.approx(sampler_x_mean, rel=1e-4) == 1.3215959309269466e-05
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.0017857934386989697
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.0003543900315124393
    assert pytest.approx(sampler_y_mean, rel=1e-4) == -2.1371774069270798e-05
    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.019783911947098665
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.9995454897522926


