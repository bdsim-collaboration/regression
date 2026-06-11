import pybdsim
import os
import pytest

def test():

    os.chdir(os.path.dirname(__file__))

    # === Setup filenames ===
    base_name = "field_1d"
    template_name = "field_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    # === Test parameters ===
    params = {
        'FIELD_TYPE': 'bmap1d',
        'FIELD_FORMAT': 'bdsim1d',
        'FIELD_FILE': '1dexample.dat',
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

    print("assert sampler_number == ", sampler_number)

    print("assert pytest.approx(sampler_x_sigma, rel=1e-4) == ", sampler_x_sigma)
    print("assert pytest.approx(sampler_xp_sigma, rel=1e-4) == ", sampler_xp_sigma)
    print("assert pytest.approx(sampler_x_mean, rel=1e-4) == ", sampler_x_mean)

    print("assert pytest.approx(sampler_y_sigma, rel=1e-4) == ", sampler_y_sigma)
    print("assert pytest.approx(sampler_yp_sigma, rel=1e-4) == ", sampler_yp_sigma)
    print("assert pytest.approx(sampler_y_mean, rel=1e-4) == ", sampler_y_mean)

    print("assert pytest.approx(sampler_energy_sigma, rel=1e-4) == ", sampler_energy_sigma)
    print("assert pytest.approx(sampler_energy_mean, rel=1e-4) == ", sampler_energy_mean)

    assert sampler_number == 5000
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 0.005376337095341065
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.0020095770501316365
    assert pytest.approx(sampler_x_mean, rel=1e-4) == 0.009675747349709855
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.013080372518385144
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.007200602146071492
    assert pytest.approx(sampler_y_mean, rel=1e-4) == -0.4976539144992828
    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.019783911947098665
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.9995454897522926