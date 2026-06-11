import pybdsim
import os
import pytest

def test():

    os.chdir(os.path.dirname(__file__))

    # === Setup filenames ===
    base_name = "field_3d"
    template_name = "field_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    # === Test parameters ===
    params = {
        'FIELD_TYPE': 'bmap3d',
        'FIELD_FORMAT': 'bdsim3d',
        'FIELD_FILE': '3dexample.dat',
        'FIELD_LENGTH': '1.0',
        'BEAM_ENERGY': '10.0'
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
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 0.014271013981376675
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.007416345178360901
    assert pytest.approx(sampler_x_mean, rel=1e-4) == 0.5799057118058205
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.004409168245179408
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.002207146467824476
    assert pytest.approx(sampler_y_mean, rel=1e-4) == 0.226743698656559
    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.19783912688382108
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 9.995454892158508



