import pybdsim
import os
import pytest

FIELD_CONFIGS = [
    {
        'name': '1d',
        'field_type': 'bmap1d',
        'field_format': 'bdsim1d',
        'field_file': '1dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'expected_samples': 5000,
        'expected_x_sigma': 0.005376337095341065,
        'expected_xp_sigma': 0.0020095770501316365,
        'expected_x_mean': 0.009675747349709855,
        'expected_y_sigma': 0.013080372518385144,
        'expected_yp_sigma': 0.007200602146071492,
        'expected_y_mean': -0.4976539144992828,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '2d',
        'field_type': 'bmap2d',
        'field_format': 'bdsim2d',
        'field_file': '2dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'expected_samples': 5000,
        'expected_x_sigma': 0.018419807674184555,
        'expected_xp_sigma': 0.012275748730702854,
        'expected_x_mean': 0.00010576244183662311,
        'expected_y_sigma': 0.0028903558879167703,
        'expected_yp_sigma': 0.0032919311481443892,
        'expected_y_mean': 1.0439962574434957e-05,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '3d',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',  # Higher energy needed for strong 3D fields
        'expected_samples': 5000,
        'expected_x_sigma': 0.014271013981376675,
        'expected_xp_sigma': 0.007416345178360901,
        'expected_x_mean': 0.5799057118058205,
        'expected_y_sigma': 0.004409168245179408,
        'expected_yp_sigma': 0.002207146467824476,
        'expected_y_mean': 0.226743698656559,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': '4d',
        'field_type': 'bmap4d',
        'field_format': 'bdsim4d',
        'field_file': '4dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'expected_samples': 5000,
        'expected_x_sigma': 0.0017761945113131978,
        'expected_xp_sigma': 0.00035136689967708905,
        'expected_x_mean': 1.3215959309269466e-05,
        'expected_y_sigma': 0.0017857934386989697,
        'expected_yp_sigma': 0.0003543900315124393,
        'expected_y_mean': -2.1371774069270798e-05,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
]


@pytest.mark.parametrize("config", FIELD_CONFIGS, ids=[c['name'] for c in FIELD_CONFIGS])
def test_field_tracking(config):

    os.chdir(os.path.dirname(__file__))

    base_name = f"field_{config['name']}"
    template_name = "field_map_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    params = {
        'FIELD_TYPE': config['field_type'],
        'FIELD_FORMAT': config['field_format'],
        'FIELD_FILE': config['field_file'],
        'FIELD_LENGTH': config['field_length'],
        'BEAM_ENERGY': config['beam_energy']
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

    assert sampler_number == config['expected_samples']
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == config['expected_x_sigma']
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == config['expected_xp_sigma']
    assert pytest.approx(sampler_x_mean, rel=1e-4) == config['expected_x_mean']
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == config['expected_y_sigma']
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == config['expected_yp_sigma']
    assert pytest.approx(sampler_y_mean, rel=1e-4) == config['expected_y_mean']
    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == config['expected_energy_sigma']
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == config['expected_energy_mean']