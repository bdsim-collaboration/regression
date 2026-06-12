import pybdsim
import os
import pytest

FIELD_OFFSET_CONFIGS = [
    {
        'name': 'x_offset',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_x_offset': '0.1',
        'field_y_offset': '0.0',
        'field_z_offset': '0.0',
        'expected_samples': 5000,
        'expected_x_sigma': 0.013835341613313606,
        'expected_xp_sigma': 0.0072268511635527445,
        'expected_x_mean': 0.568285196685791,
        'expected_y_sigma': 0.00478033296891485,
        'expected_yp_sigma': 0.0023960091258065264,
        'expected_y_mean': 0.23385638386309146,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': 'y_offset',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_x_offset': '0.0',
        'field_y_offset': '0.1',
        'field_z_offset': '0.0',
        'expected_samples': 5000,
        'expected_x_sigma': 0.01438563542490261,
        'expected_xp_sigma': 0.007424089790444462,
        'expected_x_mean': 0.577986688876152,
        'expected_y_sigma': 0.005241687258850288,
        'expected_yp_sigma': 0.002609618738995901,
        'expected_y_mean': 0.25584089472591875,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': 'z_offset',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_x_offset': '0.0',
        'field_y_offset': '0.0',
        'field_z_offset': '0.1',
        'expected_samples': 5000,
        'expected_x_sigma': 0.013459971494015222,
        'expected_xp_sigma': 0.007452970502992219,
        'expected_x_mean': 0.5479260095000267,
        'expected_y_sigma': 0.004177156229900975,
        'expected_yp_sigma': 0.002217752745742983,
        'expected_y_mean': 0.2142822514116764,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
]

@pytest.mark.parametrize("config", FIELD_OFFSET_CONFIGS, ids=[c['name'] for c in FIELD_OFFSET_CONFIGS])
def test_field_offset(config):

    os.chdir(os.path.dirname(__file__))

    base_name = f"field_{config['name']}"
    template_name = "field_map_offset.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    params = {
        'FIELD_TYPE': config['field_type'],
        'FIELD_FORMAT': config['field_format'],
        'FIELD_FILE': config['field_file'],
        'FIELD_LENGTH': config['field_length'],
        'BEAM_ENERGY': config['beam_energy'],
        'FIELD_X_OFFSET': config['field_x_offset'],
        'FIELD_Y_OFFSET': config['field_y_offset'],
        'FIELD_Z_OFFSET': config['field_z_offset'],
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

    pass