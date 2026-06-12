import pybdsim
import os
import pytest
import pandas as pd

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

from pprint import pformat

print(
    "FIELD_OFFSET_CONFIGS = "
    + pformat(
        pd.DataFrame(FIELD_OFFSET_CONFIGS).to_dict(orient="list"),
        sort_dicts=False,
        width=120
    )
)

FIELD_OFFSET_CONFIGS = {
    'name': ['x_offset', 'y_offset', 'z_offset'],
    'FIELD_TYPE': ['bmap3d', 'bmap3d', 'bmap3d'],
    'FIELD_FORMAT': ['bdsim3d', 'bdsim3d', 'bdsim3d'],
    'FIELD_FILE': ['3dexample.dat', '3dexample.dat', '3dexample.dat'],
    'FIELD_LENGTH': ['1.0', '1.0', '1.0'],
    'BEAM_ENERGY': ['10.0', '10.0', '10.0'],
    'FIELD_X_OFFSET': ['0.1', '0.0', '0.0'],
    'FIELD_Y_OFFSET': ['0.0', '0.1', '0.0'],
    'FIELD_Z_OFFSET': ['0.0', '0.0', '0.1'],
    'expected_samples': [5000, 5000, 5000],
    'expected_x_sigma': [0.013835341613313606, 0.01438563542490261, 0.013459971494015222],
    'expected_xp_sigma': [0.0072268511635527445, 0.007424089790444462, 0.007452970502992219],
    'expected_x_mean': [0.568285196685791, 0.577986688876152, 0.5479260095000267],
    'expected_y_sigma': [0.00478033296891485, 0.005241687258850288, 0.004177156229900975],
    'expected_yp_sigma': [0.0023960091258065264, 0.002609618738995901, 0.002217752745742983],
    'expected_y_mean': [0.23385638386309146, 0.25584089472591875, 0.2142822514116764]}

config_df = pd.DataFrame(FIELD_OFFSET_CONFIGS)
configs = config_df.to_dict(orient="records")
param_keys = ["FIELD_TYPE", "FIELD_FORMAT", "FIELD_FILE", "FIELD_LENGTH", "BEAM_ENERGY",
              "FIELD_X_OFFSET", "FIELD_Y_OFFSET", "FIELD_Z_OFFSET"]

simulation_results = []

@pytest.fixture(scope="module", autouse=True)
def write_results_csv():
    yield

    pd.DataFrame(simulation_results).to_csv(
        "field_offset_results.csv",
        index=False,
    )

@pytest.mark.parametrize("config", configs, ids=[c['name'] for c in configs])
def test_field_offset(config):

    os.chdir(os.path.dirname(__file__))

    base_name = f"field_{config['name']}"
    template_name = "field_map_offset.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    params = {k: config[k] for k in param_keys}

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, params)
    pybdsim.Run.Bdsim(gmad_name, base_name, ngenerate=5000, seed=1)

    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    sampler_data = data.get_sampler("d2.")

    sampler_number = len(sampler_data)
    sampler_x = sampler_data['x']
    sampler_xp = sampler_data['xp']
    sampler_y = sampler_data['y']
    sampler_yp = sampler_data['yp']

    sampler_x_sigma = sampler_x.std()
    sampler_xp_sigma = sampler_xp.std()
    sampler_x_mean = sampler_x.mean()
    sampler_y_sigma = sampler_y.std()
    sampler_yp_sigma = sampler_yp.std()
    sampler_y_mean = sampler_y.mean()

    simulation_results.append(
        {
            "name": config["name"],
            "samples": sampler_number,
            "x_sigma": sampler_x_sigma,
            "xp_sigma": sampler_xp_sigma,
            "x_mean": sampler_x_mean,
            "y_sigma": sampler_y_sigma,
            "yp_sigma": sampler_yp_sigma,
            "y_mean": sampler_y_mean,
        }
    )

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1.", "d2.", size=6, average=True)
    print(rmat)

    print(f"assert sampler_number == {sampler_number}")
    print(f"assert pytest.approx(sampler_x_sigma, rel=1e-4) == {sampler_x_sigma}")
    print(f"assert pytest.approx(sampler_xp_sigma, rel=1e-4) == {sampler_xp_sigma}")
    print(f"assert pytest.approx(sampler_x_mean, rel=1e-4) == {sampler_x_mean}")
    print(f"assert pytest.approx(sampler_y_sigma, rel=1e-4) == {sampler_y_sigma}")
    print(f"assert pytest.approx(sampler_yp_sigma, rel=1e-4) == {sampler_yp_sigma}")
    print(f"assert pytest.approx(sampler_y_mean, rel=1e-4) == {sampler_y_mean}")

    assert sampler_number == config['expected_samples']
    assert pytest.approx(sampler_x_sigma, rel=1e-4) == config['expected_x_sigma']
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == config['expected_xp_sigma']
    assert pytest.approx(sampler_x_mean, rel=1e-4) == config['expected_x_mean']
    assert pytest.approx(sampler_y_sigma, rel=1e-4) == config['expected_y_sigma']
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == config['expected_yp_sigma']
    assert pytest.approx(sampler_y_mean, rel=1e-4) == config['expected_y_mean']