import pandas as pd
import pybdsim
import os
import pytest

FIELD_TRACK_CONFIGS = {
    'name': ['1d', '2d', '3d', '4d'],
    'FIELD_TYPE': ['bmap1d', 'bmap2d', 'bmap3d', 'bmap4d'],
    'FIELD_FORMAT': ['bdsim1d', 'bdsim2d', 'bdsim3d', 'bdsim4d'],
    'FIELD_FILE': ['1dexample.dat', '2dexample.dat', '3dexample.dat', '4dexample.dat'],
    'FIELD_LENGTH': ['1.0', '1.0', '1.0', '1.0'],
    'BEAM_ENERGY': ['1.0', '1.0', '10.0', '1.0'],
    'expected_samples': [5000, 5000, 5000, 5000],
    'expected_x_sigma': [0.005376337095341065, 0.018419807674184555, 0.014271013981376675, 0.0017761945113131978],
    'expected_xp_sigma': [0.0020095770501316365, 0.012275748730702854, 0.007416345178360901, 0.00035136689967708905],
    'expected_x_mean': [0.009675747349709855, 0.00010576244183662311, 0.5799057118058205, 1.3215959309269466e-05],
    'expected_y_sigma': [0.013080372518385144, 0.0028903558879167703, 0.004409168245179408, 0.0017857934386989697],
    'expected_yp_sigma': [0.007200602146071492, 0.0032919311481443892, 0.002207146467824476, 0.0003543900315124393],
    'expected_y_mean': [-0.4976539144992828, 1.0439962574434957e-05, 0.226743698656559, -2.1371774069270798e-05],
    'expected_energy_sigma': [0.019783911947098665, 0.019783911947098665, 0.19783912688382108, 0.019783911947098665],
    'expected_energy_mean': [0.9995454897522926, 0.9995454897522926, 9.995454892158508, 0.9995454897522926],
}

config_df = pd.DataFrame(FIELD_TRACK_CONFIGS)
configs = config_df.to_dict(orient="records")
param_keys = ["FIELD_TYPE", "FIELD_FORMAT", "FIELD_FILE", "FIELD_LENGTH", "BEAM_ENERGY"]

simulation_results = []

@pytest.fixture(scope="module", autouse=True)
def write_results_csv():
    yield

    pd.DataFrame(simulation_results).to_csv(
        "field_tracking_results.csv",
        index=False,
    )

@pytest.mark.parametrize("config", configs, ids=[c['name'] for c in configs])
def test_field_tracking(config):

    os.chdir(os.path.dirname(__file__))

    base_name = f"field_{config['name']}"
    template_name = "field_map_track.tpl"
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