import pybdsim
import os
import pytest
import pandas as pd

FIELD_INTERPOLATION_CONFIGS = [
    {
        'name': '1d_linear',
        'field_type': 'bmap1d',
        'field_format': 'bdsim1d',
        'field_file': '1dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linear',
        'expected_samples': 5000,
        'expected_x_sigma': 0.005354202580646126,
        'expected_xp_sigma': 0.001999076865966775,
        'expected_x_mean': 0.009670783941715035,
        'expected_y_sigma': 0.013084708676590241,
        'expected_yp_sigma': 0.007202701938688855,
        'expected_y_mean': -0.49687208319306375,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '1d_linearmag',
        'field_type': 'bmap1d',
        'field_format': 'bdsim1d',
        'field_file': '1dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linearmag',
        'expected_samples': 5000,
        'expected_x_sigma': 0.005360657975868763,
        'expected_xp_sigma': 0.0020020719479500963,
        'expected_x_mean': 0.009683479112656596,
        'expected_y_sigma': 0.01307943167235869,
        'expected_yp_sigma': 0.0072005010788836685,
        'expected_y_mean': -0.49720062289237976,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '1d_nearest',
        'field_type': 'bmap1d',
        'field_format': 'bdsim1d',
        'field_file': '1dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'nearest',
        'expected_samples': 5000,
        'expected_x_sigma': 0.0018571098686965446,
        'expected_xp_sigma': 0.0004469812469452971,
        'expected_x_mean': 0.009439205919555388,
        'expected_y_sigma': 0.01307691510397638,
        'expected_yp_sigma': 0.007198931645154899,
        'expected_y_mean': -0.49778312349915504,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '2d_linear',
        'field_type': 'bmap2d',
        'field_format': 'bdsim2d',
        'field_file': '2dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linear',
        'expected_samples': 5000,
        'expected_x_sigma': 0.015804539022049573,
        'expected_xp_sigma': 0.010316564235008652,
        'expected_x_mean': 8.499772373153859e-05,
        'expected_y_sigma': 0.0021579719729259575,
        'expected_yp_sigma': 0.0023403134007770198,
        'expected_y_mean': 1.8904866732009396e-05,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '2d_linearmag',
        'field_type': 'bmap2d',
        'field_format': 'bdsim2d',
        'field_file': '2dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linearmag',
        'expected_samples': 5000,
        'expected_x_sigma': 0.01969824900808272,
        'expected_xp_sigma': 0.013199531233735034,
        'expected_x_mean': 0.00010222111806351676,
        'expected_y_sigma': 0.0038697662117538945,
        'expected_yp_sigma': 0.004129952429779972,
        'expected_y_mean': 4.843972614912673e-06,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '2d_nearest',
        'field_type': 'bmap2d',
        'field_format': 'bdsim2d',
        'field_file': '2dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'nearest',
        'expected_samples': 5000,
        'expected_x_sigma': 0.003369653184855525,
        'expected_xp_sigma': 0.0017251188722198771,
        'expected_x_mean': 3.561127408404445e-05,
        'expected_y_sigma': 0.021222188911178278,
        'expected_yp_sigma': 0.021194220357673686,
        'expected_y_mean': -0.00018837782195287218,
        'expected_energy_sigma': 0.019783911947098665,
        'expected_energy_mean': 0.9995454897522926,
    },
    {
        'name': '2d_cubic',
        'field_type': 'bmap2d',
        'field_format': 'bdsim2d',
        'field_file': '2dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'cubic',
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
        'name': '3d_linear',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_interpolator': 'linear',
        'expected_samples': 5000,
        'expected_x_sigma': 0.014324369577942097,
        'expected_xp_sigma': 0.007515015279365734,
        'expected_x_mean': 0.5660426369667053,
        'expected_y_sigma': 0.0044860079849747825,
        'expected_yp_sigma': 0.002270031791075037,
        'expected_y_mean': 0.22305434429049492,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': '3d_linearmag',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_interpolator': 'linearmag',
        'expected_samples': 5000,
        'expected_x_sigma': 0.014327077885869297,
        'expected_xp_sigma': 0.007516079445364859,
        'expected_x_mean': 0.5661105772137642,
        'expected_y_sigma': 0.004486432469495695,
        'expected_yp_sigma': 0.002270183744965747,
        'expected_y_mean': 0.22307491642832755,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': '3d_nearest',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_interpolator': 'nearest',
        'expected_samples': 5000,
        'expected_x_sigma': 0.015180510433518032,
        'expected_xp_sigma': 0.007908333993363166,
        'expected_x_mean': 0.5767620910525322,
        'expected_y_sigma': 0.00595327484532122,
        'expected_yp_sigma': 0.003382907443600101,
        'expected_y_mean': 0.22346383558809757,
        'expected_energy_sigma': 0.19783912688382108,
        'expected_energy_mean': 9.995454892158508,
    },
    {
        'name': '3d_cubic',
        'field_type': 'bmap3d',
        'field_format': 'bdsim3d',
        'field_file': '3dexample.dat',
        'field_length': '1.0',
        'beam_energy': '10.0',
        'field_interpolator': 'cubic',
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
    # 4D field - test 'linear', 'linearmag', 'nearest', 'cubic'
    {
        'name': '4d_linear',
        'field_type': 'bmap4d',
        'field_format': 'bdsim4d',
        'field_file': '4dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linear',
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
    {
        'name': '4d_linearmag',
        'field_type': 'bmap4d',
        'field_format': 'bdsim4d',
        'field_file': '4dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'linearmag',
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
    {
        'name': '4d_nearest',
        'field_type': 'bmap4d',
        'field_format': 'bdsim4d',
        'field_file': '4dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'nearest',
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
    {
        'name': '4d_cubic',
        'field_type': 'bmap4d',
        'field_format': 'bdsim4d',
        'field_file': '4dexample.dat',
        'field_length': '1.0',
        'beam_energy': '1.0',
        'field_interpolator': 'cubic',
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

FIELD_INTERPOLATION_CONFIGS = {
    'name': ['1d_linear', '1d_linearmag', '1d_nearest',
             '2d_linear', '2d_linearmag', '2d_nearest', '2d_cubic',
             '3d_linear', '3d_linearmag', '3d_nearest', '3d_cubic',
             '4d_linear', '4d_linearmag', '4d_nearest', '4d_cubic'],
    'FIELD_TYPE': ['bmap1d', 'bmap1d', 'bmap1d',
                   'bmap2d', 'bmap2d', 'bmap2d', 'bmap2d',
                   'bmap3d', 'bmap3d', 'bmap3d', 'bmap3d',
                   'bmap4d', 'bmap4d', 'bmap4d', 'bmap4d'],
    'FIELD_FORMAT': ['bdsim1d', 'bdsim1d', 'bdsim1d',
                     'bdsim2d', 'bdsim2d','bdsim2d', 'bdsim2d',
                     'bdsim3d', 'bdsim3d', 'bdsim3d', 'bdsim3d',
                     'bdsim4d', 'bdsim4d', 'bdsim4d', 'bdsim4d'],
    'FIELD_FILE': ['1dexample.dat', '1dexample.dat', '1dexample.dat',
                   '2dexample.dat', '2dexample.dat', '2dexample.dat', '2dexample.dat',
                   '3dexample.dat', '3dexample.dat', '3dexample.dat', '3dexample.dat',
                   '4dexample.dat', '4dexample.dat', '4dexample.dat', '4dexample.dat'],
    'FIELD_LENGTH': ['1.0', '1.0', '1.0',
                     '1.0', '1.0', '1.0', '1.0',
                     '1.0', '1.0', '1.0', '1.0',
                     '1.0', '1.0', '1.0', '1.0'],
    'BEAM_ENERGY': ['1.0', '1.0', '1.0',
                    '1.0', '1.0', '1.0', '1.0',
                    '10.0', '10.0', '10.0', '10.0',
                    '1.0', '1.0', '1.0', '1.0'],
    'FIELD_INTERPOLATOR': ['linear', 'linearmag', 'nearest',
                           'linear', 'linearmag', 'nearest', 'cubic',
                           'linear', 'linearmag', 'nearest', 'cubic',
                           'linear', 'linearmag', 'nearest', 'cubic'],
    'expected_samples': [5000, 5000, 5000,
                         5000, 5000, 5000, 5000,
                         5000, 5000, 5000, 5000,
                         5000, 5000, 5000, 5000],
 'expected_x_sigma': [0.005354202580646126, 0.005360657975868763, 0.0018571098686965446,
                      0.015804539022049573, 0.01969824900808272, 0.003369653184855525, 0.018419807674184555,
                      0.014324369577942097, 0.014327077885869297, 0.015180510433518032, 0.014271013981376675,
                      0.0017761945113131978, 0.0017761945113131978, 0.0017761945113131978, 0.0017761945113131978],
 'expected_xp_sigma': [0.001999076865966775, 0.0020020719479500963, 0.0004469812469452971,
                       0.010316564235008652, 0.013199531233735034, 0.0017251188722198771, 0.012275748730702854,
                       0.007515015279365734, 0.007516079445364859, 0.007908333993363166, 0.007416345178360901,
                       0.00035136689967708905, 0.00035136689967708905, 0.00035136689967708905, 0.00035136689967708905],
 'expected_x_mean': [0.009670783941715035, 0.009683479112656596, 0.009439205919555388,
                     8.499772373153859e-05, 0.00010222111806351676, 3.561127408404445e-05, 0.00010576244183662311,
                     0.5660426369667053, 0.5661105772137642, 0.5767620910525322, 0.5799057118058205,
                     1.3215959309269466e-05, 1.3215959309269466e-05, 1.3215959309269466e-05, 1.3215959309269466e-05],
 'expected_y_sigma': [0.013084708676590241, 0.01307943167235869, 0.01307691510397638,
                      0.0021579719729259575, 0.0038697662117538945, 0.021222188911178278, 0.0028903558879167703,
                      0.0044860079849747825, 0.004486432469495695, 0.00595327484532122, 0.004409168245179408,
                      0.0017857934386989697, 0.0017857934386989697, 0.0017857934386989697, 0.0017857934386989697],
 'expected_yp_sigma': [0.007202701938688855, 0.0072005010788836685, 0.007198931645154899,
                       0.0023403134007770198, 0.004129952429779972, 0.021194220357673686, 0.0032919311481443892,
                       0.002270031791075037, 0.002270183744965747, 0.003382907443600101, 0.002207146467824476,
                       0.0003543900315124393, 0.0003543900315124393, 0.0003543900315124393, 0.0003543900315124393],
 'expected_y_mean': [-0.49687208319306375, -0.49720062289237976, -0.49778312349915504,
                     1.8904866732009396e-05, 4.843972614912673e-06, -0.00018837782195287218, 1.0439962574434957e-05,
                     0.22305434429049492, 0.22307491642832755, 0.22346383558809757, 0.226743698656559,
                     -2.1371774069270798e-05, -2.1371774069270798e-05, -2.1371774069270798e-05, -2.1371774069270798e-05]
}

config_df = pd.DataFrame(FIELD_INTERPOLATION_CONFIGS)
configs = config_df.to_dict(orient="records")
param_keys = ["FIELD_TYPE", "FIELD_FORMAT", "FIELD_FILE", "FIELD_LENGTH", "BEAM_ENERGY", "FIELD_INTERPOLATOR"]

simulation_results = []

@pytest.fixture(scope="module", autouse=True)
def write_results_csv():
    yield

    pd.DataFrame(simulation_results).to_csv(
        "field_interpolation_results.csv",
        index=False,
    )

@pytest.mark.parametrize("config", configs, ids=[c['name'] for c in configs])
def test_field_interpolation(config):

    os.chdir(os.path.dirname(__file__))

    base_name = f"field_{config['name']}"
    template_name = ("field_map_interpolator.tpl")
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
