import pytest
import pybdsim
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "placement_sampler"
    template_name = base_name+".tpl"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"
    optics_name   = base_name+"_optics.root"
    combine_name  = base_name+"_combine.root"
    
    params = {
        'LENGTH': '1.0',
        'BEAM_ENERGY' : '1'
    }

    ngenerate = 500

    pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,params)
    pybdsim.Run.Bdsim(gmad_name,base_name,ngenerate,1)
    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    sampler_data = data.get_sampler("s1.")

    # extract the relevant data
    sampler_number = len(sampler_data)
    sampler_x = sampler_data['x']
    sampler_xp = sampler_data['xp']
    sampler_y = sampler_data['y']
    sampler_yp = sampler_data['yp']
    sampler_energy = sampler_data['energy']

    # calculate values that should remain consistent with every run
    sampler_x_sigma = sampler_x.std()
    sampler_x_mean = sampler_x.mean()
    sampler_xp_sigma = sampler_xp.std()
    sampler_xp_mean = sampler_xp.mean()

    sampler_y_sigma = sampler_y.std()
    sampler_y_mean = sampler_y.mean()
    sampler_yp_sigma = sampler_yp.std()
    sampler_yp_mean = sampler_yp.mean()

    sampler_energy_sigma = sampler_energy.std()
    sampler_energy_mean = sampler_energy.mean()

    # print the values used to validate the run
    print("assert sampler_number == ", sampler_number)

    print("assert pytest.approx(sampler_x_sigma, rel=1e-4) == ", sampler_x_sigma)
    print("assert pytest.approx(sampler_x_mean, rel=1e-4) == ", sampler_x_mean)
    print("assert pytest.approx(sampler_xp_sigma, rel=1e-4) == ", sampler_xp_sigma)
    print("assert pytest.approx(sampler_xp_mean, rel=1e-4) == ", sampler_xp_mean)

    print("assert pytest.approx(sampler_y_sigma, rel=1e-4) == ", sampler_y_sigma)
    print("assert pytest.approx(sampler_y_mean, rel=1e-4) == ", sampler_y_mean)
    print("assert pytest.approx(sampler_yp_sigma, rel=1e-4) == ", sampler_yp_sigma)
    print("assert pytest.approx(sampler_yp_mean, rel=1e-4) == ", sampler_yp_mean)

    print("assert pytest.approx(sampler_energy_sigma, rel=1e-4) == ", sampler_energy_sigma)
    print("assert pytest.approx(sampler_energy_mean, rel=1e-4) == ", sampler_energy_mean)

    assert sampler_number == 6

    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 0.04245527509351293
    assert pytest.approx(sampler_x_mean, rel=1e-4) == -0.01258978076900045
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.410648573003683
    assert pytest.approx(sampler_xp_mean, rel=1e-4) == 0.34474750608205795

    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.06614192545287227
    assert pytest.approx(sampler_y_mean, rel=1e-4) == -0.006785242973516385
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.34335081175825494
    assert pytest.approx(sampler_yp_mean, rel=1e-4) == 0.07147048320621252

    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.004263307653161516
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.003288653703445258
