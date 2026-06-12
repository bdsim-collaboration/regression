import pytest
import pybdsim
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "element_sampler"
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

    # run and load the pybdsim data
    pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,params)
    pybdsim.Run.Bdsim(gmad_name,base_name,ngenerate,1)
    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    sampler_data = data.get_sampler("d2.")

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

    assert pytest.approx(sampler_x_sigma, rel=1e-4) == 1.5877514529773704
    assert pytest.approx(sampler_x_mean, rel=1e-4) == -0.8595296616355578
    assert pytest.approx(sampler_xp_sigma, rel=1e-4) == 0.20538206969024522
    assert pytest.approx(sampler_xp_mean, rel=1e-4) == -0.19051110558211803

    assert pytest.approx(sampler_y_sigma, rel=1e-4) == 0.6645647357818376
    assert pytest.approx(sampler_y_mean, rel=1e-4) == 0.2016523089259863
    assert pytest.approx(sampler_yp_sigma, rel=1e-4) == 0.41745956967563386
    assert pytest.approx(sampler_yp_mean, rel=1e-4) == -0.10976425806681316

    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 1.77073810789213e-05
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.00013629104554032287

    # these assertions will throw an error if there is a bug or change in the code
    # no error will show if the assertions are successful

    # copy and paste the printed out values if there is a change in the code/update
    # (only when confirmed it is not a bug)
