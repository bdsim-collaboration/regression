import pytest
import pybdsim
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "spherical_sampler"
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
    sampler_data = data.get_samplers("s1.")

    sampler_number = len(sampler_data)
    sampler_theta = sampler_data['theta']
    sampler_thetap = sampler_data['thetap']
    sampler_phi = sampler_data['phi']
    sampler_phip = sampler_data['phip']
    #sampler_rp = sampler_data['rp']
    sampler_energy = sampler_data['totalEnergy']

    sampler_theta_sigma = sampler_theta.std()
    sampler_theta_mean = sampler_theta.mean()
    sampler_thetap_sigma = sampler_thetap.std()
    sampler_thetap_mean = sampler_thetap.mean()

    sampler_phi_sigma = sampler_phi.std()
    sampler_phi_mean = sampler_phi.mean()
    sampler_phip_sigma = sampler_phip.std()
    sampler_phip_mean = sampler_phip.mean()

    sampler_energy_sigma = sampler_energy.std()
    sampler_energy_mean = sampler_energy.mean()

    print("assert sampler_number == ", sampler_number)

    print("assert pytest.approx(sampler_theta_sigma, rel=1e-4) == ", sampler_theta_sigma)
    print("assert pytest.approx(sampler_theta_mean, rel=1e-4) == ", sampler_theta_mean)
    print("assert pytest.approx(sampler_thetap_sigma, rel=1e-4) == ", sampler_thetap_sigma)
    print("assert pytest.approx(sampler_thetap_mean, rel=1e-4) == ", sampler_thetap_mean)

    print("assert pytest.approx(sampler_phi_sigma, rel=1e-4) == ", sampler_phi_sigma)
    print("assert pytest.approx(sampler_phi_mean, rel=1e-4) == ", sampler_phi_mean)
    print("assert pytest.approx(sampler_phip_sigma, rel=1e-4) == ", sampler_phip_sigma)
    print("assert pytest.approx(sampler_phip_mean, rel=1e-4) == ", sampler_phip_mean)

    print("assert pytest.approx(sampler_energy_sigma, rel=1e-4) == ", sampler_energy_sigma)
    print("assert pytest.approx(sampler_energy_mean, rel=1e-4) == ", sampler_energy_mean)

    assert sampler_number == 2119

    assert pytest.approx(sampler_theta_sigma, rel=1e-4) == 0.07266751873464505
    assert pytest.approx(sampler_theta_mean, rel=1e-4) == 3.078928793303637
    assert pytest.approx(sampler_thetap_sigma, rel=1e-4) == 1.0886758064426896
    assert pytest.approx(sampler_thetap_mean, rel=1e-4) == 1.369520605494738

    assert pytest.approx(sampler_phi_sigma, rel=1e-4) == 1.846997304969656
    assert pytest.approx(sampler_phi_mean, rel=1e-4) == 0.0529459856373973
    assert pytest.approx(sampler_phip_sigma, rel=1e-4) == 1.2538820617519535
    assert pytest.approx(sampler_phip_mean, rel=1e-4) == 0.04129473779128511

    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.42443211593789076
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.23627107728795502
