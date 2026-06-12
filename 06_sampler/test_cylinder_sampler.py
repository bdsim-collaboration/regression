import pytest
import pybdsim
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "cylinder_sampler"
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
    sampler_data = data.get_samplerc("s1.")

    sampler_number = len(sampler_data)
    sampler_z = sampler_data['z']
    sampler_zp = sampler_data['zp']
    sampler_rp = sampler_data['rp']
    sampler_phi = sampler_data['phi']
    sampler_phip = sampler_data['phip']
    sampler_energy = sampler_data['totalEnergy']

    sampler_z_sigma = sampler_z.std()
    sampler_z_mean = sampler_z.mean()
    sampler_zp_sigma = sampler_zp.std()
    sampler_zp_mean = sampler_zp.mean()

    sampler_rp_sigma = sampler_rp.std()
    sampler_rp_mean = sampler_rp.mean()

    sampler_phi_sigma = sampler_phi.std()
    sampler_phi_mean = sampler_phi.mean()
    sampler_phip_sigma = sampler_phip.std()
    sampler_phip_mean = sampler_phip.mean()

    sampler_energy_sigma = sampler_energy.std()
    sampler_energy_mean = sampler_energy.mean()

    print("assert sampler_number == ", sampler_number)

    print("assert pytest.approx(sampler_z_sigma, rel=1e-4) == ", sampler_z_sigma)
    print("assert pytest.approx(sampler_z_mean, rel=1e-4) == ", sampler_z_mean)
    print("assert pytest.approx(sampler_zp_sigma, rel=1e-4) == ", sampler_zp_sigma)
    print("assert pytest.approx(sampler_zp_mean, rel=1e-4) == ", sampler_zp_mean)

    print("assert pytest.approx(sampler_rp_sigma, rel=1e-4) == ", sampler_rp_sigma)
    print("assert pytest.approx(sampler_rp_mean, rel=1e-4) == ", sampler_rp_mean)

    print("assert pytest.approx(sampler_phi_sigma, rel=1e-4) == ", sampler_phi_sigma)
    print("assert pytest.approx(sampler_phi_mean, rel=1e-4) == ", sampler_phi_mean)
    print("assert pytest.approx(sampler_phip_sigma, rel=1e-4) == ", sampler_phip_sigma)
    print("assert pytest.approx(sampler_phip_mean, rel=1e-4) == ", sampler_phip_mean)

    print("assert pytest.approx(sampler_energy_sigma, rel=1e-4) == ", sampler_energy_sigma)
    print("assert pytest.approx(sampler_energy_mean, rel=1e-4) == ", sampler_energy_mean)

    assert sampler_number == 7

    assert pytest.approx(sampler_z_sigma, rel=1e-4) == 0.14109729326665602
    assert pytest.approx(sampler_z_mean, rel=1e-4) == -0.3276261846934046
    assert pytest.approx(sampler_zp_sigma, rel=1e-4) == 0.46784390793909336
    assert pytest.approx(sampler_zp_mean, rel=1e-4) == 0.24842584771769388

    assert pytest.approx(sampler_rp_sigma, rel=1e-4) == 0.28632609711225854
    assert pytest.approx(sampler_rp_mean, rel=1e-4) == 0.7943912936108453

    assert pytest.approx(sampler_phi_sigma, rel=1e-4) == 2.2100343813521146
    assert pytest.approx(sampler_phi_mean, rel=1e-4) == -0.15624664723873138
    assert pytest.approx(sampler_phip_sigma, rel=1e-4) == 0.2794790812054533
    assert pytest.approx(sampler_phip_mean, rel=1e-4) == -0.05538616941443512

    assert pytest.approx(sampler_energy_sigma, rel=1e-4) == 0.0005348490219988756
    assert pytest.approx(sampler_energy_mean, rel=1e-4) == 0.0005639208663654115

