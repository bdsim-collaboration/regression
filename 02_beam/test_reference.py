import pytest
import pybdsim
import numpy as np
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "reference"
    template_name = base_name+".tpl"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"

    l  = 2.0 
    params = {
        'ENERGY': '1',
        'LENGTH': '2'
    }
    ngenerate = 5000

    pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,params)
    pybdsim.Run.Bdsim(gmad_name,base_name,ngenerate,1)

    data = pybdsim.DataPandas.BDSIMOutput(root_name)
    primary_data = data.get_primary()
    sampler_data=data.get_sampler('t1.')    

    primary_x = primary_data['x']
    primary_xp = primary_data['xp']
    primary_y = primary_data['y']
    primary_yp = primary_data['yp']  
    primary_t = primary_data['T']
    primary_energy = primary_data['energy']
    number_particles=len(primary_x)

    cov_xxp = np.cov(primary_x,primary_xp)
    cov_yyp = np.cov(primary_y,primary_yp)

    sigma_x_generated = np.sqrt(cov_xxp[0][0])
    sigma_xp_generated = np.sqrt(cov_xxp[1][1])
    sigma_y_generated = np.sqrt(cov_yyp[0][0])
    sigma_yp_generated = np.sqrt(cov_yyp[1][1])
    
    emittance_x_generated = np.sqrt(cov_xxp[0][0]*cov_xxp[1][1]-cov_xxp[0][1]*cov_xxp[1][0])
    emittance_y_generated = np.sqrt(cov_yyp[0][0]*cov_yyp[1][1]-cov_yyp[0][1]*cov_yyp[1][0])

    sampler_number = len(sampler_data)
    sampler_x = sampler_data['x']
    sampler_xp = sampler_data['xp']
    sampler_y = sampler_data['y']
    sampler_yp = sampler_data['yp']
    sampler_energy = sampler_data['energy']

    sampler_cov_xxp = np.cov(sampler_x,sampler_xp)
    sampler_cov_yyp = np.cov(sampler_y,sampler_yp)

    sigma_x_sampler = np.sqrt(sampler_cov_xxp[0][0])
    sigma_xp_sampler = np.sqrt(sampler_cov_xxp[1][1])
    sigma_y_sampler = np.sqrt(sampler_cov_yyp[0][0])
    sigma_yp_sampler = np.sqrt(sampler_cov_yyp[1][1])

    emittance_x_sampler = np.sqrt(sampler_cov_xxp[0][0]*sampler_cov_xxp[1][1]-sampler_cov_xxp[0][1]*sampler_cov_xxp[1][0])
    emittance_y_sampler = np.sqrt(sampler_cov_yyp[0][0]*sampler_cov_yyp[1][1]-sampler_cov_yyp[0][1]*sampler_cov_yyp[1][0])

    assert(number_particles == ngenerate)
    assert(0 == pytest.approx(sigma_x_generated,abs=1e-3))
    assert(0 == pytest.approx(sigma_x_generated,abs=1e-3))
    assert(0 == pytest.approx(emittance_x_generated,abs=1e-3))
    assert(0 == pytest.approx(emittance_y_generated,abs=1e-3))
    assert(0 == pytest.approx(sigma_xp_generated,abs=1e-3))
    assert(0 == pytest.approx(sigma_yp_generated,abs=1e-3))

    assert(sampler_number == ngenerate)
    assert(0 == pytest.approx(sigma_x_sampler,abs=1e-3))
    assert(0 == pytest.approx(sigma_x_sampler,abs=1e-3))
    assert(0 == pytest.approx(emittance_x_sampler,abs=1e-3))
    assert(0 == pytest.approx(emittance_y_sampler,abs=1e-3))
    assert(0 == pytest.approx(sigma_xp_sampler,abs=1e-3))
    assert(0 == pytest.approx(sigma_yp_sampler,abs=1e-3))






    



