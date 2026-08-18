import pytest
import pybdsim
import os
import uproot
import numpy as np
import scipy.constants as const


def getSamplerData(f, samplerName, variable, mask=None) :
    d =  f[samplerName+variable].array(library='np') 
    elem = d[0]
    if np.ndim(elem) == 0:
        d_out=np.array(d)
    else:
        d_out =  np.array([x[0] for x in d if x.size>0])
    if mask is not None:
        d_out = d_out[mask]
    d_out_mean=np.mean(d_out)
    
    return d_out_mean, d_out


def calculate_number_of_particles(f, samplerName, particle_ID=0):   
    _,parent_ID= getSamplerData(f, samplerName, 'parentID')
    mask=parent_ID==particle_ID
    n_out=np.sum(mask)
    
    return n_out,mask


def get_mean_projected_momentum(f, samplerName, pi='zp', mask=None): 
    _,p_out=getSamplerData(f, samplerName, 'p', mask)
    _,pi_frac=getSamplerData(f, samplerName, pi, mask)
    pi_out=np.multiply(pi_frac,p_out)
    pi_out_mean=np.mean(pi_out)
    
    return pi_out_mean
    
    
def get_theoretical_beam_parameters(Ek_in, p_type, sim_env):
    c=const.c
    e=const.e
    if p_type=="e-":
        m=const.electron_mass
    elif p_type=="proton":
        m=const.proton_mass
    gamma=Ek_in*e*1e9/(m*c**2)+1
    v=c*np.sqrt(1-1/gamma**2)
    time_out_th=sim_env["l"]/v*1e9 #s
    E_m=m*c**2
    E_tot=Ek_in*e*1e9+E_m 
    pz_in=np.sqrt(E_tot**2-E_m**2)/c
    pz_in=pz_in*c/(e*1e9) # GeV/c
    
    return E_tot, pz_in, time_out_th
    
    
@pytest.fixture
def sim_env():
    os.chdir(os.path.dirname(__file__))
    base_name = "drift2"
    return {
        "base_name": base_name,
        "template_name": base_name + ".tpl",
        "gmad_name": base_name + ".gmad",
        "root_name": base_name + ".root",
        "n_in": 2500,
        "phi_in": 0,
        "theta_in": 0,
        "l": 2 #m      
    }
    
    
class TestClass:
    @pytest.mark.parametrize(
    "Ek_in, p_type",
    [
        (1e-3, "e-"),
        (1e3, "e-"), #GeV
        (1e-4, "proton"),
        (1.0, "proton"),
    ]
    )
    def test_drift2(self, sim_env, Ek_in, p_type) :   
        E_tot, pz_in, time_out_th= get_theoretical_beam_parameters(Ek_in,p_type,sim_env)
        data = {
            'LENGTH': str(sim_env["l"]),
            'BEAM_ENERGY' : str(Ek_in), #GeV
            'P_TYPE': p_type
        }      
        pybdsim.Run.RenderGmadJinjaTemplate(sim_env["template_name"],sim_env["gmad_name"],data)
        pybdsim.Run.Bdsim(sim_env["gmad_name"],sim_env["base_name"],sim_env["n_in"],1)
        with uproot.open(sim_env["root_name"]) as file:
            samplerName="d1."
            f=file['Event'][samplerName]
            n_out,mask=calculate_number_of_particles(f,samplerName,particle_ID=0) #only primaries
            theta_out_mean,_=getSamplerData(f,samplerName,"theta",mask)
            pz_out_mean=get_mean_projected_momentum(f,samplerName,pi='zp',mask=mask)
            phi_out_mean,_=getSamplerData(f,samplerName,"phi",mask)
            Ek_out_mean,_=getSamplerData(f,samplerName,"kineticEnergy",mask)
            t_out_mean,_=getSamplerData(f,samplerName,"T",mask)
            s_sampler_mean,_=getSamplerData(f,samplerName,"S",mask)

        assert Ek_out_mean==pytest.approx(Ek_in,rel=1e-9)
        assert phi_out_mean==pytest.approx(sim_env["phi_in"], rel=1e-9)
        assert theta_out_mean==pytest.approx(sim_env["theta_in"],rel=1e-9)
        assert sim_env["n_in"]==n_out
        assert pz_out_mean==pytest.approx(pz_in,1e-5)
        assert t_out_mean==pytest.approx(time_out_th,1e-6)
        assert s_sampler_mean==pytest.approx(sim_env["l"],1e-9)
