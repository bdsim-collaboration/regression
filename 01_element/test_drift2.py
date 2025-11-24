import pytest
import pybdsim
import os
import uproot
import numpy as np
import scipy.constants as const
from scipy import stats

def getSamplerData(f,samplerName, variable, mask) :
    d =  f[samplerName+variable].array(library='np') 
    elem = d[0]
    if np.ndim(elem) == 0:
        d_out=np.array(d)
    else:
        d_out =  np.array([x[0] for x in d if x.size>0])
    d_out = d_out[mask]
    
    return d_out

def get_p_value(mu0,x, rtol=1e-6,atol=1e-12):
    mean = np.mean(x)
    sigma = np.std(x, ddof=1)
    test_value=0
    if abs(mean) > atol:
        if np.abs(sigma/mean) <rtol:
            test_value=1
    else:
        if sigma< atol:
            test_value=1
            
    if test_value:
        if np.isclose(mean, mu0, rtol):
            return 0.0   # H0 perfect
        else:
            return 1.0   # H0 clearly wrong
            
    t_stat, p_value = stats.ttest_1samp(x, mu0)

    return p_value
    
def get_theoretical_beam_parameters(ek_in,p_type,sim_env):
    c=const.c
    e=const.e
    if p_type=="e-":
        m=const.electron_mass
    elif p_type=="proton":
        m=const.proton_mass
    gamma=ek_in*e*1e9/(m*c**2)+1
    v=c*np.sqrt(1-1/gamma**2)
    time_out_th=sim_env["l"]/v*1e9 #s
    e_m=m*c**2
    E_tot=ek_in*e*1e9+e_m 
    pz_rel=np.sqrt(E_tot**2-e_m**2)/c
    pz_rel=pz_rel*c/(e*1e9) # GeV/c
    
    return E_tot, pz_rel, time_out_th
    
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

    alpha_test=0.05
    
    @pytest.mark.parametrize(
    "ek_in, p_type",
    [
        (1e-3, "e-"),
        (1e3, "e-"), #GeV
        (1e-4, "proton"),
        (1.0, "proton"),
    ]
    )
    
    def test_drift2(self, sim_env, ek_in, p_type) :

        E_tot, pz_rel, time_out_th= get_theoretical_beam_parameters(ek_in,p_type,sim_env)
        data = {
            'LENGTH': str(sim_env["l"]),
            'BEAM_ENERGY' : str(ek_in), #GeV
            'P_TYPE': p_type
        }
        
        pybdsim.Run.RenderGmadJinjaTemplate(sim_env["template_name"],sim_env["gmad_name"],data)
        pybdsim.Run.Bdsim(sim_env["gmad_name"],sim_env["base_name"],sim_env["n_in"],1)

        with uproot.open(sim_env["root_name"]) as file:
            samplerName="d1."
            f=file['Event'][samplerName]
            a=f[samplerName+'parentID'].array(library='np')
            parent_ID= np.array([x[0] for x in a if len(x)>0])
            mask=parent_ID==0 #from initial beam only
            n_out=np.sum(mask)
            
            theta_out=getSamplerData(f,samplerName,"theta", mask)
            p_theta=get_p_value(sim_env["theta_in"],theta_out)
            
            a=f[samplerName+'p'].array(library='np')
            p_out= np.array([x[0] for x in a if len(x)>0])
            a=f[samplerName+'zp'].array(library='np')
            pz_frac= np.array([x[0] for x in a if len(x)>0])
            pz_out=np.multiply(pz_frac,p_out)
            pz_out=pz_out[mask]
            p_pz=get_p_value(pz_rel,pz_out)
    
            phi_out=getSamplerData(f,samplerName,"phi", mask)
            p_phi=get_p_value(sim_env["phi_in"],phi_out)
            
            ek_out=getSamplerData(f,samplerName,"kineticEnergy", mask)
            p_ek=get_p_value(ek_in,ek_out)
            
            t_out=getSamplerData(f,samplerName,"T", mask)
            p_t=get_p_value(time_out_th,t_out)
            
            s_sampler=getSamplerData(f,samplerName,"S", mask)
            p_l=get_p_value(sim_env["l"],s_sampler)

        assert ( p_ek< self.alpha_test)
        assert (p_phi< self.alpha_test)
        assert (p_theta< self.alpha_test)
        assert (sim_env["n_in"]==n_out)
        assert (p_pz< self.alpha_test)
        assert ( p_t< self.alpha_test)
        assert ( p_l< self.alpha_test)

