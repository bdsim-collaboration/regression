import pytest
import pybdsim
import os
import uproot
import numpy as np
import scipy.constants as const

class TestClass:
    
    def test_LowEnergy_electron(self) :
        os.chdir(os.path.dirname(__file__))
        base_name     = "drift2"
        template_name = base_name+".tpl"
        gmad_name     = base_name+".gmad"
        root_name     = base_name+".root"

        l  = 2.0 #m
        ek_in=1e-3 #GeV
        phi_in=0
        theta_in=0
        n_in=2500
        p_type="e-"
        c=const.c
        e=const.e
        m=const.electron_mass
        e_m=m*c**2
        E_tot=ek_in*e*1e9+e_m 
        pz_rel=np.sqrt(E_tot**2-e_m**2)/c
        pz_rel=pz_rel*c/(e*1e9) # GeV/c
        data = {
            'LENGTH': str(l),
            'BEAM_ENERGY' : str(ek_in), #GeV
            'P_TYPE': p_type
        }
        
        pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)
        pybdsim.Run.Bdsim(gmad_name,base_name,n_in,1)

        with uproot.open(root_name) as file:
            name="d1."
            f=file['Event'][name]
            a=f[name+'parentID'].array(library='np')
            parent_ID= np.array([x[0] for x in a if len(x)>0])
            mask=parent_ID==0 #from initial beam only
            n_out=np.sum(mask)
            
            a=f[name+'theta'].array(library='np')
            theta_out= np.array([x[0] for x in a if x.size>0])
            theta_out=theta_out[mask]
            theta_out=np.mean(theta_out)
            
            a=f[name+'p'].array(library='np')
            p_out= np.array([x[0] for x in a if len(x)>0])
            a=f[name+'zp'].array(library='np')
            pz_frac= np.array([x[0] for x in a if len(x)>0])
            pz_out=np.multiply(pz_frac,p_out)
            pz_out=pz_out[mask]
            pz_out=np.mean(pz_out)

            
            a=f[name+'phi'].array(library='np')
            phi_out= np.array([x[0] for x in a if x.size>0])
            phi_out=phi_out[mask]
            phi_out=np.mean(phi_out)
            
            a=f[name+'kineticEnergy'].array(library='np') # GeV
            ek_out= np.array([x[0] for x in a if len(x)>0])
            ek_out=ek_out[mask]
            ek_out=np.mean(ek_out)  

        assert (pytest.approx(ek_out)==ek_in) & (pytest.approx(phi_out)==phi_in) & (pytest.approx(theta_out)==theta_in) & (n_in==n_out) & (pytest.approx(pz_out)==pz_rel)
       

        
    def test_HighEnergy_electron(self) :

        os.chdir(os.path.dirname(__file__))        
        base_name     = "drift2"
        template_name = base_name+".tpl"
        gmad_name     = base_name+".gmad"
        root_name     = base_name+".root"

        l  = 2.0 
        ek_in=1e3 #GeV
        phi_in=0
        theta_in=0
        n_in=2500
        p_type="e-"
        c=const.c
        e=const.e
        m=const.electron_mass   
        e_m=m*c**2
        E_tot=ek_in*e*1e9+ e_m
        pz_rel=np.sqrt(E_tot**2-e_m**2)/c
        pz_rel=pz_rel*c/(e*1e9)
        data = {
            'LENGTH': str(l),
            'BEAM_ENERGY' : str(ek_in), #GeV
            'P_TYPE': p_type
        }

        pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)
        pybdsim.Run.Bdsim(gmad_name,base_name,n_in,1)

        with uproot.open(root_name) as file:
            name="d1."
            f=file['Event'][name]
            a=f[name+'parentID'].array(library='np')
            parent_ID= np.array([x[0] for x in a if len(x)>0])
            mask=parent_ID==0 #from initial beam only
            n_out=np.sum(mask)
            
            a=f[name+'theta'].array(library='np')
            theta_out= np.array([x[0] for x in a if x.size>0])
            theta_out=theta_out[mask]
            theta_out=np.mean(theta_out)
            
            a=f[name+'phi'].array(library='np')
            phi_out= np.array([x[0] for x in a if x.size>0])
            phi_out=phi_out[mask]
            phi_out=np.mean(phi_out)
            
            a=f[name+'p'].array(library='np')
            p_out= np.array([x[0] for x in a if len(x)>0])
            a=file['Event'][name][name+'zp'].array(library='np')
            pz_frac= np.array([x[0] for x in a if len(x)>0])
            pz_out=np.multiply(pz_frac,p_out)
            pz_out=pz_out[mask]
            pz_out=np.mean(pz_out)
            
            a=f[name+'kineticEnergy'].array(library='np') # in GeV
            ek_out= np.array([x[0] for x in a if len(x)>0])
            ek_out=ek_out[mask]
            ek_out=np.mean(ek_out)  

        assert (pytest.approx(ek_out)==ek_in) & (pytest.approx(phi_out)==phi_in) & (pytest.approx(theta_out)==theta_in) & (n_in==n_out) & (pytest.approx(pz_out)==pz_rel)
        
        
    def test_LowEnergy_proton(self) :
        os.chdir(os.path.dirname(__file__))       
        base_name     = "drift2"
        template_name = base_name+".tpl"
        gmad_name     = base_name+".gmad"
        root_name     = base_name+".root"

        l  = 2.0 
        ek_in=1e-4 #GeV
        phi_in=0
        theta_in=0
        n_in=2500
        p_type="proton"
        c=const.c
        e=const.e
        m=const.proton_mass
        e_m=m*c**2
        E_tot=ek_in*e*1e9+ e_m
        pz_rel=np.sqrt(E_tot**2-e_m**2)/c
        pz_rel=pz_rel*c/(e*1e9) #GeV
        data = {
            'LENGTH': str(l),
            'BEAM_ENERGY' : str(ek_in), #GeV
            'P_TYPE': p_type
        }

        pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)      
        pybdsim.Run.Bdsim(gmad_name,base_name,n_in,1)

        with uproot.open(root_name) as file:
            name="d1."
            f=file['Event'][name]
            a=f[name+'parentID'].array(library='np')
            parent_ID= np.array([x[0] for x in a if len(x)>0])
            mask=parent_ID==0 #from initial beam only
            n_out=np.sum(mask)
            
            a=f[name+'theta'].array(library='np')
            theta_out= np.array([x[0] for x in a if x.size>0])
            theta_out=theta_out[mask]
            theta_out=np.mean(theta_out)
            
            a=f[name+'phi'].array(library='np')
            phi_out= np.array([x[0] for x in a if x.size>0])
            phi_out=phi_out[mask]
            phi_out=np.mean(phi_out)
            
            a=f[name+'p'].array(library='np')
            p_out= np.array([x[0] for x in a if len(x)>0])
            a=f[name+'zp'].array(library='np')
            pz_frac= np.array([x[0] for x in a if len(x)>0])
            pz_out=np.multiply(pz_frac,p_out)
            pz_out=pz_out[mask]
            pz_out=np.mean(pz_out)
            
            a=f[name+'kineticEnergy'].array(library='np') # in GeV
            ek_out= np.array([x[0] for x in a if len(x)>0])
            ek_out=ek_out[mask]
            ek_out=np.mean(ek_out)  

        assert (pytest.approx(ek_out)==ek_in) & (pytest.approx(phi_out)==phi_in) & (pytest.approx(theta_out)==theta_in) & (n_in==n_out) & (pytest.approx(pz_out)==pz_rel)
       

        
    def test_HighEnergy_proton(self) :

        os.chdir(os.path.dirname(__file__))        
        base_name     = "drift2"
        template_name = base_name+".tpl"
        gmad_name     = base_name+".gmad"
        root_name     = base_name+".root"

        l  = 2.0 
        ek_in=1 #GeV
        phi_in=0
        theta_in=0
        n_in=2500
        p_type="proton"
        c=const.c
        e=const.e
        m=const.proton_mass
        e_m=m*c**2
        E_tot=ek_in*e*1e9+ e_m
        pz_rel=np.sqrt(E_tot**2-e_m**2)/c
        pz_rel=pz_rel*c/(e*1e9)
        data = {
            'LENGTH': str(l),
            'BEAM_ENERGY' : str(ek_in), #GeV
            'P_TYPE': p_type
        }

        pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)
        
        pybdsim.Run.Bdsim(gmad_name,base_name,n_in,1)

        with uproot.open(root_name) as file:
            name="d1."
            f=file['Event'][name]
            a=f[name+'parentID'].array(library='np')
            parent_ID= np.array([x[0] for x in a if len(x)>0])
            mask=parent_ID==0 #from initial beam only
            n_out=np.sum(mask)
            
            a=f[name+'theta'].array(library='np')
            theta_out= np.array([x[0] for x in a if x.size>0])
            theta_out=theta_out[mask]
            theta_out=np.mean(theta_out)
            
            a=f[name+'phi'].array(library='np')
            phi_out= np.array([x[0] for x in a if x.size>0])
            phi_out=phi_out[mask]
            phi_out=np.mean(phi_out)
            
            a=f[name+'p'].array(library='np')
            p_out= np.array([x[0] for x in a if len(x)>0])
            a=f[name+'zp'].array(library='np')
            pz_frac= np.array([x[0] for x in a if len(x)>0])
            pz_out=np.multiply(pz_frac,p_out)
            pz_out=pz_out[mask]
            pz_out=np.mean(pz_out)
            
            a=f[name+'kineticEnergy'].array(library='np') # in GeV
            ek_out= np.array([x[0] for x in a if len(x)>0])
            ek_out=ek_out[mask]
            ek_out=np.mean(ek_out)  

        assert (pytest.approx(ek_out)==ek_in) & (pytest.approx(phi_out)==phi_in) & (pytest.approx(theta_out)==theta_in) & (n_in==n_out) & (pytest.approx(pz_out)==pz_rel)
        
        
