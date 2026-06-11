import pytest
import pybdsim
import numpy as np
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "gausstwiss"
    template_name = base_name+".tpl"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"
    optics_name   = base_name+"_optics.root"

    l  = 2.0 
    data = {
        'BETX': '4',
        'BETY' : '4',
        'EMITX' : '5e-7',
        'EMITY' : '5e-7'
    }
    ngenerate = 5000

    pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)
    pybdsim.Run.Bdsim(gmad_name,base_name,ngenerate,1)

    data = pybdsim.DataPandas.BDSIMOutput("root_name.root")
    primary_data = d.get_primary()
    
    primary_x = primary_data['x']
    primary_xp = primary_data['xp']
    primary_y = primary_data['y']
    primary_yp = primary_data['yp']  
    primary_t = primary_data['T']
    primary_energy = primary_data['energy']
    n=len(x)

    cov_xxp = np.cov(x,xp)
    cov_yyp = np.cov(y,yp)
    cov_et = np.cov(y,yp)


    sigma_x_generated = np.sqrt(cov_xxp[0][0])
    sigma_xp_generated = np.sqrt(cov_xxp[1][1])
    sigma_y_generated = np.sqrt(cov_yyp[0][0])
    sigma_yp_generated = np.sqrt(cov_yyp[1][1])
    sigma_x_calculated = np.sqrt(float(data['BETX'])*float(data['EMITX']))
    sigma_y_calculated = np.sqrt(float(data['BETY'])*float(data['EMITY']))

    
    emittance_x_generated = np.sqrt(cov_xxp[0][0]*cov_xxp[1][1]-cov_xxp[0][1]*cov_xxp[1][0])
    emittance_y_generated = np.sqrt(cov_yyp[0][0]*cov_yyp[1][1]-cov_yyp[0][1]*cov_yyp[1][0])
    emittance_x_input = float(data['EMITX'])
    emittance_y_input = float(data['EMITY'])


    print("assert sigma_x_calculated = ",sigma_x_calculated, "= sigma_x_generated = ", sigma_x_generated)
    print("assert sigma_y_calculated = ",sigma_y_calculated, "= sigma_y_generated = ", sigma_y_generated)
    print("assert emittance_x_input = ", emittance_x_input, "= emittance_x_calculated = ", emittance_x_generated)
    print("assert emittance_y_input = ", emittance_y_input, "= emittance_y_calculated = ", emittance_y_generated)


    #assert(n = ngenerate)




    



