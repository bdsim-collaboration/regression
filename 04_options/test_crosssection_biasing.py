import pytest
import pybdsim
import os

def test() :

    os.chdir(os.path.dirname(__file__))
    
    base_name     = "crosssection_biasing"
    template_name = base_name+".tpl"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"

    data = {'BIAS_FACTOR': '10'}

    pybdsim.Run.RenderGmadJinjaTemplate(template_name,gmad_name,data)
    pybdsim.Run.Bdsim(gmad_name,base_name,ngenerate=100,seed=1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)
    s = d.get_sampler("sampler.")
    n_muons = len(s['x'])
    weight_sum = sum(s['weight'])

    # count varies with G4 version but unbiased muon count is ~12
    assert(20 < n_muons < 50)
    assert(6 < weight_sum < 18)