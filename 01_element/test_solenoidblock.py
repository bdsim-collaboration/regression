import pytest
import pybdsim
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from jinja2 import Environment, FileSystemLoader



def check_sampler_point_mapping(points_sampler, points_validated, tolerance=1e-4):


    assert len(points_sampler) == len(points_validated), \
        f"Point count mismatch: {len(points_sampler)} vs {len(points_validated)}"
    
    
    distances = cdist(points_sampler, points_validated) 

    row_ind, col_ind = linear_sum_assignment(distances)

    matched_distances = distances[row_ind, col_ind]
    mean_distance = matched_distances.mean()
    max_distance = matched_distances.max()

    
    assert mean_distance < tolerance, \
        f"Mean distance {max_distance:.6f} exceeds threshold {tolerance}"
    
    return 
def test() :    
    os.path.dirname(__file__)
    
    base_name     = "solenoidBlock"
    template_name = base_name+".tpl"
    gmad_name     = base_name+".gmad"
    root_name     = base_name+".root"
    beam_name     = base_name+"Beam.dat"
    sampler_data_name = base_name+"Sampler.dat"
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template(template_name)
    data = {}
    blocks = {
        'beam': beam_name,
        'sampler': sampler_data_name,
        'gmad': gmad_name
    }
    for block_name, output_file in blocks.items():
        block_content = ''.join(template.blocks[block_name](template.new_context(data)))
    
        with open(output_file, 'w') as f:
            f.write(block_content)
    pybdsim.Run.Bdsim(gmad_name,base_name,21)
    d = pybdsim.Data.Load(root_name)

    sampler = pybdsim.Data.SamplerData(d, 'd1.')


    samplerX = sampler.data['x']
    samplerY = sampler.data['y']

    df_sampler = pd.DataFrame({'X': samplerX, 'Y': samplerY})
    df_validated = pd.read_csv('solenoidBlockSampler.dat', sep=' ') # Previously validated sampler
    points_sampler = df_sampler[['X', 'Y']].values
    points_validated = df_validated[['X', 'Y']].values

    check_sampler_point_mapping(points_sampler, points_validated, tolerance=1e-4)
if __name__ == "__main__":
    test()
