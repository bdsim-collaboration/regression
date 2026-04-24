import pytest
import pybdsim
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from jinja2 import Environment, FileSystemLoader


def parse_sampler_locations(gmad_file):
    sampler_locations = {}
    with open(gmad_file, 'r') as f:
        for line in f:
            match = re.match(r'^\s*(s\d+):\s*samplerplacement.*s=([-+]?\d+\.?\d*[eE]?[-+]?\d*)\*(m|mm)', line)
            if match:
                sampler_name = match.group(1)
                s_position = float(match.group(2))
                sampler_locations[sampler_name] = s_position
    return sampler_locations


def check_sampler_point_mapping(df, tolerance=5e-4, coords=['x', 'y']):
    for cellz, group in df.groupby('cellZ'):
        samplers = sorted(group['name'].unique())
        
        point_counts = group.groupby('name').size().to_dict()
        unique_counts = set(point_counts.values())
        
        assert len(unique_counts) == 1, \
            f"cellZ={cellz}: Point count mismatch: {point_counts}"
        
        sampler_data = {
            name: group[group['name'] == name][coords].values 
            for name in samplers
        }
        
        for i, sampler1 in enumerate(samplers):
            for sampler2 in samplers[i+1:]:
                points1 = sampler_data[sampler1]
                points2 = sampler_data[sampler2]
                
                distances = cdist(points1, points2, metric='euclidean')
                row_ind, col_ind = linear_sum_assignment(distances)
                matched_distances = distances[row_ind, col_ind]
                max_distance = matched_distances.max()
                
                assert max_distance <= tolerance, \
                    f"cellZ={cellz}: {sampler1} ↔ {sampler2} max distance " \
                    f"{max_distance:.2e} exceeds tolerance {tolerance:.2e}"
    
    return True


def test():    
    base_name = "muoncooler"
    template_name = base_name + ".tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"
    beam_name = base_name + "Beam.dat"
    
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template(template_name)
    data = {}
    blocks = {
        'beam': beam_name,
        'gmad': gmad_name
    }
    
    for block_name, output_file in blocks.items():
        block_content = ''.join(template.blocks[block_name](template.new_context(data)))
        with open(output_file, 'w') as f:
            f.write(block_content)
    
    pybdsim.Run.Bdsim(gmad_name, base_name, 3)
    #Sampler data is stored as z =0 and s = -1, so we need to parse the gmad to get the actual s positions of the samplers
    sampler_locations = parse_sampler_locations(gmad_name)
    
    d = pybdsim.Data.Load(root_name)
    rows = []
    for sampler_name in sampler_locations.keys():
        sampler = pybdsim.Data.SamplerData(d, sampler_name)
        x_data = sampler.data['x']
        y_data = sampler.data['y']
        s_position = sampler_locations[sampler_name]
        
        for j in range(len(x_data)):
            rows.append({
                'name': sampler_name,
                'x': x_data[j],
                'y': y_data[j],
                's': s_position
            })
    
    df = pd.DataFrame(rows)
    df['cellZ'] = df['s'] % 1600
    
    #Find the 1-to-1 mapping between the samplers at the same cellZ and check that the maximum distance between matched points is within the tolerance
    check_sampler_point_mapping(df, tolerance=5e-4)
    

if __name__ == "__main__":
    test()