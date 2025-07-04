#!/bin/python3

"""
GOAL: Calculate a list of TOTAL secondary current powers generated by each dipole 
in a list of primary current dipoles.
"""
import numpy as np
import pandas as pd
import sys
# path to the folder containing the library duneuropy.so
duneuropy_path='/home/bic/vikramn/build-release/duneuro-py/src/'

# path to the npz-archive containing the volume conductor information 
input_dir = './'

sys.path.append(duneuropy_path)

import duneuropy as dp
from multiprocessing import Pool
def calculate_power(dipole_position, dipole_moment):
    dipole = dp.Dipole3d(dipole_position, dipole_moment)

    # compute potential
    #print('Numerically computing potentials')
    meeg_driver = dp.MEEGDriver3d(driver_cfg)
    solution_storage = meeg_driver.makeDomainFunction()
    meeg_driver.solveEEGForward(dipole, solution_storage, driver_cfg)
    dissipated_power = meeg_driver.computePower(solution_storage) # compute the integral of <sigma * grad u, grad u> over the head volume
    return dissipated_power

###################
# ADAPT THIS
###################


# load data
#subject_list=['sub-0002','sub-0004','sub-0006', 'sub-0007'] #EDIT HERE!!! 
subject_list=['sub-0004','sub-0006','sub-0007']
for subject in subject_list:
    # load data
    print('Loading volume conductor data for subject:', subject)
    volume_conductor_path=input_dir+subject+'_vc_3layer.npz'
    vc_data = np.load(volume_conductor_path)

    nodes = vc_data['nodes']
    elements = vc_data['elements']
    labels = vc_data['labels']
    conductivities = vc_data['conductivities']
    print('Data loaded')
    
    source_model = 'partial_integration'
    write_output = True
    
    # create driver
    mesh_cfg = {'nodes' : nodes, 'elements' : elements}
    tensor_cfg = {'labels' : labels, 'conductivities' : conductivities}
    volume_conductor_cfg = {'grid' : mesh_cfg, 'tensors' : tensor_cfg}
    driver_cfg = {'type' : 'fitted', 'solver_type' : 'cg', 'element_type' : 'tetrahedron', 'post_process' : 'true', 'subtract_mean' : 'true'}
    solver_cfg = {'reduction' : '1e-14', 'edge_norm_type' : 'houston', 'penalty' : '20', 'scheme' : 'sipg', 'weights' : 'tensorOnly'}
    driver_cfg['solver'] = solver_cfg
    driver_cfg['volume_conductor'] = volume_conductor_cfg
    
    #print('Creating driver')
    #print('Driver created')
    
    # setting source model
    # create config
    source_model_config_partial_integration = {'type' : 'partial_integration'}
    source_model_config_venant = \
    {
      'type' : 'multipolar_venant',
      'referenceLength' : 20,
      'weightingExponent' : 1,
      'relaxationFactor' : 1e-6,
      'restrict' : True,
      'initialization' : 'closest_vertex'
    }
    source_model_config_subtraction = \
    {
      'type' : 'subtraction',
      'intorderadd' : 0,
      'intorderadd_lb' : 0
    }
    source_model_config_localized_subtraction = \
    {
      'type' : 'localized_subtraction',
      'restrict' : False,
      'initialization' : 'single_element',
      'intorderadd_eeg_patch' : 0,
      'intorderadd_eeg_boundary' : 0,
      'intorderadd_eeg_transition' : 0,
      'extensions' : 'vertex vertex'
    }
    source_model_config_database = \
    {
      'partial_integration' : source_model_config_partial_integration,
      'subtraction' : source_model_config_subtraction,
      'multipolar_venant' : source_model_config_venant,
      'localized_subtraction' : source_model_config_localized_subtraction
    }
    
    driver_cfg['source_model'] = source_model_config_database[source_model]
    
    ###################
    # ADAPT THIS
    ###################
    # set up dipole. First position, then moment.
    
    #iterate through a list of dipole positions and dipole amplitudes
    dipole_positions=np.array(pd.read_csv(input_dir+subject+'_vertices.csv', header=None))
    #print(dipole_positions)
    #dipole_moments=np.array(pd.read_csv(input_dir+subject+'_idip_rms.csv', header=None))
    dipole_moments=np.array(pd.read_csv(input_dir+subject+'_idip_0s.csv',header=None))
    #print(dipole_moments)
        
    with Pool(10) as pool:
        powers = pool.starmap(calculate_power, zip(dipole_positions, dipole_moments))
        
    pd.DataFrame(np.array(powers)).to_csv(subject+'_fem_p_0s.csv', header=False, index=False)

