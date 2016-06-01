# Import necessary modules
from os.path import join as opj
from nipype.interfaces.fsl import Merge, ImageMeants
from nipype.interfaces.spm import Level1Design
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.utility import Merge as utilMerge
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

# Specify variables
subjdir = ['IC005-1']
seeds = ['rPCC_sphere.nii',
         '44rvFI_vAt12_C2vx123.nii',
         'Boxer_DMidbrainTeg_sphere_3-5_-15_-8.nii']
nuisance_masks = ['/data/mridata/SeeleyToolbox/SeeleyFirstLevel/proc/csf_ant_post_bilateral.nii',
                  '/data/mridata/SeeleyToolbox/SeeleyFirstLevel/proc/avg152T1_white_mask.nii']
TR = 2.0

experiment_dir = '/data/mridata/jdeng/tools/first_level/nipype'
output_dir = 'get_timeseries_output'

## CREATE NODES
# For distributing subject paths
infosource = Node(IdentityInterface(fields=['subject_name', 'seed_name']),
                  name="infosource")
infosource.iterables = [('subject_name', subjdir), ('seed_name', seeds)]

templates = {'func': '/data/mridata/jdeng/tools/first_level/{subject_name}/rsfmri/processedfmri_TRCNnSFmDI/images/swua_filteredf*.nii',
             'seed': '/data/mridata/jbrown/brains/rois/{seed_name}',
             'motion': '/data/mridata/jdeng/tools/first_level/{subject_name}/rsfmri/interfmri_TRCNnSFmDI/motion_corr/rp*.txt'}
selectfiles = Node(SelectFiles(templates), name="selectfiles")

# For merging seed and nuisance mask paths and then distributing them downstream
seed_plus_nuisance = Node(utilMerge(2), name = 'seed_plus_nuisance')
seed_plus_nuisance.inputs.in2 = nuisance_masks

# 1. Obtain timeseries for seed and nuisance variables
# 1a. Merge all 3D functional images into a single 4D image
merge = Node(Merge(dimension = 't',
                   output_type = 'NIFTI',
                   tr = 2.0), name = 'merge')

# 1b. Take mean of all voxels in each roi at each timepoint
ts = MapNode(ImageMeants(), name = 'ts', iterfield = ['mask'])

# 2. Merge nuisance ts
def make_nuisance_regressors(regressors_ts_list, mot_params):
    import numpy as np
    import os
    num_timepoints = 235    # change this to not be hard-coded
    num_regressors = len(regressors_ts_list) - 1
    nuisance_regressors = np.zeros((num_timepoints, num_regressors))
    
    i = 0
    for ts in regressors_ts_list[1:]:
        nuisance_regressors[:,i] = np.loadtxt(ts)[:]
        i += 1
    
    nuisance_regressors = np.hstack((nuisance_regressors, np.loadtxt(mot_params)[:]))
    
    np.savetxt(os.path.join(os.getcwd(), 'nuisance_regressors.txt'), nuisance_regressors)
    return(os.path.join(os.getcwd(), 'nuisance_regressors.txt'))

make_nuisance_regressors = Node(Function(input_names = ['regressors_ts_list', 'mot_params'],
                                         output_names = ['out_file'],
                                         function = make_nuisance_regressors),
                                name = 'make_nuisance_regressors')

## CREATE WORKFLOW
# Create a short workflow to get the timeseries for seed + nuisance variables for each subject
get_timeseries = Workflow(name='get_timeseries')
get_timeseries.base_dir = experiment_dir

# Datasink
datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir), name="datasink")

substitutions = [('_subject_name_', '_'), ('_seed_name_', '')]
datasink.inputs.substitutions = substitutions

# Connect all components of the workflow
get_timeseries.connect([
    (infosource, selectfiles, [('subject_name', 'subject_name'), ('seed_name', 'seed_name')]),
    (selectfiles, merge, [('func', 'in_files')]),
    (merge, ts, [('merged_file', 'in_file')]),
    (selectfiles, seed_plus_nuisance, [('seed', 'in1')]),
    (seed_plus_nuisance, ts, [('out', 'mask')]),
    (ts, make_nuisance_regressors, [('out_file', 'regressors_ts_list')]),
    (selectfiles, make_nuisance_regressors, [('motion', 'mot_params')]),
    (ts, datasink, [('out_file', 'timeseries')]),
    (make_nuisance_regressors, datasink, [('out_file', 'nuisance_regressors')]),
                ])

# Visualize the workflow and run it
get_timeseries.write_graph(graph2use='flat')
get_timeseries.run('MultiProc', plugin_args={'n_procs': 16})