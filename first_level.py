# Import necessary modules
import sys
from os.path import join as opj
from nipype.interfaces.fsl import Merge, ImageMeants
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces.spm import Level1Design, EstimateModel, EstimateContrast
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.utility import Merge as utilMerge
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

# Specify variables
#subjdir = raw_input('Input the absolute paths to your subjects: ').split()
subjdir = sys.argv[1].split()
subjdir = [path.strip('/') for path in subjdir]
#seed_paths = raw_input('Input the absolute paths to the seeds: ').split()
seed_paths = sys.argv[2].split()
nuisance_masks = ['/data/mridata/SeeleyToolbox/SeeleyFirstLevel/proc/csf_ant_post_bilateral.nii',
                  '/data/mridata/SeeleyToolbox/SeeleyFirstLevel/proc/avg152T1_white_mask.nii']
TR = 2.0

experiment_dir = '/data/mridata/jdeng/tools/first_level/nipype'
output_dir = 'first_level_output'

## CREATE NODES
# For distributing subject paths
infosource = Node(IdentityInterface(fields=['subject_path', 'seed']),
                  name="infosource")
infosource.iterables = [('subject_path', subjdir), ('seed', seed_paths)]

info = dict(func = [['subject_path', 'rsfmri/processedfmri_TRCNnSFmDI/images/swua_filteredf*.nii']],
            motion = [['subject_path', 'rsfmri/processedfmri_TRCNnSFmDI/motion_params_filtered.txt']])

selectfiles = Node(DataGrabber(infields = ['subject_path'],
                              outfields = ['func', 'motion'],
                              base_directory = '/',
                              template = '%s/%s',
                              template_args = info,
                              sort_filelist = True),
                  name = 'selectfiles')

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

# 1c. - Merge nuisance ts with motion parameters to create nuisance_regressors.txt.
#     - Take temporal derivatives of nuisance_regressors.txt and append to nuisance_regressors.txt
#       to create nuisance_regressors_tempderiv.txt
#     - Square nuisance_regressors_tempderiv.txt and append to nuisance_regressors_tempderiv.txt,
#       then append seed timeseries in front to create seed_nuisance_regressors.txt
def make_regressors_files(regressors_ts_list, mot_params):
    import numpy as np
    import os
    num_timepoints = 235    # change this to not be hard-coded
    num_nuisance = len(regressors_ts_list) - 1
    
    # make nuisance_regressors.txt
    nr = np.zeros((num_timepoints, num_nuisance))
        
    i = 0
    for ts in regressors_ts_list[1:]:
        nr[:,i] = np.loadtxt(ts)[:]
        i += 1
        
    nr = np.hstack((nr, np.loadtxt(mot_params)[:]))
    np.savetxt(os.path.join(os.getcwd(), 'nuisance_regressors.txt'), nr, fmt='%.7e')
    # make nuisance_regressors_tempderiv.txt
    td = np.gradient(nr, axis=0)
    nr_td = np.hstack((nr, td))
    np.savetxt(os.path.join(os.getcwd(), 'nuisance_regressors_tempderiv.txt'), nr_td, fmt='%.7e')
    # make seed_nuisance_regressors.txt
    seed_ts = np.loadtxt(regressors_ts_list[0])[:]
    sq = np.square(nr_td)
    snr = np.hstack((seed_ts[:, np.newaxis], nr_td, sq))
    np.savetxt(os.path.join(os.getcwd(), 'seed_nuisance_regressors.txt'), snr, fmt='%.7e')
    
    # return nuisance_regressors.txt, nuisance_regressors_tempderiv.txt, and seed_nuisance_regressors.txt
    return os.path.join(os.getcwd(), 'nuisance_regressors.txt'), os.path.join(os.getcwd(), 'nuisance_regressors_tempderiv.txt'), os.path.join(os.getcwd(), 'seed_nuisance_regressors.txt')

make_regressors_files = Node(Function(input_names = ['regressors_ts_list', 'mot_params'],
                                      output_names = ['nr', 'nr_td', 'snr'],
                                      function = make_regressors_files),
                                name = 'make_regressors_files')
# 2. Build statistical model
# 2a. Create SPM.mat design matrix
def model_helper(regressors_file):
    from nipype.interfaces.base import Bunch
    import numpy as np
    regressors_file_data = np.loadtxt(regressors_file).T.tolist()
    num_regressors = len(regressors_file_data)
    condition_vector = ['seed'] + ['nuisance']*(num_regressors-1) + ['constant']
    weights_vector = [float(1)] + [float(0)]*(num_regressors-1) + [float(0)]
    contrasts = [('Condition1', 'T', condition_vector, weights_vector)]
    subject_info = [(Bunch(regressors = regressors_file_data,
                           regressor_names = ['seed'] + ['nuisance']*(num_regressors-1) + ['constant']))]
    return subject_info, contrasts

model_helper = Node(Function(input_names = ['regressors_file'],
                             output_names = ['subject_info', 'contrasts'],
                             function = model_helper),
                        name = 'model_helper')

session_info = Node(SpecifyModel(high_pass_filter_cutoff = 128,
                                    input_units = 'secs',
                                    time_repetition = 2.0),
                    name = 'session_info')

model_spec = Node(Level1Design(timing_units = 'secs',
                              interscan_interval = 2.0,
                              microtime_resolution = 16,
                              microtime_onset = 1,
                              bases = {'hrf':{'derivs': [0,0]}},
                              global_intensity_normalization = 'none',
                              mask_threshold = 0.8,
                              model_serial_correlations = 'AR(1)',
                              volterra_expansion_order = 2),
                 name = 'model_spec')

est_model = Node(EstimateModel(estimation_method = {'Classical': 1}),
                 name = 'est_model')

est_con = Node(EstimateContrast(),
               name = 'est_con')

## CREATE WORKFLOW
# Create a workflow to return the seed nuisance regressors and seed map(s) for a subject
first_level = Workflow(name='first_level')
first_level.base_dir = experiment_dir

# Datasink
datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir), name="datasink")

substitutions = [('_subject_name_', '_'), ('_seed_name_', '')]
datasink.inputs.substitutions = substitutions

# Helper functions for connections
def makelist(item):
    return [item]

# Connect all components of the workflow
first_level.connect([
    (infosource, selectfiles, [('subject_path', 'subject_path')]),
    (selectfiles, merge, [('func', 'in_files')]),
    (merge, ts, [('merged_file', 'in_file')]),
    (infosource, seed_plus_nuisance, [('seed', 'in1')]),
    (seed_plus_nuisance, ts, [('out', 'mask')]),
    (ts, make_regressors_files, [('out_file', 'regressors_ts_list')]),
    (selectfiles, make_regressors_files, [('motion', 'mot_params')]),
    (make_regressors_files, datasink, [('nr', 'timeseries'),
        ('nr_td', 'timeseries.@nr_td'),
        ('snr', 'timeseries.@snr')]),
    (make_regressors_files, model_helper, [('snr', 'regressors_file')]),
    (selectfiles, session_info, [(('func', makelist), 'functional_runs')]),
    (model_helper, session_info, [('subject_info', 'subject_info')]),
    (session_info, model_spec, [('session_info', 'session_info')]),
    (model_spec, est_model, [('spm_mat_file', 'spm_mat_file')]),
    (est_model, est_con, [('beta_images', 'beta_images'),
        ('residual_image', 'residual_image'),
        ('spm_mat_file', 'spm_mat_file')]),
    (model_helper, est_con, [('contrasts', 'contrasts')]),
    (est_con, datasink, [('spm_mat_file', 'job_stats'),
        ('con_images', 'images'),
        ('spmT_images', 'images.@T')])
                ])

# Visualize the workflow and run it
first_level.write_graph(graph2use='flat')
first_level.run(plugin='SGE', plugin_args=dict(template='/data/mridata/jdeng/tools/grid/q.sh'))