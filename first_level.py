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
             'seed': '/data/mridata/jbrown/brains/rois/{seed_name}'}
selectfiles = Node(SelectFiles(templates), name="selectfiles")

# For distributing seed and nuisance mask paths
# Merge Node to combine seed mask and nuisance mask paths
seed_plus_nuisance = Node(utilMerge(2), name = 'seed_plus_nuisance')
seed_plus_nuisance.inputs.in2 = nuisance_masks

# For distributing seed and nuisance mask paths
distributor = Node(IdentityInterface(fields=['rois']), name='distributor')

# 1. Obtain timeseries for seed and nuisance variables
# 1a. Merge all 3D functional images into a single 4D image
merge = Node(Merge(dimension = 't',
                   output_type = 'NIFTI',
                   tr = 2.0), name = 'merge')

# 1b. Take mean of all voxels in each roi at each timepoint
ts = MapNode(ImageMeants(), name = 'ts', iterfield = ['mask'])

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
    (seed_plus_nuisance, distributor, [('out', 'rois')]),
    (distributor, ts, [('rois', 'mask')]),
    (ts, datasink, [('out_file', 'timeseries')])
                ])

# Visualize the workflow and run it
get_timeseries.write_graph(graph2use='flat')
get_timeseries.run()