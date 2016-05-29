# Import necessary modules
from os.path import join as opj
from nipype.interfaces.fsl import Merge, ImageMeants
from nipype.interfaces.spm import Level1Design
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

# Specify variables
subjdir = ['IC005-1']
seeds = ['rPCC_sphere.nii',
         '44rvFI_vAt12_C2vx123.nii',
         'Boxer_DMidbrainTeg_sphere_3-5_-15_-8.nii']
nuisance_masks = ['csf_ant_post_bilateral.nii',
                  'avg152T1_white_mask.nii']
TR = 2.0

experiment_dir = '/data/mridata/jdeng/tools/first_level/nipype'
output_dir = 'get_timeseries_output'
#working_dir = 'get_timeseries_workingdir'

## CREATE NODES
# Infosource
infosource = Node(IdentityInterface(fields=['subject_name', 'seed_name', 'nuisance_name']),
                  name="infosource")
infosource.iterables = [('subject_name', subjdir), ('seed_name', seeds), ('nuisance_name', nuisance_masks)]

# SelectFiles
templates = {'func': '/data/mridata/jdeng/tools/first_level/{subject_name}/rsfmri/processedfmri_TRCNnSFmDI/images/swua_filteredf*.nii',
             'seed': '/data/mridata/jbrown/brains/rois/{seed_name}',
             'nuisance': '/data/mridata/SeeleyToolbox/SeeleyFirstLevel/proc/{nuisance_name}'}
selectfiles = Node(SelectFiles(templates), name="selectfiles")

# Join seeds and nuisance masks for use with ImageMeants MapNode later
def join_rois_func(list1, list2):
    return list1 + list2
    
join_rois = Node(Function(input_names = ['list1', 'list2'],
                          output_names = ["rois"],
                          function = join_rois_func),
                 name = 'join_rois')

# 1. Obtain timeseries for seed and nuisance variables
# 1a. Merge all 3D functional images into a single 4D image
merge = Node(Merge(dimension = 't',
                   output_type = 'NIFTI',
                   tr = 2.0), name = 'merge')

# 1b. Take mean of all voxels in each roi at each timepoint
ts = MapNode(ImageMeants(), name = 'ts', iterfield = ['mask'])

# # 1c. Take mean of all voxels in each nuisance mask at each timepoint
# ts_nuisance = Node(ImageMeants(), name = 'ts_nuisance')
# 
# # 1b. Take mean of all voxels in each seed mask at each timepoint
# ts_seed = Node(ImageMeants(), name = 'ts_seed')

## CREATE WORKFLOW
# Create a short workflow to get the timeseries for seed + nuisance variables for each subject
get_timeseries = Workflow(name='get_timeseries')
#get_timeseries.base_dir = opj(experiment_dir, working_dir)
get_timeseries.base_dir = experiment_dir

# Datasink
datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir), name="datasink")

# Connect all components of the workflow
get_timeseries.connect([
    (infosource, selectfiles, [('subject_name', 'subject_name'), ('seed_name', 'seed_name'), ('nuisance_name', 'nuisance_name')]),
    (selectfiles, merge, [('func', 'in_files')]),
    (merge, ts, [('merged_file', 'in_file')]),
    (selectfiles, join_rois, [('seed', 'list1'), ('nuisance', 'list2')]),
    (join_rois, ts, [('rois', 'mask')]),
    #(merge, ts_nuisance, [('merged_file', 'in_file')]),
    #(selectfiles, ts_nuisance, [('nuisance', 'mask')]),
    (ts, datasink, [('out_file', 'seed_timeseries')]),
    #(ts_nuisance, datasink, [('out_file', 'nuisance_timeseries')])
                ])

# Visualize the workflow and run it
get_timeseries.write_graph(graph2use='flat')
get_timeseries.run()