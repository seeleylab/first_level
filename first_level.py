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
TR = 2.0

experiment_dir = '/data/mridata/jdeng/tools/first_level/nipype'
output_dir = 'get_timeseries_output'
working_dir = 'get_timeseries_workingdir'

## CREATE NODES
# Infosource
infosource = Node(IdentityInterface(fields=['subject_path', 'seed_path']), name="infosource")
infosource.iterables = [('subject_path', subjdir), ('seed_path', seeds)]

# SelectFiles
templates = {'func': '/data/mridata/jdeng/tools/first_level/{subject_path}/rsfmri/processedfmri_TRCNnSFmDI/images/swua_filteredf*.nii',
             'seed': '/data/mridata/jbrown/brains/rois/{seed_path}'}
selectfiles = Node(SelectFiles(templates), name="selectfiles")

# 1. Obtain timeseries for seed and nuisance variables
# 1a. Merge all 3D functional images into a single 4D image
# in_files <-- selectfiles func
# out_files --> get_timeseries
merge = Node(Merge(dimension = 't',
                   output_type = 'NIFTI',
                   tr = 2.0), name = 'merge')

# 1b. Take mean of all voxels in mask at each timepoint
# in_files <-- merge, selectfiles seed
# out_files --> Datasink
ts = Node(ImageMeants(), name = 'ts')

## CREATE WORKFLOW
# Create a short workflow to get the timeseries for seed + nuisance variables for each subject
get_timeseries = Workflow(name='get_timeseries')
get_timeseries.base_dir = opj(experiment_dir, working_dir)

# Datasink
datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir), name="datasink")

# Connect all components of the workflow
get_timeseries.connect([
    (infosource, selectfiles, [('subject_path', 'subject_path'), ('seed_path', 'seed_path')]),
    (selectfiles, merge, [('func', 'in_files')]),
    (merge, ts, [('merged_file', 'in_file')]),
    (selectfiles, ts, [('seed', 'mask')]),
    (ts, datasink, [('out_file', 'timeseries')])
                ])

# Visualize the workflow and run it
get_timeseries.write_graph(graph2use='flat')
get_timeseries.run()