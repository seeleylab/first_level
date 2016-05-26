# Import necessary modules
from nipype.interfaces.spm import Level1Design
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

# Specify variables
subjdir = []
seeds = []

#What type of Node to use for inputting subject paths and how to specify
#that subjects have different paths (but all fitting a pattern:

#What type of Node to use for looping over seeds? Iterable or MapNode
#IdentityInterface.iterables(subjectID) > DataGrabber

#Infosource
infosource = Node(IdentityInterface(fields=['subject_path', 'seed_path']), name="infosource")
infosource.iterables = [('subject_path', subjdir), ('seed_path', seeds)]

# SelectFiles
templates = {'func': '{subject_path}/processedfmri_TRCNnSFmDI/swua_filteredf*.nii',
             'seed': '{seed_path}'}
selectfiles = Node(SelectFiles(templates), name="selectfiles")