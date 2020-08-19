
import os
import warnings
import copy
import numpy as np

warnings.filterwarnings('ignore')

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
from ipywidgets import widgets, interact, interactive
from IPython.display import HTML, display, Javascript
import json, sys, shutil, glob, fileupload

#Global paramters for input and output directory. These paramteres need to be changed if you want to run EPIC on local machine instead of Dockers
input_dir = '/home/kuan-hao/Documents/EPIC_Self/input/real_data/trimmed_data'
projects =  ['']
for f in [f.split(os.sep)[-2] for f in glob.glob(input_dir+"*")]:
    print(f)
    if not f.endswith("out") and not f.endswith("fa_files"):
        projects.append(f)
def f(**kwargs):
    return None

directoryName_i = widgets.SelectMultiple(
    options=projects,
    value=[projects[0]],
    description='Input',
    disabled=False
)


features_i = interactive(f, MI=False, Bayes=False, Euclidean=True, WCC=False, Jaccard=False, PCCN=False, PCC=True, Apex=False)
num_cores_i = interactive(f, num_cores="1")
clf_i = widgets.RadioButtons(
    options=["Random forest", "SVM"],
    description='Classifier',
    disabled=False
)
target_species_i = interactive(f, target_species="taxid i.e. 6239 (Worm)")
mode_i = widgets.RadioButtons(
    options=['exp', 'comb'],
    description='Mode',
    disabled=False
)


fa_source_i = widgets.RadioButtons(
    options=['GM', 'STRING', 'FILE'],
    description='FA source',
    disabled=False
)

def _handle_upload(change):
    w = change['owner']
    with open("/tmp/" + w.filename, 'wb') as f:
        f.write(w.data)
    print('Uploaded to `{}` ({:.2f} kB)'.format(
        "/tmp/" + w.filename, len(w.data) / 2**10))
    
fa_dir = input_dir+os.sep+"fa_files"

fa_files = ["No File"]

if os.path.exists(fa_dir):
    for f in [f.split(os.sep)[-1] for f in glob.glob(fa_dir + os.sep + "*")]:
            fa_files.append(f)

            
fa_file = widgets.SelectMultiple(
    options=fa_files,
    value=[fa_files[0]],
    description='FA File select',
    disabled=False
)


    
ref_file = fileupload.FileUploadWidget()
ref_file.observe(_handle_upload, names='data')
if np.array_equal(projects, ['/home/kuan-hao/Documents/EPIC_Self/EPIC/test_data/elution_profiles']):
    print '              !!!!! ERROR: No directories found !!!!\nPlease create a directory that contains all elution profile files'
else: 
    display(directoryName_i)
