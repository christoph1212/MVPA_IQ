
# -*- coding: utf-8 -*-
"""
This script defines the function that is used in 'FOOOF_python.py' to para-
meterize the total EEG spectrum into its periodic and aperiodic parts.
Written by Christoph Fruehlinger, August 2023. Change Basefolder according to 
your directory.

Input: preprocessed EEG data

Output: parameterized EEG data into its periodic and aperiodic components
"""
#%% Import relevant Packages
import importlib
import subprocess

# List of packages to check and install if necessary
packages = ['numpy', 'pandas', 'os', 'scipy', 'fooof', 'multiprocessing', 'matplotlib']

for package in packages:
    try:
        # Try importing the package
        importlib.import_module(package)       

    except ImportError:     
        # Use pip to install the package
        subprocess.check_call(['pip', 'install', package])

import numpy as np
import pandas as pd
import os
from scipy.io import savemat
from fooof import FOOOF
#import matplotlib # in case results should appear in pane


#%% Folder Setup                                
Basefolder = "../DataValidation/"
# Add Preprocdir
preproc_directory = os.path.join(Basefolder, 'Preprocessed')
fooof_directory = os.path.join(Basefolder, 'FOOOF_Data')
os.makedirs(fooof_directory, exist_ok=True)


# Set paths and create folders if necessary
Pre_EO_Path = os.path.join(preproc_directory, 'Pre', 'EyesOpen')

Pre_EO_Plots = os.path.join(fooof_directory, 'Plots', 'Pre', 'EyesOpen')
os.makedirs(Pre_EO_Plots, exist_ok=True)

Pre_EO_summary = os.path.join(fooof_directory, 'Model_Summary', 'Pre', 
                                 'EyesOpen')
os.makedirs(Pre_EO_summary, exist_ok=True)

Pre_EO_periodic = os.path.join(fooof_directory, 'Periodic', 'Pre', 
                                  'EyesOpen')
os.makedirs(Pre_EO_periodic, exist_ok=True)

Pre_EO_aperiodic = os.path.join(fooof_directory, 'Aperiodic', 'Pre', 
                                   'EyesOpen')
os.makedirs(Pre_EO_aperiodic, exist_ok=True)

Pre_EO_aperiodic_param = os.path.join(fooof_directory, 'Aperiodic', 'Parameter', 'Pre', 
                                   'EyesOpen')
os.makedirs(Pre_EO_aperiodic_param, exist_ok=True)

##
Post_EO_Path = os.path.join(preproc_directory, 'Post', 'EyesOpen')

Post_EO_Plots = os.path.join(fooof_directory, 'Plots', 'Post', 'EyesOpen')
os.makedirs(Post_EO_Plots, exist_ok=True)

Post_EO_summary = os.path.join(fooof_directory, 'Model_Summary', 'Post', 
                                 'EyesOpen')
os.makedirs(Post_EO_summary, exist_ok=True)

Post_EO_periodic = os.path.join(fooof_directory, 'Periodic', 'Post', 
                                  'EyesOpen')
os.makedirs(Post_EO_periodic, exist_ok=True)

Post_EO_aperiodic = os.path.join(fooof_directory, 'Aperiodic', 'Post', 
                                   'EyesOpen')
os.makedirs(Post_EO_aperiodic, exist_ok=True)

Post_EO_aperiodic_param = os.path.join(fooof_directory, 'Aperiodic', 'Parameter', 'Post', 
                                   'EyesOpen')
os.makedirs(Post_EO_aperiodic_param, exist_ok=True)

input_paths = [Pre_EO_Path, Post_EO_Path]
plot_paths = [Pre_EO_Plots, Post_EO_Plots]
summary_paths = [Pre_EO_summary, Post_EO_summary]
periodic_path = [Pre_EO_periodic, Post_EO_periodic]
aperiodic_path = [Pre_EO_aperiodic, Post_EO_aperiodic]
parameter_path = [Pre_EO_aperiodic_param, Post_EO_aperiodic_param]
conditionnames = ['pre_Average_EO', 'post_Average_EO']

#%% FOOOF
def fooof_func(condition):
    # List of files and participant IDs
    files = os.listdir(input_paths[condition])
    Participants = [file.split('_')[1] for file in 
                    os.listdir(input_paths[condition])]
    model_summary = {}                                                             
    
    # Loop through participants
    for run, Participant in enumerate(Participants):
	# Test if file exist
        file_path1 = os.path.join(periodic_path[condition], f'rest_{Participant}_{conditionnames[condition]}_periodic.csv')
        file_path2 = os.path.join(aperiodic_path[condition], f'rest_{Participant}_{conditionnames[condition]}_aperiodic.csv') 
        file_path3 = os.path.join(parameter_path[condition], f'rest_{Participant}_{conditionnames[condition]}_aperiodic_parameter.csv')
           
        if (os.path.isfile(file_path1) and os.path.isfile(file_path2) and os.path.isfile(file_path3)):
            continue

        print(f'Condition {condition+1}/2. Running Participant {run+1}/{len(Participants)}')

        # Load files and transform into numpy array (easier to handle)

        data = pd.read_csv(os.path.join(input_paths[condition], files[run]), 
                           header=None)
        data_new = data.to_numpy()
        # Create empty df for periodic and aperiodic and add EEG channels
        periodic = pd.DataFrame(np.zeros(np.shape(data_new)))                      
        periodic.iloc[:,0] = data.iloc[:,0]                                        

        aperiodic = pd.DataFrame(np.zeros(np.shape(data_new)))                     
        aperiodic.iloc[:,0] = data.iloc[:,0]                                       

        aperiodic_params = pd.DataFrame(np.zeros((len(data_new), 3)))             
        aperiodic_params.iloc[:,0] = data.iloc[:,0]


        # Set frequency range for toolbox
        freqs = np.arange(0.5, 30.5, 0.5)
        # Loop through each channel in file

        for i in range(0, len(data_new)):
            psd = np.squeeze(data_new[i,1:])
            # Settings for toolbox and initiate parameterization
            fm = FOOOF(aperiodic_mode = 'fixed', peak_width_limits = [1, 12])
            fm.fit(freqs, psd, [0.5, 30])
            #fm.print_results() # uncomment if results should appear in pane
            # Get power and save as df
            periodic.loc[i, 1:] = fm._spectrum_flat 
            aperiodic.loc[i, 1:] = fm._ap_fit  
            aperiodic_params.loc[i,1:] = fm.aperiodic_params_
            
        fooof_results = fm.get_results()                                           
        fooof_results_dict = fooof_results._asdict()                               
        model_summary[Participant + '_modelsum'] = fooof_results_dict              

        # Model Summary
        savemat(os.path.join(summary_paths[condition], Participant +               
                conditionnames[condition] +'_modelsum.mat'), 
                fooof_results_dict)
        # Aperiodic
        file_path = os.path.join(aperiodic_path[condition],                        
                    Participant + '_' + conditionnames[condition] +'_aperiodic.csv')
        aperiodic.to_csv(file_path, index=False, header=False)                     
        
        file_path = os.path.join(parameter_path[condition], 
                    Participant + '_' + conditionnames[condition] +'_aperiodic_parameter.csv')
        aperiodic_params.to_csv(file_path, index=False, header=False)
        # Periodic
        file_path = os.path.join(periodic_path[condition], 
                    Participant + '_' + conditionnames[condition] +'_periodic.csv')      
        periodic.to_csv(file_path, index=False, header=False)                      
# %%