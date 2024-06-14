"""
This script was written by Christoph Fruehlinger. It parameterizes the total EEG
power spectrum created in the Preprocessing script using the 'fitting oscilla-
tions & one over f' (FOOOF) toolbox (Donoghue et al., 2020).

Reference:
Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, Noto T, 
Lara AH, Wallis JD, Knight RT, Shestyuk A, & Voytek B (2020). Parameterizing 
neural power spectra into periodic and aperiodic components. Nature 
Neuroscience, 23, 1655-1665. DOI: 10.1038/s41593-020-00744-x
"""
# Import relevant packages. If you want to show plots in pane and save as fig, 
# import matplotlib packages below and uncomment all 'fm.plot'commands in 
# script.

import importlib
import subprocess

# List of packages to check and install if necessary
packages = ['os', 'multiprocessing', 'time']

for package in packages:
    try:
        # Try importing the package
        importlib.import_module(package)
        print(f'{package} is already installed.')
    except ImportError:
        print(f'{package} is not installed. Installing...')
        
        # Use pip to install the package
        subprocess.check_call(['pip', 'install', package])
        
        print(f'{package} has been installed.')

import time
import os
from multiprocessing.pool import ThreadPool 
#path = '/home/bay2875/Scripts_DataRequests/MVPA_Personality/'
path = './'
os.environ['PATH'] += ':'+path
from fooof_func import fooof_func


start_time = time.time()

if __name__ == '__main__':
    pool = ThreadPool()
    inputs = [0, 1, 2, 3]
    output = pool.map(fooof_func, inputs)

print('\nDone')

end_time = time.time()

execution_time = end_time - start_time
hours = int(execution_time // 3600)  # Convert seconds to hours
minutes = int((execution_time % 3600) // 60)  # Convert to minutes
seconds = int(execution_time % 60)


print(f"\nScript executed in {hours} hours, {minutes} minutes and {seconds} " 
      "seconds")