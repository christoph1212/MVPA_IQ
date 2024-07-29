# :construction: This Repository is not complete yet! :construction:

# Welcome!
This is the GitHub Repository for the study "Predicting fluid and crystallized Intelligence from resting-state EEG". This project is preregistered on [OSF](https://doi.org/10.17605/OSF.IO/JB64Z). 

You will find all the code in the respective folders. To reproduce our findings (analysis is not finished yet), simply clone the repository and request the data from the authors.

## Analysis Pipeline
Make sure to stay in the directory of the script you run (e.g. Matlab script $\rightarrow$ Matlab directory), otherwise the script will have problems finding the other files.
1. Run the script `Run_preproc.m` - This preprocesses the EEG data.
2. Switch to the Python folder and run `FOOOF_python.py` - This parameterizes the total EEG signal into its periodic and aperiodic signal components and aperiodic parameters. We used the [FOOOF Toolbox](https://fooof-tools.github.io/fooof/) (Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A., Wallis, J. D., Knight, R. T., Shestyuk, A. Y. & Voytek, B. (2020). Parameterizing neural power spectra into periodic and aperiodic components. *Nature Neuroscience, 23(12)*, 1655–1665. https://doi.org/10.1038/s41593-020-00744-x).
3. Switch back to the Matlab directory and run `sort_files.m` - This ensures that the EEG channels are all in the same order.
4. Run `import_data_all_conds.m` - This reads in the EEG and behavioral data and creates files that are used for subsequent MVPA.
5. To start the Support Vector Regressions, start `Run_SVR` - We use the Decision Decoding Toolbox [DDTBox](https://github.com/DDTBOX/DDTBOX) (Bode, S., Feuerriegel, D., Bennett, D., & Alday, P. M. (2019). The Decision Decoding ToolBOX (DDTBOX) – A multivariate pattern analysis toolbox for event-related potentials. *Neuroinformatics, 17(1)*, 27–42. https://doi.org/10.1007/s12021-018-9375-z).
   - We adapted the file `prepare_my_vectors_erp.m` by implementing a `parfor`-Loop. This enables to run the decoding process on several cores which in turn improves overall run time (new file name: `prepare_my_vectors_erp_parfor.m`).
   - This script will take quite a while especially when running 1,000 permutation tests.

Further instruction will follow as soon as they are available.
