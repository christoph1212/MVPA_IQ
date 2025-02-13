# Welcome!
This is the GitHub Repository for the study "Sleepiness but neither fluid nor crystallized intelligence can be predicted from resting-state EEG – Evidence from the large scale CoScience EEG-Personality Project". This project is preregistered on [OSF](https://doi.org/10.17605/OSF.IO/JB64Z). 

You will find all the code in the respective folders. To reproduce our findings (analysis is not finished yet), simply clone the repository and request the data from the authors.

## Analysis Pipeline
Make sure to stay in the directory of the script you run (e.g. Matlab script $\rightarrow$ Matlab directory), otherwise the script will have problems finding the other files.
1. Run the script `Run_preproc.m` - This preprocesses the EEG data.
2. Switch to the Python folder and run `FOOOF_python.py` - This parameterizes the total EEG signal into its periodic and aperiodic signal components and aperiodic parameters. We used the [FOOOF Toolbox](https://fooof-tools.github.io/fooof/) (Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A., Wallis, J. D., Knight, R. T., Shestyuk, A. Y. & Voytek, B. (2020). Parameterizing neural power spectra into periodic and aperiodic components. *Nature Neuroscience, 23(12)*, 1655–1665. https://doi.org/10.1038/s41593-020-00744-x).
3. Switch back to the Matlab directory and run `sort_files.m` - This ensures that the EEG channels are all in the same order.
4. Run `import_data_all_conds.m` - This reads in the EEG and behavioral data and creates files that are used for subsequent MVPA.

   4.1 You can now run `Correlations_Plots.m` and `Partial_Correlation_Plots.m` to create correlation plots such as Figures 2 and 3 in the manuscript.
6. To start the Support Vector Regressions, start `Run_SVR` - We use the Decision Decoding Toolbox [DDTBox](https://github.com/DDTBOX/DDTBOX) (Bode, S., Feuerriegel, D., Bennett, D., & Alday, P. M. (2019). The Decision Decoding ToolBOX (DDTBOX) – A multivariate pattern analysis toolbox for event-related potentials. *Neuroinformatics, 17(1)*, 27–42. https://doi.org/10.1007/s12021-018-9375-z).
   - We adapted the file `prepare_my_vectors_erp.m` by implementing a `parfor`-Loop. This enables to run the decoding process on several cores which in turn improves overall run time (new file name: `prepare_my_vectors_erp_parfor.m`).
   - This script will take quite a while especially when running 1,000 permutation tests. The data from the permutations will be used in `Cluster_Correction.m` to perform exploratory cluster corrections on the MVPA results.
7. The results from the MVPA and exploratory cluster correction will be plotted in `R`. So switch to the R folder and run `Graphs_MVPA_output.R`. This script will create among other Figure 1 and S2-S5.
8. In the R folder, there are also scripts for hypothesis testing and descriptive analysis.
