function runDDT_param(i_trait, i_signal, i_cond, cross_val_steps, n_rep_cross_val,  permut_rep, i_param, Folder) 
% 
% This script was created using the original script from the study from 
% Jach et al. (2020) and was modified for the CoScience datasets. Change 
% the information about Basefolder according to your directories. This will
% perform SVR on the aperiodic parameter data created in the 
% import_data_all_conds_combined() script. To run SVR on the spectral power
% signals use runDDT().
%
% Important: to run the cross-validation and permutation in parallel, you
% have to have use prepare_my_vectors_erp_parfor() in decoding_erp(). Both 
% functions should be in the DDTBOX-master folder. Check line 600 in 
% decoding_erp() for the correct function.
%
% Input:
% (1) i_trait = integer 1 to 5 according to the Big 5 traits
% (2) i_signal = integer according to the EEG aperiodic signal
% (3) i_cond = integer 1 to 4 according to the conditions
% (5) cross_val_steps = number of cross-validation steps
% (6) n_rep_cross_val = number of cross-validation repetitions
% (7) permut_rep = number of permutation repetitions
% (8) Folder = Inputfolder containing files readable for DDTBox
%
% Output:
% (1) RESULT = structure array containing the results from SVR
%
% EXAMPLE_run_decoding_analyses.m
%
% This script is used for configuring and running decoding analyses in DDTBOX.  
% A brief explanation of each configurable setting is described below.
% More information on the analysis options in this script, as well as 
% a tutorial on how to run MVPA in DDTBOX, can be found in the DDTBOX wiki, 
% available at: https://github.com/DDTBOX/DDTBOX/wiki
%
% Please make copies of this script for your own projects.
% 
% This script calls decoding_erp.m
%
%
% Copyright (c) 2013-2017, Daniel Feuerriegel and contributors 
% 
% This file is part of DDTBOX.
%
% DDTBOX is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

fprintf("The Input for runDDT is: i_trait %s,\n i_param %s,\n i_signal %s,\n i_cond %s,\n cross_val_steps %s,\n n_rep_cross_val %s,\n permut_rep %s\n", ...
    num2str(i_trait), num2str(i_param), num2str(i_signal), num2str(i_cond), num2str(cross_val_steps), num2str(n_rep_cross_val), num2str(permut_rep))

%% Get path and add to file (DDF & EEGLAB)
% note: this command won't work unless you run the entire script at once -
% i.e., not in command window and not as a 'section'.
% Important - add path to where EEGlab is housed on your computer.
traits = ["Conscientiousness", "Extraversion", "NegativeEmotionality", "OpenMindedness", "Agreeableness"];
signalsEEG = ["Total", "Aperiodic", "Periodic"];
parameters = ["Offset", "Exponent"];
Eye_conditions = [fullfile("Pre", "EyesOpen"), fullfile("Pre", "EyesClosed"), fullfile("Post", "EyesOpen"), fullfile("Post", "EyesClosed")];

% For Cluster inputs are usually given as character
if ischar(i_trait) | isstring(i_trait); i_trait = str2double(i_trait); end
if ischar(i_param) | isstring(i_param); i_param = str2double(i_param); end
if ischar(i_signal) | isstring(i_signal); i_signal = str2double(i_signal); end
if ischar(i_cond) | isstring(i_cond); i_cond = str2double(i_cond); end
if ischar(cross_val_steps) | isstring(cross_val_steps); cross_val_steps = str2double(cross_val_steps); end
if ischar(n_rep_cross_val) | isstring(n_rep_cross_val); n_rep_cross_val = str2double(n_rep_cross_val); end
if ischar(permut_rep) | isstring(permut_rep); permut_rep = str2double(permut_rep); end



% Adjust Folders Here
Basefolder = pwd;        
Datafolder = fullfile(Basefolder, 'Data', filesep);    
addpath(genpath(fullfile(Basefolder, 'DDTBOX-master')));
         
        %% Select Subject Datasets and Discrimination Groups (dcgs)


        % Set the subject datasets on which to perform MVPA
        sbj_todo = [1];%:size(eeg_sorted_cond{1,1},3)];
        
        % Enter the discrimination groups (dcgs) for decoding analyses. 
        % Each discrimination group should be in a separate cell entry.
        % Decoding analyses will be run for all dcgs listed here.
        % e.g. dcgs_for_analyses{1} = [1];
        % e.g. dcgs_for_analyses{2} = [3];
        % Two discrimination groups can be entered when performing cross-condition decoding.
        % (SVM trained using the first entry/dcg, tested on the second entry/dcg)
        % e.g. dcgs_for_analyses{1} = [1, 2];
        dcgs_for_analyses = {};
        dcgs_for_analyses{1} = [1];
        %dcgs_for_analyses{2} = [2];
        
        % Perform cross-condition decoding? 
        % 0 = No / 1 = Yes
        
        cross = 0;

        %% Filepaths and Locations of Subject Datasets
        
        % Enter the name of the study (for labeling saved decoding results files)
        study_name = sprintf('MVPA_%s_%s_%s_%s_%s', signalsEEG(i_signal), parameters(i_param), traits(i_trait), [split(Eye_conditions(i_cond), filesep)]');                            
        
        % Base directory path (where single subject EEG datasets and channel locations files are stored)
        bdir = fullfile(Datafolder, Folder, 'EEG_sorted_cond', signalsEEG(i_signal), 'Parameters', parameters(i_param), Eye_conditions(i_cond));
        
        if ~(isfolder(fullfile(Datafolder , 'MVPA_Output_Parameters', filesep)))                
            mkdir(fullfile(Datafolder , 'MVPA_Output_Parameters', filesep))
        end
            % output_dir = [fullfile(Datafolder , 'MVPA_Output', filesep)];
        output_dir = fullfile(Datafolder, 'MVPA_Output_Parameters', filesep);
          
        % Filepaths of single subject datasets (relative to the base directory)
        sbj_code = {fullfile(Datafolder, Folder, "EEG_sorted_cond", signalsEEG(i_signal), 'Parameters', parameters(i_param), Eye_conditions(i_cond), "eeg_sorted_cond")}; % participants EEG
        sbj_code_bfi = {fullfile(Datafolder, Folder, "Personality_scores", traits(i_trait), Eye_conditions(i_cond), "eeg_sorted_cond_regress_sorted_cond")}; % participants BFI
        
        % Automatically calculates number of subjects from the number of data files
        nsbj = size(sbj_code, 1);
        
        % MATLAB workspace name for single subject data arrays and structures
        data_struct_name = 'eeg_sorted_cond'; % Data arrays for use with DDTBOX must use this name as their MATLAB workspace variable name

        %% EEG Dataset Information

        nchannels = 59; % Number of channels
        sampling_rate = 1000; % Data sampling rate in Hz
        pointzero = 0; % Corresponds to the time of the event of interest (e.g. stimulus presentation) relative to the start of the epoch (in ms)
        
        % For plotting single subject temporal decoding results 
        % (not required if performing spatial or spatiotemporal decoding)
        % channel_names_file = ''; % Name of the .mat file containing channel labels and channel locations
        % channellocs = ''; % Path of the directory containing channel information file

        %% Condition and Discrimination Group (dcg) Information

        % Label each condition / category
        % Usage: cond_labels{condition number} = 'Name of condition';
        % Example: cond_labels{1} = 'Correct Responses';
        % Example: cond_labels{2} = 'Error Responses';
        % Condition label {X} corresponds to data in column X of the single subject
        % data arrays.
        cond_labels = {};
        cond_labels{1} = traits(i_trait);                                           
                
        % Discrimination groups
        % Enter the condition numbers of the conditions used in classification analyses.
        % Usage: dcg{discrimination group number} = [condition 1, condition 2];
        % Example: dcg{1} = [1, 2]; to use conditions 1 and 2 for dcg 1
        
        % If performing support vector regression, only one condition number is
        % needed per dcg.
        % SVR example: dcg{1} = [1]; to perform SVR on data from condition 1
        dcg = {};
        dcg{1} = [1];
        
        % Support Vector Regression (SVR) condition labels
        % Enter the array entry containing condition labels for each discrimination
        % group number. The SVR_labels array contains multiple cells, each
        % containing a list of SVR condition labels.
        % Usage: svr_cond_labels{dcg} = [cell number in SVR_labels];
        % Example: svr_cond_labels{1} = [2]; to use SVR labels in cell 2 for dcg 1
        svr_cond_labels = {};
        svr_cond_labels{1} = [1];
                      
        % Label each discrimination group
        % Usage: dcg_labels{Discrimination group number} = 'Name of discrimination group'
        % Example: dcg_labels{1} = 'Correct vs. Error Responses';
        dcg_labels = {};
        dcg_labels{1} = sprintf('Level_of_%s_%s', traits(i_trait), parameters(i_param));
        
        % This section automaticallly fills in various parameters related to dcgs and conditions 
        ndcg = size(dcg, 2);
        nclasses = size(dcg{1}, 2);      
        ncond = size(cond_labels, 2);

        %% Multivariate Classification/Regression Parameters

        analysis_mode = 3; % ANALYSIS mode (1 = SVM classification with LIBSVM / 2 = SVM classification with LIBLINEAR / 3 = SVR with LIBSVM)
        stmode = 3; % SPACETIME mode (1 = spatial / 2 = temporal / 3 = spatio-temporal)
        avmode = 1; % AVERAGE mode (1 = no averaging; use single-trial data / 2 = use run-averaged data). Note: Single trials needed for SVR
        window_width_ms = 1; % Width of sliding analysis window in ms
        step_width_ms = 1; % Step size with which sliding analysis window is moved through the trial
        zscore_convert = 0; % Convert data into z-scores before decoding? 0 = no / 1 = yes
        %cross_val_steps = 10; % How many cross-validation steps (if no runs available)?
        %n_rep_cross_val = 10; % How many repetitions of full cross-validation with re-ordered data?
        perm_test = 1; % Run decoding using permuted condition labels? 0 = no / 1 = yes
        %permut_rep = 5; % How many repetitions of full cross-validation for permuted labels analysis?        
        
        % Feature weights extraction
        feat_weights_mode = 1; % Extract feature weights? 0 = no / 1 = yes
        
        % Single subject decoding results plotting
        display_on = 0; % Display single subject decoding performance results? 0 = no / 1 = yes
        perm_disp = 0; % Display the permuted labels decoding results in figure? 0 = no / 1 = yes

        %% Copy All Settings Into the cfg Structure
        % No user input required in this section
        cfg = struct();
        cfg.bdir = bdir;
        cfg.output_dir = output_dir;
        cfg.sbj_code = sbj_code;
        cfg.nsbj = nsbj;
        cfg.data_struct_name = data_struct_name;
        cfg.nchannels = nchannels;
        %cfg.channel_names_file = channel_names_file;
        %cfg.channellocs = channellocs;
        cfg.sampling_rate = sampling_rate;
        cfg.pointzero = pointzero;
        cfg.cond_labels = cond_labels;
        cfg.dcg = dcg;
        cfg.dcg_labels = dcg_labels;
        cfg.svr_cond_labels = svr_cond_labels;
        cfg.ndcg = ndcg;
        cfg.nclasses = nclasses;
        cfg.ncond = ncond;
        cfg.study_name = study_name;
        cfg.cross = cross;
        cfg.analysis_mode = analysis_mode;
        cfg.stmode = stmode;
        cfg.avmode = avmode;
        cfg.window_width_ms = window_width_ms;
        cfg.step_width_ms = step_width_ms;
        cfg.zscore_convert = zscore_convert;
        cfg.perm_test = perm_test;
        cfg.cross_val_steps = cross_val_steps;
        cfg.n_rep_cross_val = n_rep_cross_val;
        cfg.permut_rep = permut_rep;
        cfg.feat_weights_mode = feat_weights_mode;
        cfg.display_on = display_on;
        cfg.perm_disp = perm_disp;

        %% Run the Decoding Analyses For Specified Subjects and dcgs
        
        for dcg_set = 1:length(dcgs_for_analyses)
            
            dcg_todo = [];
            dcg_todo = dcgs_for_analyses{dcg_set};
                
            for sbj = sbj_todo
        
                % Save subject and dcg numbers into the configuration settings
                % structure
                cfg.sbj = sbj;
                cfg.dcg_todo = dcg_todo;
                
                % Set subject-specific filepaths for opening and saving files
                cfg.data_open_name = [strcat(sbj_code{sbj}, ".mat")];
                cfg.data_save_name = [strcat(sbj_code{sbj}, "_data.mat")];
                cfg.regress_label_name = [strcat(sbj_code_bfi{sbj}, '.mat')]; % Filepath for regression labels file
        
                % Run the decoding analyses
                decoding_erp(cfg);

            end % of for sbj
            
        end % of for dcg_set
end % function
  