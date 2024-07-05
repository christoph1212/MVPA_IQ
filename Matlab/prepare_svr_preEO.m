function prepare_svr_preEO(i_behav, behav, folders, Data_sorted_pre_EO, Data_sorted_pre_EO_male, Data_sorted_pre_EO_female, OutputFolderBehav)
%% Prepare behavioral data file for DDTBox for Pre Eyes Open
% Select the relevant trait for Pre Eyes Open
temp_pre_EO = Data_sorted_pre_EO.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO};
pre_EO_file_name = fullfile(OutputFolderBehav, folders(i_behav), 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, folders(i_behav), 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, folders(i_behav), 'Pre', 'EyesOpen'))
end

save(pre_EO_file_name, 'SVR_labels');

% Repeat for Male
temp_pre_EO_male = Data_sorted_pre_EO_male.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO_male};
pre_EO_male_file_name = fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Pre', 'EyesOpen'))
end

save(pre_EO_male_file_name, 'SVR_labels');

% Repeat for Female
temp_pre_EO_female = Data_sorted_pre_EO_female.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO_female};
pre_EO_female_file_name = fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Pre', 'EyesOpen'))
end

save(pre_EO_female_file_name, 'SVR_labels');