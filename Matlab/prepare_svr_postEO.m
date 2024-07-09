function prepare_svr_postEO(i_behav, behav, folders, Data_sorted_post_EO, Data_sorted_post_EO_male, Data_sorted_post_EO_female, OutputFolderBehav)
%% Prepare behavioral data file for DDTBox for post Eyes Open
temp_post_EO = Data_sorted_post_EO.(behav{i_behav}); 

SVR_labels = {temp_post_EO};
post_EO_file_name = fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesOpen'))
end  

save(post_EO_file_name, 'SVR_labels');

% Repeat for Male
temp_post_EO_male = Data_sorted_post_EO_male.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EO_male};
post_EO_male_file_name = fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesOpen'))
end

save(post_EO_male_file_name, 'SVR_labels');

% Repeat for Female
temp_post_EO_female = Data_sorted_post_EO_female.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EO_female};
post_EO_female_file_name = fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesOpen'))
end

save(post_EO_female_file_name, 'SVR_labels');