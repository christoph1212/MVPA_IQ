function prepare_svr_postEC(i_behav, behav, folders, Data_sorted_post_EC, Data_sorted_post_EC_male, Data_sorted_post_EC_female, OutputFolderBehav)
%% Prepare behavioral data file for DDTBox for Eyes Closed
temp_post_EC = Data_sorted_post_EC.(behav{i_behav}); 

SVR_labels = {temp_post_EC};
post_EC_file_name = fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesClosed', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesClosed')))
    mkdir(fullfile(OutputFolderBehav, 'Full_Sample', folders(i_behav), 'Post', 'EyesClosed'))
end  

save(post_EC_file_name, 'SVR_labels');

% Repeat for Male
temp_post_EC_male = Data_sorted_post_EC_male.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EC_male};
post_EC_male_file_name = fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesClosed', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesClosed')))
    mkdir(fullfile(OutputFolderBehav, 'Male', folders(i_behav), 'Post', 'EyesClosed'))
end

save(post_EC_male_file_name, 'SVR_labels');

% Repeat for Female
temp_post_EC_female = Data_sorted_post_EC_female.(behav{i_behav});

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EC_female};
post_EC_female_file_name = fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesClosed', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesClosed')))
    mkdir(fullfile(OutputFolderBehav, 'Female', folders(i_behav), 'Post', 'EyesClosed'))
end

save(post_EC_female_file_name, 'SVR_labels');