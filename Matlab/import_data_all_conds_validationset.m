%%
% This file imports (1) EEG data and (2) Behavioral data from the 
% validation setand stores it ina readable format for the DDTBox. This 
% script was created using the original scripts from the study from Jach et
% al. (2020) and was modified for the datasets from Ociepka et al. (2022).
%
% Input:
% (1) Ordered preprocessed EEG data files for the three (?) signals and
%     aperiodic parameters
% (2) Fluid Intelligence data for each participant
%
% Output:
% (1) 'eeg_sorted_cond.mat' containing the EEG data readable for the DDTBox
%     for each condition, signal type and aperiodic parameter
% (2) a list of participants that were included in the file
% (3) 'eeg_sorted_cond_regress_sorted_cond.mat' containing gf data 
%     readable for the DDTBox for each trait
%
% Author: Christoph Fruehlinger (July, 2024)

%% EEG Data
% Set Folders
Datafolder = fullfile('../DataValidation');
InputFolder = fullfile(Datafolder, 'Preprocessed_ordered_data');
OutputFolderEEG = fullfile(Datafolder, 'Ready_for_DDTBOX', 'EEG_sorted_cond');
BehavFolder = fullfile('../ValidationSet'); 


% Load Sociodemographics to sort by Gender
columns = {'ID', 'RAVEN', 'TAO', 'FIG', 'NUMB', 'CATTEL', 'PAPER', 'PLEC'};
Data = readtable(fullfile(BehavFolder, 'IAF_and_beh_data.csv'), 'VariableNamingRule','preserve');
Data = Data(:,columns);
Data = renamevars(Data, 'PLEC', 'Gender');
Data.Gender = string(Data.Gender);
Data.Gender(strcmp(Data.Gender, "K")) = "female";
Data.Gender(strcmp(Data.Gender, "M")) = "male";
[~, score, ~] = pca(table2array(Data(:, 2:7)));
Data.gf_factor_score = score(:,1);


% List of all Participants per condition
subs_pre_EO = {};
subs_post_EO = {};

eeg_sorted_cond_female = [];
eeg_sorted_cond_male = [];

for i_signal = ["Total", "Periodic"]
    for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Post", "EyesOpen")]
        filepath = fullfile(InputFolder, i_signal, i_cond);
        Participants = dir(filepath);
        Participants = Participants(~ismember({Participants.name}, {'.', '..'}));
        Participants = {Participants.name};
        

        % List of Participants with complete data (EEG + Behavioral)
        completePatterns = {};

        for subj = 1:length(Participants)
            % select Participant and get sub-ID
            i_Participant = Participants{subj};
            pattern = split(i_Participant, '_');
            pattern = pattern{2};

            % Append to list complete data
            if ismember(str2double(pattern), Data.ID)
                completePatterns{end+1} = i_Participant;
            end

        end % for subj

        idxFem = 1;
        idxMale = 1;

        % Store EEG data in a readable format for DDTBox
        for i_file = 1:length(completePatterns)
            fprintf('\nRunning Participant %s/%s, %s\n', num2str(i_file), num2str(length(completePatterns)), i_cond)
            temp = csvread(fullfile(filepath, completePatterns{i_file}), 1,1);            

            % Create matrix of frequency bands by channels  by participants
            double(:,:,i_file) = [temp(:,:)'];

            eeg_sorted_cond = {double};
            
            % Create gender-specific subsets
            pattern = split(i_Participant, '_');
            pattern = pattern{2};
                
            idx = find(ismember(Data.ID, str2double(pattern)));
            if strcmp(Data.Gender(idx), "female")       
                double_fem(:,:,idxFem) = [temp(:,:)'];
                eeg_sorted_cond_female = {double_fem};
                idxFem = idxFem + 1;
            else          
                double_male(:,:,idxMale) = [temp(:,:)'];
                eeg_sorted_cond_male = {double_male};
                idxMale = idxMale + 1;
            end
        
        end % for i_file
        clear double double_male double_fem
        
        % Create list of Participants
        sub_IDs = split(completePatterns, '_');
        if strcmp(sub_IDs(:,:,3), "pre")
            subs_pre_EO = sub_IDs(:,:,2);
        else
            subs_post_EO = sub_IDs(:,:,2);
        end

        % Save data in location 'Ready for DDTBox'
        if ~(isfolder(fullfile(OutputFolderEEG, i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, i_signal, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond))
        end
        save(fullfile(OutputFolderEEG, i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');
        save(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond_female');
        save(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond_male');
    end % for i_cond
end % for i_signal

subs_female = Data.ID(strcmp(Data.Gender, "female"));
subs_male = Data.ID(strcmp(Data.Gender, "male"));

subs_pre_EO = subs_pre_EO';
subs_post_EO = subs_post_EO';

subs_pre_EO_female = subs_pre_EO(ismember(str2double(subs_pre_EO), subs_female));
subs_post_EO_female = subs_post_EO(ismember(str2double(subs_post_EO), subs_female));

subs_pre_EO_male = subs_pre_EO(ismember(str2double(subs_pre_EO), subs_male));
subs_post_EO_male = subs_post_EO(ismember(str2double(subs_post_EO), subs_male));

% Save Participant List
writecell(subs_pre_EO, fullfile(Datafolder, 'subs_pre_EO_val.txt'));
writecell(subs_post_EO, fullfile(Datafolder, 'subs_post_EO_val.txt'));
writecell(subs_pre_EO_female, fullfile(Datafolder, 'subs_pre_EO_female_val.txt'));
writecell(subs_post_EO_female, fullfile(Datafolder, 'subs_post_EO_female_val.txt'));
writecell(subs_pre_EO_male, fullfile(Datafolder, 'subs_pre_EO_male_val.txt'));
writecell(subs_post_EO_male, fullfile(Datafolder, 'subs_post_EO_male_val.txt'));

%% Repeat for Aperiodic-Parameters (different Folder Strucutre)

fprintf("\nFinished with Power Spectra.\nStarting Parameters.\n")
clear eeg_sorted_cond_male eeg_sorted_cond_female eeg_sorted_cond temp

InputFolder = fullfile(Datafolder, 'Preprocessed_ordered_data', 'Aperiodic', 'Parameters');
OutputFolderEEG = fullfile(Datafolder, 'Ready_for_DDTBOX', 'EEG_sorted_cond');

eeg_sorted_cond_female = [];
eeg_sorted_cond_male = [];

for i_param = ["Exponent", "Offset"]
    for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Post", "EyesOpen")]
        filepath = fullfile(InputFolder, i_param, i_cond);
        Participants = dir(filepath);
        Participants = Participants(~ismember({Participants.name}, {'.', '..'}));
        Participants = {Participants.name};

        % List of Participants with complete data (EEG + Behavioral)
        completePatterns = {};

        for subj = 1:length(Participants)
            % select Participant and get sub-ID
            i_Participant = Participants{subj};
            pattern = split(i_Participant, '_');
            pattern = pattern{1};

            % If condition for Pre Post
            if ismember(pattern, Data.ID)
                completePatterns{end+1} = i_Participant;
            end

        end % for subj

        idxFem = 1;
        idxMale = 1;

        for i_file = 1:length(completePatterns)
            fprintf('\nRunning Participant %s/%s, %s, %s\n', num2str(i_file), num2str(length(completePatterns)), i_param, i_cond)
            params = csvread(fullfile(filepath, completePatterns{i_file}), 1,1);

            params_all(:,:,i_file) = [params(:)'];               

            eeg_sorted_cond = {params_all};

            % Create gender-specific subsets
            tokens = regexp(completePatterns{i_file}, '^(rest_[^_]+)', 'tokens');           
            participantID = tokens{1}{1};

            idx = find(strcmp(Data.ID, participantID));
            if strcmp(Data.Gender(idx), "female")       
                params_fem(:,:,idxFem) = [params(:,:)'];
                eeg_sorted_cond_female = {params_fem};
                idxFem = idxFem + 1;
            else          
                params_male(:,:,idxMale) = [params(:,:)'];
                eeg_sorted_cond_male = {params_male};
                idxMale = idxMale + 1;
            end

        end % for i_file
        clear params_all params_fem params_male

        % Save data in location 'Ready for DDTBox'
        if ~(isfolder(fullfile(OutputFolderEEG, 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        save(fullfile(OutputFolderEEG, 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');
        save(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond_female');
        save(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond_male');
    end % for i_cond
end % for i_param

%% Behavioral Data

fprintf("\nFinished with EEG Data.\nStarting Behavioral Data.\n")

% Create Folders
OutputFolderBehav = fullfile(Datafolder, 'Ready_for_DDTBOX', 'Behavioral_scores');
if ~(isfolder(fullfile(OutputFolderBehav)))
    mkdir(fullfile(OutputFolderBehav))
end

% Sort Participants in Behav Data to match EEG order
[~, sortedIndices] = ismember(str2double(subs_pre_EO), Data.ID);
Data_sorted_pre_EO = Data(sortedIndices, :);

[~, sortedIndices] = ismember(str2double(subs_pre_EO_male), Data.ID);
Data_sorted_pre_EO_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(str2double(subs_pre_EO_female), Data.ID);
Data_sorted_pre_EO_female = Data(sortedIndices, :);

[~, sortedIndices] = ismember(str2double(subs_post_EO), Data.ID);
Data_sorted_post_EO = Data(sortedIndices, :);

[~, sortedIndices] = ismember(str2double(subs_post_EO_male), Data.ID);
Data_sorted_post_EO_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(str2double(subs_post_EO_female), Data.ID);
Data_sorted_post_EO_female = Data(sortedIndices, :);



%% Prepare behavioral data file for DDTBox for Pre Eyes Open
% Select the relevant trait for Pre Eyes Open
temp_pre_EO = Data_sorted_pre_EO.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO};
pre_EO_file_name = fullfile(OutputFolderBehav, 'Fluid', 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Fluid', 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Fluid', 'Pre', 'EyesOpen'))
end

save(pre_EO_file_name, 'SVR_labels');

% Repeat for Male
temp_pre_EO_male = Data_sorted_pre_EO_male.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO_male};
pre_EO_male_file_name = fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Pre', 'EyesOpen'))
end

save(pre_EO_male_file_name, 'SVR_labels');

% Repeat for Female
temp_pre_EO_female = Data_sorted_pre_EO_female.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_pre_EO_female};
pre_EO_female_file_name = fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Pre', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Pre', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Pre', 'EyesOpen'))
end

save(pre_EO_female_file_name, 'SVR_labels');


%% Repeat for Post Eyes Open
% Select the relevant trait for Post Eyes Open
temp_post_EO = Data_sorted_post_EO.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EO};
post_EO_file_name = fullfile(OutputFolderBehav, 'Fluid', 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Fluid', 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Fluid', 'Post', 'EyesOpen'))
end

save(post_EO_file_name, 'SVR_labels');

% Repeat for Male
temp_post_EO_male = Data_sorted_post_EO_male.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EO_male};
post_EO_male_file_name = fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Male', 'Fluid', 'Post', 'EyesOpen'))
end

save(post_EO_male_file_name, 'SVR_labels');

% Repeat for Female
temp_post_EO_female = Data_sorted_post_EO_female.gf_factor_score;

% Add files to DDTBox-readable format & save
SVR_labels = {temp_post_EO_female};
post_EO_female_file_name = fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Post', 'EyesOpen', 'eeg_sorted_cond_regress_sorted_cond');

if ~(isfolder(fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Post', 'EyesOpen')))
    mkdir(fullfile(OutputFolderBehav, 'Female', 'Fluid', 'Post', 'EyesOpen'))
end

save(post_EO_female_file_name, 'SVR_labels');


disp("Finished")