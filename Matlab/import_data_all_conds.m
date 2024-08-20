%%
% This file imports (1) EEG data and (2) Behavioral data and stores it in
% a readable format for the DDTBox. This script was created using the
% original scripts from the study from Jach et al. (2020) and was modified
% for the CoScience datasets.
%
% Input:
% (1) Ordered preprocessed EEG data files for the three signals and
%     aperiodic parameters
% (2) Behavioral data for each participant
%
% Output:
% (1) 'eeg_sorted_cond.mat' containing the EEG data readable for the DDTBox
%     for each condition, signal type and aperiodic parameter
% (2) a list of participants that were included in the file
% (3) 'eeg_sorted_cond_regress_sorted_cond.mat' containing Behavioral data 
%     readable for the DDTBox for each trait
%
% Author: Christoph Fruehlinger (August, 2023)

%%
clear
clc

%% EEG Data
% Set Folders
Datafolder = fullfile('../Data');
InputFolder = fullfile(Datafolder, 'Preprocessed_ordered_data');
OutputFolderEEG = fullfile(Datafolder, 'Ready_for_DDTBOX', 'EEG_sorted_cond');
BehavFolder = fullfile(Datafolder, 'BehavioralData'); 

% Create Table with Self-Report and IQ Data
disp('Preparing Behavioral Data...')
Behav_data = prepare_beh(BehavFolder);
disp('Done.')

% Load Sociodemographics to sort by Gender
GenderData = readtable(fullfile(Datafolder, "SocioDemographics.txt"));
GenderData = GenderData(:,1:3);
GenderData.Gender = string(GenderData.Gender);
GenderData.Gender(strcmp(GenderData.Gender, "1")) = "female";
GenderData.Gender(strcmp(GenderData.Gender, "2")) = "male";

Data = innerjoin(GenderData, Behav_data, "Keys","ID");

% Exclude outliers
Data = exclude_outliers(Data);

% Calculate z-scores for crystallized intelligence within gender
Data.gc_score_z = nan(height(Data), 1);
Data.gc_score_z(strcmp(Data.Gender, 'male')) = zscore(Data(strcmp(Data.Gender, "male"), :).gc_score);
Data.gc_score_z(strcmp(Data.Gender, 'female')) = zscore(Data(strcmp(Data.Gender, "female"), :).gc_score);

% Save Data Table for Descriptive Analysis
writetable(Data, fullfile(BehavFolder, "Behavioral_Data.csv"));

% List of all Participants per condition
subs_pre_EO = {};
subs_pre_EC = {};
subs_post_EO = {};
subs_post_EC = {};

eeg_sorted_cond_female = [];
eeg_sorted_cond_male = [];

for i_signal = ["Total", "Aperiodic", "Periodic"]
    for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Pre", "EyesClosed"), fullfile("Post", "EyesOpen"), fullfile("Post", "EyesClosed")]
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
            pattern = pattern{1};

            % Append to list complete data
            if ismember(pattern, Data.ID)
                completePatterns{end+1} = i_Participant;
            end

        end % for subj

        idxFem = 1;
        idxMale = 1;

        % Store EEG data in a readable format for DDTBox
        for i_file = 1:length(completePatterns)
            fprintf('\nRunning Participant %s/%s, %s, %s\n', num2str(i_file), num2str(length(completePatterns)), i_signal, i_cond)
            temp = csvread(fullfile(filepath, completePatterns{i_file}), 1,1);            

            % Create matrix of frequency bands by channels  by participants
            double(:,:,i_file) = [temp(:,:)'];

            eeg_sorted_cond = {double};
            
            % Create gender-specific subsets
            tokens = regexp(completePatterns{i_file}, '^(sub-[^_]+)', 'tokens');           
            participantID = tokens{1}{1};
                
            idx = find(strcmp(Data.ID, participantID));
            if strcmp(Data.Gender(idx), "female")       
                double_fem(:,:,idxFem) = [temp(:,:)'];
                eeg_sorted_cond_female = {double_fem};
                idxFem = idxFem + 1;
            elseif strcmp(Data.Gender(idx), "male")
                double_male(:,:,idxMale) = [temp(:,:)'];
                eeg_sorted_cond_male = {double_male};
                idxMale = idxMale + 1;
            end
        
        end % for i_file
        clear double double_male double_fem
        
        % Create list of Participants
        sub_IDs = split(completePatterns, '_');
        if strcmp(sub_IDs(:,:,2), "pre")
            if strcmp(sub_IDs(:,:,4), 'EO')
                subs_pre_EO = sub_IDs(:,:,1);
            else
                subs_pre_EC = sub_IDs(:,:,1);
            end
        else
            if strcmp(sub_IDs(:,:,4), 'EO')
                subs_post_EO = sub_IDs(:,:,1);
            else
                subs_post_EC = sub_IDs(:,:,1);
            end
        end

        % Save data in location 'Ready for DDTBox'
        if ~(isfolder(fullfile(OutputFolderEEG, 'Full_Sample', i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Full_Sample', i_signal, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond))
        end
        save(fullfile(OutputFolderEEG, 'Full_Sample', i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');

        eeg_sorted_cond = eeg_sorted_cond_female;
        save(fullfile(OutputFolderEEG, 'Female', i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');

        eeg_sorted_cond = eeg_sorted_cond_male;
        save(fullfile(OutputFolderEEG, 'Male', i_signal, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');
    end % for i_cond
end % for i_signal

subs_female = Data.ID(strcmp(Data.Gender, "female"));
subs_male = Data.ID(strcmp(Data.Gender, "male"));

subs_pre_EO = subs_pre_EO';
subs_pre_EC = subs_pre_EC';
subs_post_EO = subs_post_EO';
subs_post_EC = subs_post_EC';

subs_pre_EO_female = subs_pre_EO(ismember(subs_pre_EO, subs_female));
subs_pre_EC_female = subs_pre_EC(ismember(subs_pre_EC, subs_female));
subs_post_EO_female = subs_post_EO(ismember(subs_post_EO, subs_female));
subs_post_EC_female = subs_post_EC(ismember(subs_post_EC, subs_female));

subs_pre_EO_male = subs_pre_EO(ismember(subs_pre_EO, subs_male));
subs_pre_EC_male = subs_pre_EC(ismember(subs_pre_EC, subs_male));
subs_post_EO_male = subs_post_EO(ismember(subs_post_EO, subs_male));
subs_post_EC_male = subs_post_EC(ismember(subs_post_EC, subs_male));

% Save Participant List
writecell(subs_pre_EO, fullfile(Datafolder, 'subs_pre_EO.txt'));
writecell(subs_pre_EC, fullfile(Datafolder, 'subs_pre_EC.txt'));
writecell(subs_post_EO, fullfile(Datafolder, 'subs_post_EO.txt'));
writecell(subs_post_EC, fullfile(Datafolder, 'subs_post_EC.txt'));
writecell(subs_pre_EO_female, fullfile(Datafolder, 'subs_pre_EO_female.txt'));
writecell(subs_pre_EC_female, fullfile(Datafolder, 'subs_pre_EC_female.txt'));
writecell(subs_post_EO_female, fullfile(Datafolder, 'subs_post_EO_female.txt'));
writecell(subs_post_EC_female, fullfile(Datafolder, 'subs_post_EC_female.txt'));
writecell(subs_pre_EO_male, fullfile(Datafolder, 'subs_pre_EO_male.txt'));
writecell(subs_pre_EC_male, fullfile(Datafolder, 'subs_pre_EC_male.txt'));
writecell(subs_post_EO_male, fullfile(Datafolder, 'subs_post_EO_male.txt'));
writecell(subs_post_EC_male, fullfile(Datafolder, 'subs_post_EC_male.txt'));

%% Repeat for Aperiodic-Parameters (different Folder Strucutre)

fprintf("\nFinished with Power Spectra.\nStarting Parameters.\n")
clear eeg_sorted_cond_male eeg_sorted_cond_female eeg_sorted_cond temp

InputFolder = fullfile(Datafolder, 'Preprocessed_ordered_data', 'Aperiodic', 'Parameters');
OutputFolderEEG = fullfile(Datafolder, 'Ready_for_DDTBOX', 'EEG_sorted_cond');

eeg_sorted_cond_female = [];
eeg_sorted_cond_male = [];

for i_param = ["Exponent", "Offset"]
    for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Pre", "EyesClosed"), fullfile("Post", "EyesOpen"), fullfile("Post", "EyesClosed")]
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
            tokens = regexp(completePatterns{i_file}, '^(sub-[^_]+)', 'tokens');           
            participantID = tokens{1}{1};

            idx = find(strcmp(Data.ID, participantID));
            if strcmp(Data.Gender(idx), "female")       
                params_fem(:,:,idxFem) = [params(:,:)'];
                eeg_sorted_cond_female = {params_fem};
                idxFem = idxFem + 1;
            elseif strcmp(Data.Gender(idx), "male")
                params_male(:,:,idxMale) = [params(:,:)'];
                eeg_sorted_cond_male = {params_male};
                idxMale = idxMale + 1;
            end

        end % for i_file
        clear params_all params_fem params_male

        % Save data in location 'Ready for DDTBox'
        if ~(isfolder(fullfile(OutputFolderEEG, 'Full_Sample', 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Full_Sample', 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        if ~(isfolder(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond)))
            mkdir(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond))
        end
        save(fullfile(OutputFolderEEG, 'Full_Sample', 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');

        eeg_sorted_cond = eeg_sorted_cond_female;
        save(fullfile(OutputFolderEEG, 'Female', 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');

        eeg_sorted_cond = eeg_sorted_cond_male;
        save(fullfile(OutputFolderEEG, 'Male', 'Aperiodic', 'Parameters', i_param, i_cond, 'eeg_sorted_cond'), 'eeg_sorted_cond');
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
[~, sortedIndices] = ismember(subs_pre_EO, Data.ID);
Data_sorted_pre_EO = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_pre_EO_male, Data.ID);
Data_sorted_pre_EO_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_pre_EO_female, Data.ID);
Data_sorted_pre_EO_female = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_pre_EC, Data.ID);
Data_sorted_pre_EC = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_pre_EC_male, Data.ID);
Data_sorted_pre_EC_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_pre_EC_female, Data.ID);
Data_sorted_pre_EC_female = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EO, Data.ID);
Data_sorted_post_EO = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EO_male, Data.ID);
Data_sorted_post_EO_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EO_female, Data.ID);
Data_sorted_post_EO_female = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EC, Data.ID);
Data_sorted_post_EC = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EC_male, Data.ID);
Data_sorted_post_EC_male = Data(sortedIndices, :);

[~, sortedIndices] = ismember(subs_post_EC_female, Data.ID);
Data_sorted_post_EC_female = Data(sortedIndices, :);

% all Traits
behav = {'gf_score', 'gc_score_z', 'Pre_Task_Sleepiness', 'Post_Task_Sleepiness'};

% List folder names to put traits in
folders = ["Fluid", "Crystallized", "Pre_Sleepiness", "Post_Sleepiness"];


for i_behav = 1:length(behav)
    
    if i_behav < 3 % Fluid + Crystallized
        prepare_svr_preEO(i_behav, behav, folders, Data_sorted_pre_EO, Data_sorted_pre_EO_male, Data_sorted_pre_EO_female, OutputFolderBehav)
    
        prepare_svr_preEC(i_behav, behav, folders, Data_sorted_pre_EC, Data_sorted_pre_EC_male, Data_sorted_pre_EC_female, OutputFolderBehav)
    
        prepare_svr_postEO(i_behav, behav, folders, Data_sorted_post_EO, Data_sorted_post_EO_male, Data_sorted_post_EO_female, OutputFolderBehav)
    
        prepare_svr_postEC(i_behav, behav, folders, Data_sorted_post_EC, Data_sorted_post_EC_male, Data_sorted_post_EC_female, OutputFolderBehav)

    elseif i_behav == 3 % Pre_Task_Sleepiness

        prepare_svr_preEO(i_behav, behav, folders, Data_sorted_pre_EO, Data_sorted_pre_EO_male, Data_sorted_pre_EO_female, OutputFolderBehav)
    
        prepare_svr_preEC(i_behav, behav, folders, Data_sorted_pre_EC, Data_sorted_pre_EC_male, Data_sorted_pre_EC_female, OutputFolderBehav)
    
    else % Post_Task_Sleepiness

        prepare_svr_postEO(i_behav, behav, folders, Data_sorted_post_EO, Data_sorted_post_EO_male, Data_sorted_post_EO_female, OutputFolderBehav)
    
        prepare_svr_postEC(i_behav, behav, folders, Data_sorted_post_EC, Data_sorted_post_EC_male, Data_sorted_post_EC_female, OutputFolderBehav)
    
    end

end % for i_trait

disp("Finished")