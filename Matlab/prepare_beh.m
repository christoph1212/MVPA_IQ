function behav_data = prepare_beh(Inputfolder)
% This function creates a table with the combined behavioural data.
% Sleepiness scores are calculated by averaging the z-standardized scores
% from tiredness and exhaustion for the self-report measurements that are
% closest to the interested EEG rest-measurements. 
%
% Input: Folderpath to behavioural data
% Output: Table with combined behavioural data

% Get Participant IDs
Participants = dir(fullfile(Inputfolder, 'task-IST'));
Participants = {Participants(~ismember({Participants.name}, {'.', '..'})).name};

% Initiate Table
behav_data = table();
behav_data.ID = Participants';

% Loop through 1st Self-Report (Closest to first EEG rest measurement)
for i_sub = 1:length(Participants)

    file = dir(fullfile(Inputfolder, 'task-SR', Participants{i_sub}, 'beh', '*1_beh.csv'));

    if isempty(file)
        behav_data.Tired_Pre(i_sub) = NaN;
        behav_data.Exhausted_Pre(i_sub) = NaN;
        continue
    end

    sr = readtable(fullfile(file.folder, file.name), 'VariableNamingRule','preserve');

    if isempty(sr)
        behav_data.Tired_Pre(i_sub) = NaN;
        behav_data.Exhausted_Pre(i_sub) = NaN;
        continue
    end

    behav_data.Tired_Pre(i_sub) = sr.("tired.response");
    behav_data.Exhausted_Pre(i_sub) = sr.("exhausted.response");

end


% Loop through 4th Self-Report (Closest to second EEG rest measurement)
for i_sub = 1:length(Participants)

    file = dir(fullfile(Inputfolder, 'task-SR', Participants{i_sub}, 'beh', '*4_beh.csv'));
    
    if isempty(file)
        behav_data.Tired_Post(i_sub) = NaN;
        behav_data.Exhausted_Post(i_sub) = NaN;
        continue
    end

    sr = readtable(fullfile(file.folder, file.name), 'VariableNamingRule','preserve');

    if isempty(sr)
        behav_data.Tired_Post(i_sub) = NaN;
        behav_data.Exhausted_Post(i_sub) = NaN;
        continue
    end

    behav_data.Tired_Post(i_sub) = sr.("tired.response");
    behav_data.Exhausted_Post(i_sub) = sr.("exhausted.response");

end


% Loop through Fluid Intelligence Data
fluid_cor = readtable(fullfile(Inputfolder, "IST_fluid_A.xlsx"));           % Check if Version is correct!
fluid_cor = table2array(fluid_cor(:,2));

for i_sub = 1:length(Participants)

    file = dir(fullfile(Inputfolder, 'task-IST', Participants{i_sub}, 'beh', '*Fluid_beh.csv'));
    
    if isempty(file)
        behav_data.gf_score(i_sub) = NaN;
        continue
    end
    
    opts = detectImportOptions(fullfile(file.folder, file.name));
    opts.SelectedVariableNames = "ratings";
    gf = readtable(fullfile(file.folder, file.name), opts);
    gf = table2array(gf(2:end-1, :));

    if isempty(gf)
        behav_data.gf_score(i_sub) = NaN;
        continue
    end

    behav_data.gf_score(i_sub) = sum(strcmp(fluid_cor, gf));

end

% Loop through Crystallized Intelligence Data
cryst_cor = readtable(fullfile(Inputfolder, "IST_crystallized_A.xlsx"));    % Check if Version is correct!
cryst_cor = table2array(cryst_cor(:,2));

for i_sub = 1:length(Participants)

    file = dir(fullfile(Inputfolder, 'task-IST', Participants{i_sub}, 'beh', '*Crystallized_beh.csv'));
    
    if isempty(file)
        behav_data.gc_score(i_sub) = NaN;
        continue
    end
    
    opts = detectImportOptions(fullfile(file.folder, file.name));
    if ~any(strcmp('ratings',opts.VariableNames))
        behav_data.gc_score(i_sub) = NaN;
        continue
    end
    opts.SelectedVariableNames = "ratings";
    gc = readtable(fullfile(file.folder, file.name), opts);
    gc = table2array(gc(2:85, :));

    if isempty(gc)
        behav_data.gc_score(i_sub) = NaN;
        continue
    end

    behav_data.gc_score(i_sub) = sum(strcmp(cryst_cor, gc));

end


% Create Sleepiness Scores
behav_data = behav_data(~any(ismissing(behav_data),2),:);
behav_data.Pre_Task_Sleepiness = mean([zscore(behav_data.Tired_Pre), zscore(behav_data.Exhausted_Pre)],2);
behav_data.Post_Task_Sleepiness = mean([zscore(behav_data.Tired_Post), zscore(behav_data.Exhausted_Post)],2);
