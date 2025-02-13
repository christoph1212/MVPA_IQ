%%
% This file imports the channel order from a file and orders the remaining 
% files according to the template. This script was created using the
% original script from the study from Jach et al. (2020) and was modified
% for the CoScience datasets. Change the information about Datafolder 
% according to your directories.
%
% Input: 
% (1) Preprocessed EEG data files for the three signal types and aperiodic
%     parameters created by the Preprocessing and FOOOF script
%
% Output:
% (1) Ordered preprocessed EEG dat files for the three signals and
%     aperiodic parametersgit
%
% Author: Christoph Fruehlinger (August, 2023)

% Set Folders
DataFolder = '../DataValidation/';
InputFolder = fullfile(DataFolder, 'Preprocessed');
OutputFolder = fullfile(DataFolder, 'Preprocessed_ordered_data');

% Create a Template to Order Files
addpath(DataFolder);
Template = readtable(fullfile(InputFolder, 'Pre', 'EyesOpen', 'rest_1_pre_Average_EO.csv'));
Template = Template(:,1);


% Load Files and order Channels according to the Template
for i_signal = ["Total", "Periodic"]
    for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Post", "EyesOpen")]
        if i_signal == "Total"
            filepath = fullfile(InputFolder, i_cond);
        else
            filepath = fullfile(DataFolder, 'FOOOF_Data', i_signal, i_cond);
        end
    
        files = dir(filepath);
        files = files(~ismember({files.name}, {'.', '..'}));
        files = {files.name};
        for file = 1:length(files)
            current_file = cell2mat(files(file));
            loaded_file = readtable(fullfile(filepath, current_file));
            % ? loaded_file = loaded_file(1:64, 1:61);
            [~, order_indices] = ismember(loaded_file(:,1), Template);
            [~, sorted_indices] = sort(order_indices);
            sorted_data = loaded_file(sorted_indices, :);
            [~,filename,~] = fileparts(current_file);
    
            if ~(isfolder(fullfile(OutputFolder, i_signal, i_cond)))
            mkdir(fullfile(OutputFolder, i_signal, i_cond))
            end
    
            if i_signal == "Total"
                writetable(sorted_data, fullfile(OutputFolder, i_signal, i_cond, [filename '_new.csv']));
            else
                writetable(sorted_data, fullfile(OutputFolder, i_signal, i_cond, ['rest_' filename '_new.csv']));
            end

        end
    end
end


% Repeat for Aperiodic-Paramters

OutputFolder = fullfile(DataFolder, 'Preprocessed_ordered_data', 'Aperiodic', 'Parameters');

for i_cond = [fullfile("Pre", "EyesOpen"), fullfile("Pre", "EyesClosed"), fullfile("Post", "EyesOpen"), fullfile("Post", "EyesClosed")]

    filepath = fullfile(DataFolder, 'FOOOF_Data', "Aperiodic", "Parameter", i_cond);

    files = dir(filepath);
    files = files(~ismember({files.name}, {'.', '..'}));
    files = {files.name};
    for file = 1:length(files)
        current_file = cell2mat(files(file));
        loaded_file = readtable(fullfile(filepath, current_file));
        [~, order_indices] = ismember(loaded_file(:,1), Template);
        [~, sorted_indices] = sort(order_indices);
        sorted_data = loaded_file(sorted_indices, :);
        [~,filename,~] = fileparts(current_file);

        sorted_data_offset = sorted_data(:,1:2);
        sorted_data_exponent = sorted_data(:,[1,3]);

        if ~(isfolder(fullfile(OutputFolder, 'Offset', i_cond)))
            mkdir(fullfile(OutputFolder, 'Offset', i_cond))
        end

        if ~(isfolder(fullfile(OutputFolder, 'Exponent', i_cond)))
            mkdir(fullfile(OutputFolder, 'Exponent', i_cond))
        end

        writetable(sorted_data_offset, fullfile(OutputFolder, 'Offset', i_cond, ['rest_' filename '_offset_new.csv']));
        writetable(sorted_data_exponent, fullfile(OutputFolder, 'Exponent', i_cond, ['rest_' filename '_exponent_new.csv']));
    end

end

chanlocs_all = readtable('Chanlocs_Biosemi.txt');
[~, sort_idx] = ismember(upper(Template.Var1), [chanlocs_all.labels]);
valid_idx = sort_idx > 0;
sorted_chanlocs_all = chanlocs_all(sort_idx(valid_idx), :);
sorted_chanlocs_all = table2struct(sorted_chanlocs_all);

save('chanlocs_sorted_validationset.mat', 'sorted_chanlocs_all')