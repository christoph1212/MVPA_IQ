% --------------------------------------------------------
% --------- Prepare Decoding Results for graphing --------
% --------------------------------------------------------

% Description: 
%   This short script extracts decoding accuracy from MATLAB RESULTS files
%   and averages over repetitions of 10fold cross validation, 
%   and puts it into a single dataset so that we can graph variability 
%   in CV repetitions.

% Began 13 January 2019
% Last changes 26 May 2020

% Written by Hayley Jach


%%  housekeeping
clear all
clc;

%% Add paths to MVPA results data 

ResultsFolder = fullfile('../', 'Results');
addpath(ResultsFolder);

samples = dir(ResultsFolder);
samples = {samples(~ismember({samples.name}, {'.', '..'})).name};

%% Extract results & average

for sample = 1:length(samples)
    files = dir(fullfile(ResultsFolder, samples{sample}, '*.mat'));
    files = files(~ismember({files.name}, {'.', '..'}));
    for file = files' % for the number of files in the folder
        load(fullfile(ResultsFolder, samples{sample}, file.name));
        % pa = prediction accuracy
        pa(:,:,:) = RESULTS.prediction_accuracy{1}(:,:,:); % take out of the cell for easier averaging
        perm_pa(:,:,:) = RESULTS.perm_prediction_accuracy{1}(:,:,:); % ditto for perm
        ten_10_folds = nanmean(pa, 3); % average over repetitions of 10-fold CV
        perm_ten_10_folds = nanmean(perm_pa, 3); % average over permuted reps
        comb_ten_10_folds = [ten_10_folds, perm_ten_10_folds]; % combine real & perm decoding into a DF
        savename = [ResultsFolder, filesep, samples{sample}, filesep, file.name, '_averageCV.csv'];
        csvwrite(savename, comb_ten_10_folds);

        splits = strsplit(file.name, '_');

        if strcmp(splits{4}, 'Parameters')
            warning off
            format long
            fprintf("%s %s %s %s %s: r = %.2f\n", splits{2}, splits{5}, splits{6}, splits{7}, splits{8}, mean(comb_ten_10_folds(:,1:10)))
            clear pa perm_pa ten_10_folds perm_ten_10_folds comb_ten_10_folds
        elseif strcmp(splits{5}, 'Parameters')
            warning off
            format long
            fprintf("%s %s %s %s %s: r = %.2f\n", splits{2}, splits{6}, splits{7}, splits{8}, splits{9}, mean(comb_ten_10_folds(:,1:10)))
            clear pa perm_pa ten_10_folds perm_ten_10_folds comb_ten_10_folds
        else
            clear pa perm_pa ten_10_folds perm_ten_10_folds comb_ten_10_folds
            continue
        end

        clear pa perm_pa ten_10_folds perm_ten_10_folds comb_ten_10_folds

    end % for file 
end % for sample
