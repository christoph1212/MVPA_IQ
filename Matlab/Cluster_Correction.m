%%  housekeeping
clear;
clc;
close all;

%% Add paths to MVPA results data 

ResultsFolder = fullfile('../', 'Results');
addpath(ResultsFolder);

samples = dir(fullfile(ResultsFolder, 'MVPA*'));
samples = {samples(~ismember({samples.name}, {'.', '..'})).name};

i_row = 1;
p_val_struct = struct('File', {}, 'Bins', {});

p_val = 0.01;

%% Extract results & average
for sample = 1:length(samples)
    files = dir(fullfile(ResultsFolder, samples{sample}, '*.mat'));
    files = files(~ismember({files.name}, {'.', '..'}) & ~contains({files.name}, 'Parameter'));
    for file = files' % for the number of files in the folder
        load(fullfile(ResultsFolder, samples{sample}, file.name));
        % pa = prediction accuracy
        pa(:,:,:) = RESULTS.prediction_accuracy{1}(:,:,:); % take out of the cell for easier averaging
        perm_pa(:,:,:) = RESULTS.perm_prediction_accuracy{1}(:,:,:); % ditto for perm
        ten_10_folds = nanmean(pa, 3); % average over repetitions of 10-fold CV
               
        % Get mean across CV repetitions
        mean_ten_10_folds = mean(ten_10_folds, 2);
        mean_perm_ten_10_folds = squeeze(nanmean(perm_pa,2));
        
        % initialize empty array for max cluster sizes
        max_perm_cluster_sizes = zeros(size(mean_perm_ten_10_folds, 2),1);

        % loop through each permutation matrix
        for i_perm = 1:size(mean_perm_ten_10_folds, 2)
            
            mean_perm = mean_perm_ten_10_folds(:,i_perm);
            
            % convert Fisher-z-transformed coefficients back to r and
            % compute p-value (one-sided; only positive values show
            % association between EEG and behavioral data)
            if sample == 1
                if any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 364;
                elseif any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesClosed'))
                    n = 360;
                elseif any(strcmp(split(file.name, '_'), 'Post')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 362;
                else
                    n = 360;
                end
            elseif sample == 2
                if any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 749;
                elseif any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesClosed'))
                    n = 742;
                elseif any(strcmp(split(file.name, '_'), 'Post')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 742;
                else
                    n = 738;
                end
            else 
                if any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 377;
                elseif any(strcmp(split(file.name, '_'), 'Pre')) && any(strcmp(split(file.name, '_'), 'EyesClosed'))
                    n = 374;
                elseif any(strcmp(split(file.name, '_'), 'Post')) && any(strcmp(split(file.name, '_'), 'EyesOpen'))
                    n = 372;
                else
                    n = 370;
                end
            end
            
            r = (exp(2*mean_perm) - 1) ./ (exp(2*mean_perm) + 1);
    
            t = (r .* sqrt(n-2)) ./ (sqrt(1-r.^2));
    
            perm_p_values = 1 - tcdf(t, n-2);
            
            % create a p-value mask
            perm_p_value_mask = perm_p_values < p_val;
            % compute clusters and save max cluster size
            perm_clusters = bwconncomp(perm_p_value_mask);
    
            if numel(perm_clusters.PixelIdxList) > 0
            
                perm_tempclustsizes = cellfun(@length, perm_clusters.PixelIdxList);
            
                max_perm_cluster_sizes(i_perm) = max(perm_tempclustsizes);
            
            end
        end

        % estimate cluster threshold
        cluster_thresh = prctile(max_perm_cluster_sizes, 100-(100*p_val));
        if cluster_thresh < 2
            cluster_thresh = 2;
        end

        % Initiate empty p-val matrix
        p_values = zeros(size(mean_ten_10_folds));
        % calculate p-value and add to matrix
        for bin = 1:length(mean_ten_10_folds)
            p_values(bin) = sum(mean_perm_ten_10_folds(bin, :) > mean_ten_10_folds(bin)) / 1000;
        end
        % create a p-value mask
        p_value_mask = p_values < p_val;
        
        % calculate clusters
        clusters = bwconncomp(p_value_mask);
        n_clusters = clusters.NumObjects;
        cluster_size = cellfun(@length, clusters.PixelIdxList);
        % apply cluster threshold to data
        for i = 1:n_clusters
            if numel(clusters.PixelIdxList{i}) < cluster_thresh
                p_values(clusters.PixelIdxList{i}) = 1;
            end
        end

        % find significant p-values
        if isempty(find(p_values < p_val, 1))
            rel_bin = 0;
        else
            rel_bin = find(p_values < p_val)';
        end
        % save data in struct
        p_val_struct(i_row).File = {erase(file.name, '.mat')};

        p_val_struct(i_row).Bins = rel_bin;
        i_row = i_row + 1;

    end
end

p_val_table = struct2table(p_val_struct);
writetable(p_val_table, fullfile(ResultsFolder,'p_values_Cluster_corrected.csv'))