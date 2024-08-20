function Data_trimmed = exclude_outliers(Data)
% This functions reads in the behavioral data and excludes outliers for
% fluid and crystallized intelligence.
%
% Input: Data (Table)
% Output: Data_trimmed (Table) - outlier corrected

gc_cutoff = mean(Data.gc_score) - 3.3 * std(Data.gc_score);

gc_outliers = Data(Data.gc_score < gc_cutoff,:);

gf_cutoff = mean(Data.gf_score) - 3.3 * std(Data.gf_score);

gf_outliers = Data(Data.gf_score < gf_cutoff,:);

outliers = outerjoin(gc_outliers, gf_outliers, 'MergeKeys', true);

Data_trimmed = Data(~ismember(Data.ID, outliers.ID), :);

end