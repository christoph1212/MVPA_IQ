function [corr_matrix, mean_corr_matrix, p] = run_correlation(beh, eeg, signal)

arguments
    beh (1,:) cell
    eeg (1,:) cell
    signal (1,:) char {mustBeMember(signal,{'spectral','parameters'})}
end

if strcmpi(signal, 'parameters')
    
    corr_matrix = zeros(59,1);
    mean_corr_matrix = [];
    p = zeros(59,1);
            
    for i = 1:59
            
            [A, p(i)] = corr(beh{1}, ...
                squeeze(eeg{1,1}(1, i, :)), 'type', 'Spearman', 'Rows', 'pairwise');
    
            corr_matrix(i) = A(1);
    
    end

elseif strcmpi(signal, 'spectral')

    frequency_bands = {1:8; 8:16; 16:26; 26:60};

    corr_matrix = zeros(60,59);
    mean_corr_matrix = zeros(4,59);
    p = zeros(4,59);
    
    for i = 1:length(frequency_bands)
        
        for j = 1:59
            
            [A, p(i, j)] = corr(beh{1}, ...
                squeeze(mean(eeg{1,1}(frequency_bands{i}, j, :), 1)), 'type', 'Spearman', 'Rows', 'pairwise');
    
            mean_corr_matrix(i, j) = A(1);

        end
    
    end

    for i = 1:60
        
        for j = 1:59
            
            [A] = corr(beh{1}, ...
                squeeze(eeg{1,1}(i, j, :)), 'type', 'Spearman', 'Rows', 'pairwise');
    
           corr_matrix(i, j) = A(1);

        end
    
    end

end

end