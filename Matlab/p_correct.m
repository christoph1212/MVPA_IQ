function p_corrected_matices = p_correct(varargin)

p_values = cellfun(@(x) x(:), varargin, 'UniformOutput', false);
p_values = vertcat(p_values{:});

p_values_corrected = mafdr(p_values, 'BHFDR', true);

p_corrected_matices = cell(size(varargin));
idx = 1;
for i = 1:length(varargin)
    matrix_size = numel(varargin{i});
    p_corrected_matices{i} = reshape(p_values_corrected(idx:idx + matrix_size - 1), size(varargin{i}));
    idx = idx + matrix_size;
end