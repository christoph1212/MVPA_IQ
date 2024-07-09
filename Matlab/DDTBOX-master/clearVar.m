function clearVar(variable)
% This function is used to clear variables in the workspace. This function
% can be used in parfor loops and overcome the problem with the original
% clear() function. Input must be string.

    evalin('caller', ['clear ' variable]);
end