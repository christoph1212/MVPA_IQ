% This Run-Script initialises the Preprocessing Pipeline.
%
% Written by: Christoph Fruehlinger (March, 2024)

clear
clc

%% Preprocessing

% We used the two parameters for parallelisation of the preprocessing. You
% don't need to change the default settings.

NrSubsets = '1';
ThisSubset = '1';

preproc(NrSubsets, ThisSubset)
