% This Run-Script initialises the the MVPA Pipeline.
%
% Written by: Christoph Fruehlinger (March, 2024)

clear
clc

addpath(genpath('DDTBOX-master'))

% Cross Validation and Permutation parameters
cross_val_steps = 10;
n_rep_cross_val = 10;
permut_rep      = 1000;

Basefolder = '../';
Folder = "Ready_for_DDTBOX";

% Run MVPA for spectral signal
for i_behav = 1:4                       % gf, gc, Pre/Post Sleepiness
    
    if i_behav == 1 || i_behav == 2
        samples = 1:3;                  % Full, male, female sample
    else
        samples = 1;                    % only full sample
    end

    if i_behav == 3
        cond = 1:2;
    elseif i_behav == 4
        cond = 3:4;
    else
        cond = 1:4;
    end

    for i_signal = 1:4                  % EEG signal types

        for i_cond = cond                % Pre/Post EO/EC

            for i_sample = 1:length(samples)

                runDDT(i_behav, i_signal, i_cond, i_sample, cross_val_steps, n_rep_cross_val, permut_rep, Folder, Basefolder)

            end

        end

    end

end
