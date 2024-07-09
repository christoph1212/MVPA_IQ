% This script performs Support Vector Regression (1) on the spectral power
% signals and (2) on the aperiodic parameters. 

% Settings for SVR
cross_val_steps = 10;
n_rep_cross_val = 10;
permut_rep = 1000;
Folder = 'Ready_for_DDTBOX';

% Run SVR on spectral power signals
for i_trait = 1:5
    for i_signal = 1:3
        for i_cond = 1:4
            runDDT(i_trait, i_signal, i_cond, cross_val_steps, n_rep_cross_val,  permut_rep, Folder)
        end
    end
end

clear i_signal
clear i_cond
clear i_trait

% 2 = Aperiodic
i_signal = 2;
% Run SVR on aperiodic parameters
for i_trait = 1:5
    for i_cond = 1:4
        runDDT_param(i_trait, i_signal, i_cond, cross_val_steps, n_rep_cross_val,  permut_rep, Folder)
    end
end