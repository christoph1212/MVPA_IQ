%% Calculate Spearman Correlations between Brain Power and Behavioral Data

% Created by Christoph Fruehlinger, 2024 using a script from Daniel
% Feuerriegel and Hayley Jach.

%% Housekeeping

clear
close 
clc

% Important - add your path to EEGLAB
% (or figures will not work)
eeglab
clc

%%  Load behavioral data

% ------ Sleepiness -------
% Sleepiness pre EO
load('../Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Pre_Sleepiness/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat')
sleepiness_pre_EO = SVR_labels;

% Sleepiness pre EC
load('../Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Pre_Sleepiness/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat')
sleepiness_pre_EC = SVR_labels;

% Sleepiness post EO
load('../Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Post_Sleepiness/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat')
sleepiness_post_EO = SVR_labels;

% Sleepiness post EC
load('../Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Post_Sleepiness/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat')
sleepiness_post_EC = SVR_labels;

%%  Load channel locations file as in DDTBOX folder
load('DDTBOX-master/chanlocs_sorted.mat')

%% Start Loop
signal = ["Total", "Periodic", "Offset", "Exponent"];
sample = ["Full_Sample", "Male", "Female"];

for i_sample = 1:length(sample)
    %% remaining behavioral data
    % ------- Fluid -------
    % Fluid pre EO
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Fluid/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    fluid_pre_EO = SVR_labels;
    
    % Fluid pre EC
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Fluid/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    fluid_pre_EC = SVR_labels;
    
    % Fluid post EO
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Fluid/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    fluid_post_EO = SVR_labels;
    
    % Fluid post EC
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Fluid/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    fluid_post_EC = SVR_labels;
    
    % ------- Crystallized -------
    % Crystallized pre EO
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Crystallized/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    crystallized_pre_EO = SVR_labels;
    
    % Crystallized pre EC
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Crystallized/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    crystallized_pre_EC = SVR_labels;
    
    % Crystallized post EO
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Crystallized/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    crystallized_post_EO = SVR_labels;
    
    % Crystallized post EC
    load(sprintf('../Data/Ready_for_DDTBOX/Behavioral_scores/%s/Crystallized/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat', sample(i_sample)))
    crystallized_post_EC = SVR_labels;
    
    for i_signal = 1:length(signal)

        fprintf('Plotting: %s, %s Signal\n', sample(i_sample), signal(i_signal))
        
        %%  EEG data
        if strcmp(signal(i_signal), "Total") || strcmp(signal(i_signal), "Periodic")

            % Pre eyes open 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/%s/Pre/EyesOpen/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            pre_eyes_open = eeg_sorted_cond;
            
            % Pre eyes closed
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/%s/Pre/EyesClosed/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            pre_eyes_closed = eeg_sorted_cond;
            
            % Post eyes open 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/%s/Post/EyesOpen/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            post_eyes_open = eeg_sorted_cond;
                        
            % Post eyes closed 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/%s/Post/EyesClosed/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            post_eyes_closed = eeg_sorted_cond;
            
            signal_option = "spectral";

        else

            % Pre eyes open 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/Aperiodic/Parameters/%s/Pre/EyesOpen/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            pre_eyes_open = eeg_sorted_cond;
                        
            % Pre eyes closed
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/Aperiodic/Parameters/%s/Pre/EyesClosed/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            pre_eyes_closed = eeg_sorted_cond;
                        
            % Post eyes open 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/Aperiodic/Parameters/%s/Post/EyesOpen/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            post_eyes_open = eeg_sorted_cond;
            
            % Post eyes closed 
            load(sprintf('../Data/Ready_for_DDTBOX/EEG_sorted_cond/%s/Aperiodic/Parameters/%s/Post/EyesClosed/eeg_sorted_cond.mat', sample(i_sample), signal(i_signal)))
            post_eyes_closed = eeg_sorted_cond;                    
            
            signal_option = "parameters";

        end

        %% Run Correlations
        if strcmp(signal(i_signal), "Total") || strcmp(signal(i_signal), "Periodic")            
            %% Fluid correlations
            % --------- Pre eyes open ------------
            [corr_fluid_pre_EO, mean_corr_fluid_pre_EO, p_fluid_pre_EO] = run_correlation(fluid_pre_EO, pre_eyes_open, signal_option);
                    
            % ------ Pre eyes closed --------        
            [corr_fluid_pre_EC, mean_corr_fluid_pre_EC, p_fluid_pre_EC] = run_correlation(fluid_pre_EC, pre_eyes_closed, signal_option);
            
            % ----- Post eyes open --------        
            [corr_fluid_post_EO, mean_corr_fluid_post_EO, p_fluid_post_EO] = run_correlation(fluid_post_EO, post_eyes_open, signal_option);
           
            % ----- Post eyes closed -----        
            [corr_fluid_post_EC, mean_corr_fluid_post_EC, p_fluid_post_EC] = run_correlation(fluid_post_EC, post_eyes_closed, signal_option);
            
            % Correct for multiple comparisons
            p_corrected_fluid = p_correct(p_fluid_pre_EO, p_fluid_pre_EC, p_fluid_post_EO, p_fluid_post_EC);
            
            % list conds
            mean_conds_fluid = {};
            mean_conds_fluid{1} = mean_corr_fluid_pre_EO;
            mean_conds_fluid{2} = mean_corr_fluid_pre_EC;
            mean_conds_fluid{3} = mean_corr_fluid_post_EO;
            mean_conds_fluid{4} = mean_corr_fluid_post_EC;
            mean_conds_fluid{5} = p_corrected_fluid{1};
            mean_conds_fluid{6} = p_corrected_fluid{2};
            mean_conds_fluid{7} = p_corrected_fluid{3};
            mean_conds_fluid{8} = p_corrected_fluid{4};

            %% Crystallized correlations
            % --------- Pre eyes open ------------
            [corr_crystallized_pre_EO, mean_corr_crystallized_pre_EO, p_crystallized_pre_EO]  = run_correlation(crystallized_pre_EO, pre_eyes_open, signal_option);
            
            % ------ Pre eyes closed --------        
            [corr_crystallized_pre_EC, mean_corr_crystallized_pre_EC, p_crystallized_pre_EC] = run_correlation(crystallized_pre_EC, pre_eyes_closed, signal_option);
            
            % ----- Post eyes open --------                
            [corr_crystallized_post_EO, mean_corr_crystallized_post_EO, p_crystallized_post_EO] = run_correlation(crystallized_post_EO, post_eyes_open, signal_option);
            
            % ----- Post eyes closed -----        
            [corr_crystallized_post_EC, mean_corr_crystallized_post_EC, p_crystallized_post_EC] = run_correlation(crystallized_post_EC, post_eyes_closed, signal_option);
            
            % Correct for multiple comparisons
            p_corrected_crystallized = p_correct(p_crystallized_pre_EO, p_crystallized_pre_EC, p_crystallized_post_EO, p_crystallized_post_EC);

            % list conds
            mean_conds_crystallized = {};
            mean_conds_crystallized{1} = mean_corr_crystallized_pre_EO;
            mean_conds_crystallized{2} = mean_corr_crystallized_pre_EC;
            mean_conds_crystallized{3} = mean_corr_crystallized_post_EO;
            mean_conds_crystallized{4} = mean_corr_crystallized_post_EC;
            mean_conds_crystallized{5} = p_corrected_crystallized{1};
            mean_conds_crystallized{6} = p_corrected_crystallized{2};
            mean_conds_crystallized{7} = p_corrected_crystallized{3};
            mean_conds_crystallized{8} = p_corrected_crystallized{4};

            if i_sample == 1 
                %% Sleepiness
                % ------ Pre eyes open -------- 
                [corr_sleepiness_pre_EO, mean_corr_sleepiness_pre_EO, p_sleepiness_pre_EO] = run_correlation(sleepiness_pre_EO, pre_eyes_open, signal_option);
                
                % ------ Pre eyes closed --------            
                [corr_sleepiness_pre_EC, mean_corr_sleepiness_pre_EC, p_sleepiness_pre_EC] = run_correlation(sleepiness_pre_EC, pre_eyes_closed, signal_option);
                
                % ----- Post eyes open --------                    
                [corr_sleepiness_post_EO, mean_corr_sleepiness_post_EO, p_sleepiness_post_EO] = run_correlation(sleepiness_post_EO, post_eyes_open, signal_option);
                
                % ----- Post eyes closed -----            
                [corr_sleepiness_post_EC, mean_corr_sleepiness_post_EC, p_sleepiness_post_EC] = run_correlation(sleepiness_post_EC, post_eyes_closed, signal_option);
                
                % Correct for multiple comparisons
                p_corrected_sleepiness = p_correct( ...
                    p_sleepiness_pre_EO, p_sleepiness_pre_EC, p_sleepiness_post_EO, p_sleepiness_post_EC);

                mean_conds_sleepiness = {};
                mean_conds_sleepiness{1} = mean_corr_sleepiness_pre_EO;
                mean_conds_sleepiness{2} = mean_corr_sleepiness_pre_EC;
                mean_conds_sleepiness{3} = mean_corr_sleepiness_post_EO;
                mean_conds_sleepiness{4} = mean_corr_sleepiness_post_EC;
                mean_conds_sleepiness{5} = p_corrected_sleepiness{1};
                mean_conds_sleepiness{6} = p_corrected_sleepiness{2};
                mean_conds_sleepiness{7} = p_corrected_sleepiness{3};
                mean_conds_sleepiness{8} = p_corrected_sleepiness{4};

            end

        else            
            %% Fluid correlations
            % ------ Pre eyes open --------
            [corr_fluid_pre_EO, ~, p_fluid_pre_EO] = run_correlation(fluid_pre_EO, pre_eyes_open, signal_option);
                    
            % ------ Pre eyes closed --------        
            [corr_fluid_pre_EC, ~, p_fluid_pre_EC] = run_correlation(fluid_pre_EC, pre_eyes_closed, signal_option);
            
            % ----- Post eyes open --------        
            [corr_fluid_post_EO, ~, p_fluid_post_EO] = run_correlation(fluid_post_EO, post_eyes_open, signal_option);
           
            % ----- Post eyes closed -----        
            [corr_fluid_post_EC, ~, p_fluid_post_EC] = run_correlation(fluid_post_EC, post_eyes_closed, signal_option);
        
            % Correct for multiple comparisons
            p_corrected_fluid = p_correct(p_fluid_pre_EO, p_fluid_pre_EC, p_fluid_post_EO, p_fluid_post_EC);

            % list conds
            mean_conds_fluid = {};
            mean_conds_fluid{1} = p_corrected_fluid{1};
            mean_conds_fluid{2} = p_corrected_fluid{2};
            mean_conds_fluid{3} = p_corrected_fluid{3};
            mean_conds_fluid{4} = p_corrected_fluid{4};           

            %% Crystallized correlations
            % --------- Pre eyes open ------------
            [corr_crystallized_pre_EO, ~, p_crystallized_pre_EO]  = run_correlation(crystallized_pre_EO, pre_eyes_open, signal_option);
            
            % ------ Pre eyes closed --------        
            [corr_crystallized_pre_EC, ~, p_crystallized_pre_EC] = run_correlation(crystallized_pre_EC, pre_eyes_closed, signal_option);
            
            % ----- Post eyes open --------                
            [corr_crystallized_post_EO, ~, p_crystallized_post_EO] = run_correlation(crystallized_post_EO, post_eyes_open, signal_option);
            
            % ----- Post eyes closed -----        
            [corr_crystallized_post_EC, ~, p_crystallized_post_EC] = run_correlation(crystallized_post_EC, post_eyes_closed, signal_option);
            
            % Correct for multiple comparisons
            p_corrected_crystallized = p_correct(p_crystallized_pre_EO, p_crystallized_pre_EC, p_crystallized_post_EO, p_crystallized_post_EC);

            % list conds
            mean_conds_crystallized = {};
            mean_conds_crystallized{1} = p_corrected_crystallized{1};
            mean_conds_crystallized{2} = p_corrected_crystallized{2};
            mean_conds_crystallized{3} = p_corrected_crystallized{3};
            mean_conds_crystallized{4} = p_corrected_crystallized{4};

            if i_sample == 1
                %% Sleepiness
                % ------ Pre eyes open -------- 
                [corr_sleepiness_pre_EO, ~, p_sleepiness_pre_EO] = run_correlation(sleepiness_pre_EO, pre_eyes_open, signal_option);
            
                % ------ Pre eyes closed --------            
                [corr_sleepiness_pre_EC, ~, p_sleepiness_pre_EC] = run_correlation(sleepiness_pre_EC, pre_eyes_closed, signal_option);
                
                % ----- Post eyes open --------                    
                [corr_sleepiness_post_EO, ~, p_sleepiness_post_EO] = run_correlation(sleepiness_post_EO, post_eyes_open, signal_option);
                
                % ----- Post eyes closed -----            
                [corr_sleepiness_post_EC, ~, p_sleepiness_post_EC] = run_correlation(sleepiness_post_EC, post_eyes_closed, signal_option);
    
                % Correct for multiple comparisons
                p_corrected_sleepiness = p_correct(p_sleepiness_pre_EO, p_sleepiness_pre_EC, p_sleepiness_post_EO, p_sleepiness_post_EC);

                % list all conds
                mean_conds_sleepiness = {};
                mean_conds_sleepiness{1} = p_corrected_sleepiness{1};
                mean_conds_sleepiness{2} = p_corrected_sleepiness{2};
                mean_conds_sleepiness{3} = p_corrected_sleepiness{3};
                mean_conds_sleepiness{4} = p_corrected_sleepiness{4};

            end

        end

        % list conds   
        conds_fluid = {};
        conds_fluid{1} = corr_fluid_pre_EO;
        conds_fluid{2} = corr_fluid_pre_EC;
        conds_fluid{3} = corr_fluid_post_EO;
        conds_fluid{4} = corr_fluid_post_EC;

        % list all conds
        conds_crystallized = {};
        conds_crystallized{1} = corr_crystallized_pre_EO;
        conds_crystallized{2} = corr_crystallized_pre_EC;
        conds_crystallized{3} = corr_crystallized_post_EO;
        conds_crystallized{4} = corr_crystallized_post_EC;

        % list all conds
        conds_sleepiness = {};
        conds_sleepiness{1} = corr_sleepiness_pre_EO;
        conds_sleepiness{2} = corr_sleepiness_pre_EC;
        conds_sleepiness{3} = corr_sleepiness_post_EO;
        conds_sleepiness{4} = corr_sleepiness_post_EC;

        cond_labels = ["Pre-task EO", "Pre-task EC", ...
            "Post-task EO", "Post-task EC"];   
       
        %% Plot it out 
               
        if i_sample == 1

            outputfolder = "FullSample";

        else

            outputfolder = sample(i_sample);
            
        end

        %% Plots for each condition and frequency band
               
        % Fluid
        fluid_plot = make_topoplot(conds_fluid, mean_conds_fluid, cond_labels, i_signal, sorted_data);
    
        exportgraphics(fluid_plot, sprintf('../Results/Plots_%s/Fluid_Topoplot_%s.png', outputfolder, signal(i_signal)), 'Resolution', 1000)

        % Crystallized
        crystallized_plot = make_topoplot(conds_crystallized, mean_conds_crystallized, cond_labels, i_signal, sorted_data);
        
        exportgraphics(crystallized_plot, sprintf('../Results/Plots_%s/Crystallized_Topoplot_%s.png', outputfolder, signal(i_signal)), 'Resolution', 1000)
        
        if i_sample == 1 
            % Sleepiness
            sleepiness_plot = make_topoplot(conds_sleepiness, mean_conds_sleepiness, cond_labels, i_signal, sorted_data);
            
            exportgraphics(sleepiness_plot, sprintf('../Results/Plots_%s/Sleepiness_Topoplot_%s.png', outputfolder, signal(i_signal)), 'Resolution', 1000)               

        end
    
    end
end

disp("Done")
