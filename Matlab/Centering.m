

% Folder Setup
Datafolder = '../Data/';


% Load Lab File 
labs = readtable(fullfile(Datafolder, "Labs.txt"), "FileType","text",'Delimiter', '\t');
labs = labs(:,1:2);


% Load Sub List
subs_pre_EO  = readcell(fullfile(Datafolder, "subs_pre_EO.txt"));
subs_pre_EC  = readcell(fullfile(Datafolder, "subs_pre_EC.txt"));
subs_post_EO = readcell(fullfile(Datafolder, "subs_post_EO.txt"));
subs_post_EC = readcell(fullfile(Datafolder, "subs_post_EC.txt"));

% Load eeg_sorted_cond
for i_signal = ["Total", "Aperiodic/Parameters/Exponent", "Aperiodic/Parameters/Offset/", "Periodic"]
    for i_time = ["Pre", "Post"]
        for i_eye = ["EyesOpen", "EyesClosed"]            
          
            eeg_file = load(fullfile(Datafolder, "Ready_for_DDTBOX", "EEG_sorted_cond", "Full_Sample", i_signal, i_time, i_eye, "eeg_sorted_cond.mat"));
            eeg_file = eeg_file.eeg_sorted_cond;
            % 
            % Create subsets of Lab
            if strcmp(i_time, "Pre") && strcmp(i_eye, "EyesOpen")

                labs_subset = labs(ismember(labs.ID, subs_pre_EO), :);
                labs_subset.Index = (1:height(labs_subset))';

            elseif strcmp(i_time, "Pre") && strcmp(i_eye, "EyesClosed")

                labs_subset = labs(ismember(labs.ID, subs_pre_EC), :);
                labs_subset.Index = (1:height(labs_subset))';

            elseif strcmp(i_time, "Post") && strcmp(i_eye, "EyesOpen")

                labs_subset = labs(ismember(labs.ID, subs_post_EO), :);
                labs_subset.Index = (1:height(labs_subset))';

            elseif strcmp(i_time, "Post") && strcmp(i_eye, "EyesClosed")

                labs_subset = labs(ismember(labs.ID, subs_post_EC), :);
                labs_subset.Index = (1:height(labs_subset))';

            end
            % 
            % 
            % mean center files and save
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Bonn")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Bonn")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Bonn")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Dresden")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Dresden")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Dresden")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Giessen")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Giessen")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Giessen")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Riesel")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Riesel")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Riesel")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Wacker")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Wacker")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Hamburg_Wacker")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Köln")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Köln")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Köln")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Marburg")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Marburg")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Marburg")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Oldenburg")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Oldenburg")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Oldenburg")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Osnabrück")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Osnabrück")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Osnabrück")), 3);
            eeg_sorted_cond{1,1}(:,:,strcmp(labs_subset.Lab, "Würzburg")) = eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Würzburg")) - mean(eeg_file{1,1}(:,:,strcmp(labs_subset.Lab, "Würzburg")), 3);

            if ~(isfolder(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "EEG_sorted_cond", i_signal, i_time, i_eye)))
                mkdir(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "EEG_sorted_cond", i_signal, i_time, i_eye))
            end

            save(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "EEG_sorted_cond", i_signal, i_time, i_eye, "eeg_sorted_cond.mat"), "eeg_sorted_cond")

            clear eeg_sorted_cond

            % Loop through Traits and mean center

            
            for i_behav = ["Crystallized", "Fluid", "Pre_Sleepiness", "Post_Sleepiness"]
                
                if i_time == "Pre" && i_behav == "Post_Sleepiness"
                    continue
                elseif i_time == "Post" && i_behav == "Pre_Sleepiness"
                    continue
                else
    
                    beh_file = load(fullfile(Datafolder, "Ready_for_DDTBOX", "Behavioral_scores", "Full_Sample", i_behav, i_time, i_eye, "eeg_sorted_cond_regress_sorted_cond.mat"));
                    beh_file = beh_file.SVR_labels;
        
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Bonn")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Bonn")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Bonn")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Dresden")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Dresden")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Dresden")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Giessen")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Giessen")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Giessen")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Hamburg_Riesel")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Hamburg_Riesel")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Hamburg_Riesel")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Hamburg_Wacker")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Hamburg_Wacker")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Hamburg_Wacker")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Köln")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Köln")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Köln")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Marburg")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Marburg")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Marburg")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Oldenburg")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Oldenburg")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Oldenburg")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Osnabrück")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Osnabrück")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Osnabrück")),"omitnan");
                    SVR_labels{1,1}(strcmp(labs_subset.Lab, "Würzburg")) = beh_file{1,1}(strcmp(labs_subset.Lab, "Würzburg")) - mean(beh_file{1,1}(strcmp(labs_subset.Lab, "Würzburg")),"omitnan");
                    
                    if ~(isfolder(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "Behavioral_scores", "Full_Sample", i_behav, i_time, i_eye)))
                        mkdir(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "Behavioral_scores", "Full_Sample", i_behav, i_time, i_eye))
                    end
                    
                    SVR_labels{1,1} = SVR_labels{1,1}';
    
                    save(fullfile(Datafolder, "Ready_for_DDTBOX_centered", "Behavioral_scores", "Full_Sample", i_behav, i_time, i_eye, "eeg_sorted_cond_regress_sorted_cond.mat"), "SVR_labels")
        
                    clear SVR_labels

                end

            end % for i_behav                   
            
        end % for i_eye
    
    end % for i_time

end % for i_signal

