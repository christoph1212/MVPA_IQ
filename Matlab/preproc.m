function preproc(NrSubsets, ThisSubset)
% This function preprocesses the raw EEG data and saves the power spectrum
% in a .csv file within a parfor-loop. The raw data should be a .set file. 
% Change the information about Basefolder and Rootfolder according to your
% directories.
%
% Input:
% (1) NrSubsets = ###
% (2) ThisSubset = ###
%
% Output:
% (1) .csv file with the spectral power data for each participant from 0.5 to
% 30 Hz and each channel
% (2) Log file for each data file
% (3) Error file 
%
% This function was written by Christoph Frühlinger and Dr. Katharina Paul,
% last edits: February 2024

fprintf("Input Preproc- NrSubsets %s, ThisSubset %s", NrSubsets , ThisSubset)
NrSubsets = str2num(NrSubsets);
ThisSubset = str2num(ThisSubset);

% Set and create Folders
%Basefolder = '/work/bay2875/MVPA_Personality/';
Basefolder = '../';
%Rootfolder = '/work/bay2875/RawData/task-Resting/';
Rootfolder = '../Data/RawData';
Datafolder = fullfile(Basefolder, 'Data', 'Preprocessed'); 
DatafolderPreEyesClosed = fullfile(Datafolder, 'Pre', 'EyesClosed');
DatafolderPostEyesClosed = fullfile(Datafolder, 'Post', 'EyesClosed');
DatafolderPreEyesOpen = fullfile(Datafolder, 'Pre', 'EyesOpen');
DatafolderPostEyesOpen= fullfile(Datafolder, 'Post', 'EyesOpen');
LogFolder = fullfile(Datafolder, 'LogFiles');
ErrorFolder = fullfile(Datafolder, 'ErrorFiles');
if ~(isfolder(DatafolderPreEyesClosed))                
    mkdir(DatafolderPreEyesClosed)
end
if ~(isfolder(DatafolderPostEyesClosed))                
    mkdir(DatafolderPostEyesClosed)
end
if ~(isfolder(DatafolderPreEyesOpen))                
    mkdir(DatafolderPreEyesOpen)
end
if ~(isfolder(DatafolderPostEyesOpen))                
    mkdir(DatafolderPostEyesOpen)
end


if ~(isfolder(LogFolder))
    mkdir(LogFolder)
end
if ~(isfolder(ErrorFolder))
    mkdir(ErrorFolder)
end



% Get Participant IDs
Participants = dir(Rootfolder);
Participants = Participants(~ismember({Participants.name}, {'.', '..'}));
Participants = {Participants.name};     

% Make parallel across subjects
SubsetSize = ceil(length(Participants)/NrSubsets);
IndexSubs = (1:SubsetSize)+((ThisSubset-1)*SubsetSize);
IndexSubs = IndexSubs(IndexSubs<=length(Participants));
Participants = Participants(IndexSubs);

% EEG Struct for Interpolation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG_interpolation = pop_loadset('filename', ['sub-AG04EN28_task-Resting_run-1_eeg.set'], ...
    'filepath', fullfile(Rootfolder, 'sub-AG04EN28', 'eeg'));

Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', 'F1', 'F2', ...
                'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'FT7', 'FT8', 'FC5', 'FC6', 'FC3', ...
                'FC4', 'FC1', 'FC2', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'T7', 'T8', ...
                'TP7', 'TP8', 'CP5', 'CP6', 'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', ...
                'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', ...
                'O2', 'OZ', 'POZ', 'PZ', 'CPZ', 'CZ', 'FCZ', 'FZ', 'VOGabove', ...
                'VOGbelow', 'HOGl', 'HOGr'};

EEG_Channels_only = Common_Channels(ismember(Common_Channels, {EEG_interpolation.chanlocs(strcmp({EEG_interpolation.chanlocs.type}, 'EEG')).labels}));
EEG_Channels_only = EEG_interpolation.chanlocs(ismember({EEG_interpolation.chanlocs.labels}, EEG_Channels_only));


poolobj = gcp('nocreate');
delete(poolobj);
parpool()

% Loop through Participants
parfor subj = 1:ceil(length(Participants))

    i_participant = Participants{subj}; 
    % Loop through Sessions
    for i_run = 1:2
        try
            % test if not completed          
            if isfile(fullfile(DatafolderPostEyesClosed, [i_participant '_post_Average_EC.csv'])) && ...
                isfile(fullfile(DatafolderPreEyesClosed, [i_participant '_pre_Average_EC.csv'])) && ...
                isfile(fullfile(DatafolderPostEyesOpen, [i_participant '_post_Average_EO.csv'])) && ...
                isfile(fullfile(DatafolderPreEyesOpen, [i_participant '_pre_Average_EO.csv']))       
                   continue
            end

            
            % Log file setup
            Log = struct(...
                'Flatlined', [],...
                'Bad_Channels', [], ...
                'Clean_Segments_EO_T1', [], ...
                'Clean_Segments_EO_T2', [], ...
                'Clean_Segments_EC_T1', [], ...
                'Clean_Segments_EC_T2', [], ...
                'ICA_Total_Components_EO', [], ...
                'ICA_Bad_Components_EO', [], ...
                'ICA_Total_Components_EC', [], ...
                'ICA_Bad_Components_EC', [], ...
                'Interpolation', [], ...
                'Epochs_EO', [], ...
                'Artefact1_EO', [], ...
                'Artefact2_EO', [], ...
                'Artefact3_EO', [], ...
                'Epochs_EC', [], ...
                'Artefact1_EC', [], ...
                'Artefact2_EC', [], ...
                'Artefact3_EC', []);

            %% 1. Load EEGLAB and File
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            EEG = pop_loadset('filename', [i_participant '_task-Resting_run-' num2str(i_run) '_eeg.set'], ...
                'filepath', fullfile(Rootfolder, i_participant, 'eeg'));
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
            
            %% 2. Add Reference and Rereference to CZ                       
            if ~strcmp(EEG.Info_Lab.Reference,'CMS/DRL') % for Biosemi we do not recover CMS/DRL   
                    % if Reference is not already in the data, add 0 filled Channel
                    if ~ismember( upper(EEG.Info_Lab.Reference), upper({EEG.chanlocs.labels}))
                        EEG.data(end+1,:) = 0; % add 0 filled Channel
                        EEG.nbchan = size(EEG.data,1); % update Channel Number
            
                        % Chanloc information on reference electrode
                        if strcmp(EEG.Info_Lab.Reference,'FCZ')
                            RefInfo = {'FCZ',	 -0,    0.127,	0.388,	  0, 0.922,	  0,	67.2,  1,	'EEG'};
                        elseif  strcmp(EEG.Info_Lab.Reference, 'CZ')
                            RefInfo = {'CZ',	 -0,     0,	    6.12e-17,  0,  1,	  0,  90,  1,	'EEG'};
                        end
                        % Loop through Chanloc fields (labels, theta, radius etc.) and add
                        % new info
                        ChanlocFields = fields(EEG.chanlocs);
                        for ifield = 1:10
                            EEG.chanlocs(EEG.nbchan).(ChanlocFields{ifield}) = RefInfo{ifield};
                        end
                        % update EEGlab structure
                        EEG = eeg_checkset( EEG );
                    end
            end
            
            EEG = pop_reref( EEG, 'CZ');
            EEG = eeg_checkset( EEG );
            
            %% 3. Apply High-Pass, Low-Pass and Notch-Filter            
            EEG = pop_eegfiltnew(EEG, 'locutoff',0.5);
            EEG = pop_eegfiltnew(EEG, 'hicutoff',30);
            
            %% 4. Trim Channels from Dataset            
            Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', 'F1', 'F2', ...
                'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'FT7', 'FT8', 'FC5', 'FC6', 'FC3', ...
                'FC4', 'FC1', 'FC2', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'T7', 'T8', ...
                'TP7', 'TP8', 'CP5', 'CP6', 'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', ...
                'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', ...
                'O2', 'OZ', 'POZ', 'PZ', 'CPZ', 'CZ', 'FCZ', 'FZ', 'VOGabove', ...
                'VOGbelow', 'HOGl', 'HOGr'}; % 'AFZ', 'FPZ' are used as grounds in some labs
            Common_Channels = Common_Channels(ismember(Common_Channels, {EEG.chanlocs.labels})); % some labs miss e.g. VOGabove
            EEG = pop_select( EEG, 'channel',Common_Channels); 
            
            %% 5. Flat and Bad Channel Detection with FASTER            
            [EEG, ~, BUR] = clean_artifacts(EEG, 'FlatlineCriterion',5, ...
                'ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
                'BurstCriterion','off','WindowCriterion','off','BurstRejection','off', ...
                'Distance','Euclidian');
            flatlined = (BUR.nbchan - EEG.nbchan);

	        flatlinedch = find(~contains({BUR.chanlocs.labels}, {EEG.chanlocs.labels}));
            Log.Flatlined = [strjoin({BUR.chanlocs(flatlinedch).labels}, ', ')];

            
            Bad_Channel_Mask = zeros(EEG.nbchan, 1);
            BadChannels_Index_2 = [];
            Faster_estimate = channel_properties(EEG, 1:EEG.nbchan, find({EEG.chanlocs.labels} == "FZ"));
            BadChannels_Index_2 = find(min_z(Faster_estimate))';
            Bad_Channel_Mask(BadChannels_Index_2) = 1;

            Log.Bad_Channels = [strjoin({EEG.chanlocs(find(Bad_Channel_Mask)).labels}, ', ')];

            EEG = pop_select(EEG, 'channel', find(~Bad_Channel_Mask));
            EEG = eeg_checkset( EEG );
            
            
            % Go to next dataset if more than 4 channels are exluded
            if flatlined + sum(Bad_Channel_Mask) > 4
                warning('Too many channels faulty. Dataset was excluded.')
                savetofile(Log, fullfile(LogFolder, [i_participant '_run_' num2str(i_run) '_Log_excluded.mat'])) 
                %continue
            end
            
            %% 6. Separate into Eyes Open and Eyes Closed
            % Start shortly after condition onset to reduce risk of artefacts and add 60 seconds
	        Events_noBoundary = EEG.event(~strcmp({EEG.event.type}, 'boundary') );
            p = [[Events_noBoundary.latency] + EEG.srate/2; [Events_noBoundary.latency] + EEG.srate*60; Events_noBoundary.Condition]';
            points_EO = str2double([p(p(:,3) == 'Open', 1:2)]);
            points_EC = str2double([p(p(:,3) == 'Closed', 1:2)]);
            
            EEG_EO = pop_select( EEG, 'point',points_EO );
            EEG_EC = pop_select( EEG, 'point',points_EC );
            
            EEG_EO = eeg_checkset( EEG_EO );
            EEG_EC = eeg_checkset( EEG_EC );
            
            for i_cond = ["EO", "EC"]
                if i_cond == "EO"
                    EEG_work = EEG_EO;
                else 
                    EEG_work = EEG_EC;
                end
                
                try
                %% 7. Remove Bad Data with ASR                
                Clean_Segment_Mask = ones(EEG_work.pnts,1);
                Log.(strcat('Clean_Segments_', i_cond,'_T1')) = length(Clean_Segment_Mask); % KP Why is this here? Sum here equals the length, also later Why calling it clean if it is actually the total Nr?
                EEG_Channels = ismember({EEG_work.chanlocs.type}, {'EEG'});
                EEG_subset = pop_select( EEG_work, 'channel', find(EEG_Channels));
                EEG_fil = pop_eegfiltnew(EEG_subset, 'locutoff', 1, 'hicutoff', []);
                EEG_cleaned = pop_clean_rawdata(EEG_fil, 'FlatlineCriterion', 'off', ...
                             'ChannelCriterion','off','LineNoiseCriterion','off', 'Highpass','off', ...
                             'BurstCriterion',50,'WindowCriterion','off','BurstRejection','on', ...
                             'Distance','Euclidian');
                Clean_Segment_Mask = EEG_cleaned.etc.clean_sample_mask';
                retain_data_intervals = reshape(find(diff([false Clean_Segment_Mask' false])),2,[])';
                retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
                EEG_work = pop_select(EEG_work, 'point', retain_data_intervals);
                EEG_work = eeg_checkset( EEG_work );
                
                Log.(strcat('Clean_Segments_', i_cond,'_T2')) = sum(Clean_Segment_Mask);
            
                %% 8. Run ICA                    
                EEG_work_ica = pop_eegfiltnew(EEG_work, 'locutoff',1,'plotfreqz',0);
                rank_used = EEG.nbchan-1;
                rank_used = min([rank_used, rank(double(EEG_work_ica.data(:,:)))]);
                EEG_work_ica = pop_runica(EEG_work_ica,'icatype','runica','pca' ,rank_used, ...
                    'sphering', 'on');
                EEG_work.icaweights = EEG_work_ica.icaweights;
                EEG_work.icasphere = EEG_work_ica.icasphere;
                EEG_work = eeg_checkset(EEG_work, 'ica');
        
                %% 9. Remove Occular Components                    
                Temp_work = iclabel(EEG_work, 'default');
                badIC_Info = Temp_work.etc.ic_classification.ICLabel;
                badIC = find(badIC_Info.classifications(:,3)>0.7 & badIC_Info.classifications(:,1)<0.7);
                EEG_work.icaact = eeg_getica(EEG_work);
                EEG_work = pop_subcomp(EEG_work, badIC, 0);
                
                Log.(strcat('ICA_Total_Components_', i_cond)) = [size(badIC_Info.classifications, 1)];
                Log.(strcat('ICA_Bad_Components_', i_cond)) = [size(badIC_Info.classifications, 1) - size(EEG_work.icaweights, 1)];
        
                %% 10. Add previous Reference back and Rereference to common average                    
                RefInfo = {'CZ',	 -0,     0,	   0,  0,  1,	  0,  90,  1,	'EEG'};
                
                ChanlocFields = fields(EEG_work.chanlocs);
                for ifield = 1:10
                    EEG_work.chanlocs(EEG_work.nbchan+1).(ChanlocFields{ifield}) = RefInfo{ifield};
                end
                EEG_work.data(end+1,:) = 0;
                EEG_work.nbchan = size(EEG_work.data,1);
                EEG_work = eeg_checkset( EEG_work );
        
                EEG_work = pop_select( EEG_work, 'rmchantype',{'EOG'});            
                EEG_work = pop_reref( EEG_work, []);
                EEG_work = eeg_checkset( EEG_work );

                %% 11. Channel Interpolation        
                interpolate_channels = EEG_Channels_only(~ismember({EEG_Channels_only.labels}, {EEG_work.chanlocs.labels}));
                EEG_work = pop_interp(EEG_work, interpolate_channels, 'spherical');
                EEG_work = eeg_checkset( EEG_work );
                Log.Interpolation = [strjoin({interpolate_channels(:).labels}, ', ')];
        
                %% 12. Extract 2s Epochs with 1s Overlap                    
                EEG_work = eeg_regepochs(EEG_work, 1, [0 2], 0, 'X', 'on');          
                EEG_work = eeg_checkset( EEG_work );

                %% 13. Artefact Rejection   
                % Count as Artefact if:
                % - Voltage difference is >50µV or <-50µV every ms (here every 2ms) 
                %   --> 1st if Condition.
                % - Voltage difference is >200µV or <-200µV in a 100ms intervall 
                %   --> 2nd if Condition
                % - Voltage difference between min/max Voltage is < 0.5µV in a 100ms
                %   intervall
                %   --> 3rd if Condition
        
                Log.(strcat('Epochs_', i_cond)) = EEG_work.trials;
                % Artefact 1
                BadEpoch = [];
                BadEpochs_Mask = zeros(EEG_work.trials,1);                                       
                for itrial = 1:(EEG_work.trials)
                    if sum(sum(abs(diff(EEG_work.data(:,:,itrial),1,2)) > 50))
                        BadEpoch = [BadEpoch, itrial]; 
                    end
                end
                BadEpoch = unique(BadEpoch);
                BadEpochs_Mask(BadEpoch) = 1;
                Log.(strcat('Artefact1_', i_cond)) = BadEpoch;
                EEG_work = pop_select(EEG_work, 'trial', find(~BadEpochs_Mask));

                % Aftefact 2
                BadEpoch = [];
                BadEpochs_Mask = zeros(EEG_work.trials,1);                    
                window_length = EEG_work.srate*0.1;
                for itrial = 1:(EEG_work.trials)
                    for iIncr = 1:(EEG_work.pnts - window_length) 
                        if any(abs(max(EEG_work.data(:,iIncr:iIncr + window_length,itrial)) - min(EEG_work.data(:,iIncr:iIncr + window_length,itrial))) > 200)
                            BadEpoch = [BadEpoch, itrial];
                        end
                    end
                end
                BadEpoch = unique(BadEpoch);
                BadEpochs_Mask(BadEpoch) = 1;
                Log.(strcat('Artefact2_', i_cond)) = BadEpoch;
                EEG_work = pop_select(EEG_work, 'trial', find(~BadEpochs_Mask));

                % Artefact 3            
                BadEpoch = [];
                BadEpochs_Mask = zeros(EEG_work.trials,1);                    
                for itrial = 1:(EEG_work.trials)
                    for iIncr1 = 1:(EEG_work.pnts-window_length)
                        if any(abs(max(EEG_work.data(:,iIncr1:iIncr1 + window_length,itrial)) - min(EEG_work.data(:,iIncr1:iIncr1 + window_length, itrial))) < 0.5)
                            BadEpoch = [BadEpoch, itrial]; 
                        end
                    end
                end
                BadEpoch = unique(BadEpoch);
                BadEpochs_Mask(BadEpoch) = 1;
                Log.(strcat('Artefact3_', i_cond)) = BadEpoch;
                EEG_work = pop_select(EEG_work, 'trial', find(~BadEpochs_Mask));                                   
                   
                %% 14. FFT and save Files                
                nfft = 2^nextpow2(EEG_work.pnts);
                NrFreqs = nfft/2+1;
                FreqData = zeros(length({EEG_work.chanlocs.labels}), nfft, EEG_work.trials);
                Data = EEG_work.data();
                Data = Data.*[hamming(EEG_work.pnts)]';
                power = NaN(EEG_work.nbchan, NrFreqs, EEG_work.trials);
                
                for n = 1:EEG_work.nbchan
                    for itrial = 1:EEG_work.trials
                        FreqData(n, :, itrial) = fft(Data(n, :, itrial),nfft)/EEG_work.pnts;
                    end
                end

                power(:,1:NrFreqs,1:EEG_work.trials) = abs(FreqData(:,1:NrFreqs, :)).^2;
                Average = horzcat({EEG_work.chanlocs.labels}', num2cell(mean(power(:,1:60,:),3)));
        
                if i_cond == "EO"
                    if i_run == 1
                        writecell(Average, fullfile(DatafolderPreEyesOpen, [i_participant '_pre_Average_EO.csv'])); 
                    else
                        writecell(Average, fullfile(DatafolderPostEyesOpen, [i_participant '_post_Average_EO.csv']));
                    end
                else
                    if i_run == 1
                        writecell(Average, fullfile(DatafolderPreEyesClosed, [i_participant '_pre_Average_EC.csv'])); 
                    else
                        writecell(Average, fullfile(DatafolderPostEyesClosed, [i_participant '_post_Average_EC.csv']));
                    end
                end

                %% Error Management
                catch error2
                    ErrorMessage = string(error2.message);
                    for ierrors = 1:length(error2.stack)
                        ErrorMessage = strcat(ErrorMessage, '//', num2str(error2.stack(ierrors).name), ', Line: ', num2str(error2.stack(ierrors).line));
                        writematrix(ErrorMessage, fullfile(ErrorFolder, [i_participant '_run_' num2str(i_run) '_eyes_' char(i_cond) '.csv']))
                    end
                end

            end % for condition         
            savetofile(Log, fullfile(LogFolder, [i_participant '_run_' num2str(i_run) '_Log_power_butter_no_notch.mat']));

        catch error
            ErrorMessage = string(error.message);
            for ierrors = 1:length(error.stack)
                ErrorMessage = strcat(ErrorMessage, '//', num2str(error.stack(ierrors).name), ', Line: ', num2str(error.stack(ierrors).line));
            end
            writematrix(ErrorMessage, fullfile(ErrorFolder, [i_participant '_run_' num2str(i_run) '.csv']))
            savetofile(Log, fullfile(LogFolder, [i_participant '_run_' num2str(i_run) '_Log.mat']))
        end        
    end % for run
end % for parfor
