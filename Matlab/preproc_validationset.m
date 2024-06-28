function preproc_validationset(NrSubsets, ThisSubset)
% This function preprocesses the raw EEG data and saves the power spectrum
% in a .csv file within a parfor-loop. The raw data should be a .bdf file
% and stored in a folder called "ValidationSet".Change the information 
% about Basefolder and Rootfolder according to your directories. This
% script is similar to the other preprocessing function, but adapted for
% the Validation Set from Ociepka et al. (2022).
%
% Input:
% (1) NrSubsets = default 1
% (2) ThisSubset = default 1
%
% Output:
% (1) .csv file with the spectral power data for each participant from 0.5 to
% 30 Hz and each channel
% (2) Log file for each data file
% (3) Error file 
%
% This function was written by Christoph Frühlinger and Dr. Katharina Paul,
% and modified in: June 2024

fprintf("Input Preproc- NrSubsets %s, ThisSubset %s", NrSubsets , ThisSubset)
NrSubsets = str2num(NrSubsets);
ThisSubset = str2num(ThisSubset);

% Set and create Folders
Basefolder = '../';
Rootfolder = '../ValidationSet';
Datafolder = fullfile(Basefolder, 'DataValidation', 'Preprocessed'); 
DatafolderPreEyesOpen = fullfile(Datafolder, 'Pre', 'EyesOpen');
DatafolderPostEyesOpen= fullfile(Datafolder, 'Post', 'EyesOpen');
LogFolder = fullfile(Datafolder, 'LogFiles');
ErrorFolder = fullfile(Datafolder, 'ErrorFiles');

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
Participants = readtable(fullfile(Rootfolder, 'IAF_and_beh_data.csv'));
Participants = Participants.ID;   

% Make parallel across subjects
SubsetSize = ceil(length(Participants)/NrSubsets);
IndexSubs = (1:SubsetSize)+((ThisSubset-1)*SubsetSize);
IndexSubs = IndexSubs(IndexSubs<=length(Participants));
Participants = Participants(IndexSubs);

% EEG Struct for Interpolation
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG_interpolation = pop_biosig(fullfile(Rootfolder, '1. Rest', 'rest1_1.bdf'));
EEG_interpolation = pop_chanedit(EEG_interpolation, 'lookup', 'Biosemi_72.ced');

Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', 'F1', 'F2', ...
                'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'FT7', 'FT8', 'FC5', 'FC6', 'FC3', ...
                'FC4', 'FC1', 'FC2', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'T7', 'T8', ...
                'TP7', 'TP8', 'CP5', 'CP6', 'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', ...
                'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', ...
                'O2', 'OZ', 'POZ', 'PZ', 'CPZ', 'CZ', 'FCZ', 'FZ', 'VOGabove', ...
                'VOGbelow', 'HOGl', 'HOGr'};

EEG_Channels_only = Common_Channels(ismember(Common_Channels, {EEG_interpolation.chanlocs(strcmp({EEG_interpolation.chanlocs.type}, 'EEG')).labels}));
EEG_Channels_only = {EEG_interpolation.chanlocs(strcmp({EEG_interpolation.chanlocs.type}, 'EEG')).labels};
EEG_Channels_only = EEG_interpolation.chanlocs(ismember({EEG_interpolation.chanlocs.labels}, EEG_Channels_only));


poolobj = gcp('nocreate');
delete(poolobj);
parpool('local')

% Loop through Participants
parfor subj = 1:length(Participants)

    i_participant = Participants(subj); 
    % Loop through Sessions
    for i_run = 1:2
        
        try
            % test if not completed          
            if isfile(fullfile(DatafolderPostEyesOpen, ['rest_' num2str(i_participant) '_post_Average_EO.csv'])) && ...
                isfile(fullfile(DatafolderPreEyesOpen, ['rest_' num2str(i_participant) '_pre_Average_EO.csv']))       
                   continue
            end
            
            % Log file setup
            Log = struct(...
                'Flatlined', [],...
                'Bad_Channels', [], ...
                'Clean_Segments_T1', [], ...
                'Clean_Segments_T2', [], ...
                'ICA_Total_Components', [], ...
                'ICA_Bad_Components', [], ...
                'Interpolation', [], ...
                'Epochs', [], ...
                'Artefact1', [], ...
                'Artefact2', [], ...
                'Artefact3', []);

            %% 1. Load EEGLAB and File
            if i_run == 1
                RawDatafolder = fullfile(Rootfolder, '1. Rest');
            else
                RawDatafolder = fullfile(Rootfolder, '4. Rest');
            end
            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
            EEG = pop_biosig(fullfile(RawDatafolder, ['rest' num2str(i_run) '_' num2str(i_participant) '.bdf']));
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);

            EEG.chanlocs(65).labels = 'HOGl';
            EEG.chanlocs(66).labels = 'HOGr';
            EEG.chanlocs(67).labels = 'VOGabove';
            EEG.chanlocs(68).labels = 'VOGbelow';
            EEG.chanlocs(69).labels = 'MASTl';
            EEG.chanlocs(70).labels = 'MASTr';
            
            EEG = pop_chanedit(EEG, 'lookup', 'Biosemi_72.ced');
            
            %% 2. Add Reference and Rereference to CZ 
            % Reference in the Validation Set was 'common'. Thus adding
            % reference is skipped.

            % if ~strcmp(EEG.Info_Lab.Reference,'CMS/DRL') % for Biosemi we do not recover CMS/DRL   
            %         % if Reference is not already in the data, add 0 filled Channel
            %         if ~ismember( upper(EEG.Info_Lab.Reference), upper({EEG.chanlocs.labels}))
            %             EEG.data(end+1,:) = 0; % add 0 filled Channel
            %             EEG.nbchan = size(EEG.data,1); % update Channel Number
            % 
            %             % Chanloc information on reference electrode
            %             if strcmp(EEG.Info_Lab.Reference,'FCz')
            %                 RefInfo = {'FCz',	 -0,    0.127,	0.388,	  0, 0.922,	  0,	67.2,  1,	'EEG'};
            %             elseif  strcmp(EEG.Info_Lab.Reference, 'Cz')
            %                 RefInfo = {'Cz',	 -0,     0,	    6.12e-17,  0,  1,	  0,  90,  1,	'EEG'};
            %             end
            %             % Loop through Chanloc fields (labels, theta, radius etc.) and add
            %             % new info
            %             ChanlocFields = fields(EEG.chanlocs);
            %             for ifield = 1:10
            %                 EEG.chanlocs(EEG.nbchan).(ChanlocFields{ifield}) = RefInfo{ifield};
            %             end
            %             % update EEGlab structure
            %             EEG = eeg_checkset( EEG );
            %         end
            % end
            
            EEG = pop_reref( EEG, 'Cz');
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
                'O2', 'Oz', 'POz', 'Pz', 'CPz', 'Cz', 'FCz', 'Fz', 'VOGabove', ...
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
            Faster_estimate = channel_properties(EEG, 1:EEG.nbchan, find({EEG.chanlocs.labels} == "Fz"));
            BadChannels_Index_2 = find(min_z(Faster_estimate))';
            Bad_Channel_Mask(BadChannels_Index_2) = 1;

            Log.Bad_Channels = [strjoin({EEG.chanlocs(find(Bad_Channel_Mask)).labels}, ', ')];

            EEG = pop_select(EEG, 'channel', find(~Bad_Channel_Mask));
            EEG = eeg_checkset( EEG );
            
            
            % Go to next dataset if more than 4 channels are exluded
            if flatlined + sum(Bad_Channel_Mask) > 4
                warning('Too many channels faulty. Dataset was excluded.')
                savetofile(Log, fullfile(LogFolder, ['rest' num2str(i_run) '_' num2str(i_participant) '_Log_excluded.mat'])) 
                continue
            end
            
            %% 6. Separate into Eyes Open and Eyes Closed
            
	        if sum(ismember([EEG.event.type], 64513)) == 1
                eventID = find([EEG.event.type] == 64513);
                points = [EEG.event(eventID).latency, EEG.event(eventID+1).latency];       
                EEG = pop_select(EEG, 'point', points);
                EEG = eeg_checkset( EEG );
            elseif sum(ismember([EEG.event.type], 65281)) == 1
                eventID = find([EEG.event.type] == 64513)
                points = [EEG.event(eventID).latency, EEG.event(eventID+1).latency];       
                EEG = pop_select(EEG, 'point', points);
                EEG = eeg_checkset( EEG );
            else
                warning('Could not find event triggers!')
            end
                        
            %% 7. Remove Bad Data with ASR                
            Clean_Segment_Mask = ones(EEG.pnts,1);
            Log.Clean_Segments_T1 = length(Clean_Segment_Mask); 
            EEG_Channels = ismember({EEG.chanlocs.type}, {'EEG'});
            EEG_subset = pop_select( EEG, 'channel', find(EEG_Channels));
            EEG_fil = pop_eegfiltnew(EEG_subset, 'locutoff', 1, 'hicutoff', []);
            EEG_cleaned = pop_clean_rawdata(EEG_fil, 'FlatlineCriterion', 'off', ...
                         'ChannelCriterion','off','LineNoiseCriterion','off', 'Highpass','off', ...
                         'BurstCriterion',50,'WindowCriterion','off','BurstRejection','on', ...
                         'Distance','Euclidian');
            Clean_Segment_Mask = EEG_cleaned.etc.clean_sample_mask';
            retain_data_intervals = reshape(find(diff([false Clean_Segment_Mask' false])),2,[])';
            retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
            EEG = pop_select(EEG, 'point', retain_data_intervals);
            EEG = eeg_checkset( EEG );
            
            Log.Clean_Segments_T2 = sum(Clean_Segment_Mask);
        
            %% 8. Run ICA                    
            EEG_ica = pop_eegfiltnew(EEG, 'locutoff',1,'plotfreqz',0);
            rank_used = EEG.nbchan-1;
            rank_used = min([rank_used, rank(double(EEG_ica.data(:,:)))]);
            EEG_ica = pop_runica(EEG_ica,'icatype','runica','pca' ,rank_used, ...
                'sphering', 'on');
            EEG.icaweights = EEG_ica.icaweights;
            EEG.icasphere = EEG_ica.icasphere;
            EEG = eeg_checkset(EEG, 'ica');
    
            %% 9. Remove Occular Components                    
            Temp_work = iclabel(EEG, 'default');
            badIC_Info = Temp_work.etc.ic_classification.ICLabel;
            badIC = find(badIC_Info.classifications(:,3)>0.7 & badIC_Info.classifications(:,1)<0.7);
            EEG.icaact = eeg_getica(EEG);
            EEG = pop_subcomp(EEG, badIC, 0);
            
            Log.ICA_Total_Components = [size(badIC_Info.classifications, 1)];
            Log.ICA_Bad_Components = [size(badIC_Info.classifications, 1) - size(EEG.icaweights, 1)];
    
            %% 10. Add previous Reference back and Rereference to common average                    
            RefInfo = {'Cz', 'EEG',	 -0,     0,	   0,  0,  1,	  0,  90,  1};
            
            ChanlocFields = fields(EEG.chanlocs);
            for ifield = 1:10
                EEG.chanlocs(EEG.nbchan+1).(ChanlocFields{ifield}) = RefInfo{ifield};
            end
            EEG.data(end+1,:) = 0;
            EEG.nbchan = size(EEG.data,1);
            EEG = eeg_checkset( EEG );
    
            EEG = pop_select( EEG, 'rmchantype',{'EOG'});            
            EEG = pop_reref( EEG, []);
            EEG = eeg_checkset( EEG );

            %% 11. Channel Interpolation        
            interpolate_channels = EEG_Channels_only(~ismember({EEG_Channels_only.labels}, {EEG.chanlocs.labels}));
            EEG = pop_interp(EEG, interpolate_channels, 'spherical');
            EEG = eeg_checkset( EEG );
            Log.Interpolation = [strjoin({interpolate_channels(:).labels}, ', ')];
    
            %% 12. Extract 2s Epochs with 1s Overlap                    
            EEG = eeg_regepochs(EEG, 1, [0 2], 0, 'X', 'on');          
            EEG = eeg_checkset( EEG );

            %% 13. Artefact Rejection   
            % Count as Artefact if:
            % - Voltage difference is >50µV or <-50µV every ms (here every 2ms) 
            %   --> 1st if Condition.
            % - Voltage difference is >200µV or <-200µV in a 100ms intervall 
            %   --> 2nd if Condition
            % - Voltage difference between min/max Voltage is < 0.5µV in a 100ms
            %   intervall
            %   --> 3rd if Condition
    
            Log.Epochs = EEG.trials;
            % Artefact 1
            BadEpoch = [];
            BadEpochs_Mask = zeros(EEG.trials,1);                                       
            for itrial = 1:(EEG.trials)
                if sum(sum(abs(diff(EEG.data(:,:,itrial),1,2)) > 50))
                    BadEpoch = [BadEpoch, itrial]; 
                end
            end
            BadEpoch = unique(BadEpoch);
            BadEpochs_Mask(BadEpoch) = 1;
            Log.Artefact1 = BadEpoch;
            EEG = pop_select(EEG, 'trial', find(~BadEpochs_Mask));

            % Aftefact 2
            BadEpoch = [];
            BadEpochs_Mask = zeros(EEG.trials,1);                    
            window_length = round(EEG.srate*0.1);
            for itrial = 1:(EEG.trials)
                for iIncr = 1:(EEG.pnts - window_length) 
                    if any(abs(max(EEG.data(:,iIncr:iIncr + window_length,itrial)) - min(EEG.data(:,iIncr:iIncr + window_length,itrial))) > 200)
                        BadEpoch = [BadEpoch, itrial];
                    end
                end
            end
            BadEpoch = unique(BadEpoch);
            BadEpochs_Mask(BadEpoch) = 1;
            Log.Artefact2 = BadEpoch;
            EEG = pop_select(EEG, 'trial', find(~BadEpochs_Mask));

            % Artefact 3            
            BadEpoch = [];
            BadEpochs_Mask = zeros(EEG.trials,1);                    
            for itrial = 1:(EEG.trials)
                for iIncr1 = 1:(EEG.pnts-window_length)
                    if any(abs(max(EEG.data(:,iIncr1:iIncr1 + window_length,itrial)) - min(EEG.data(:,iIncr1:iIncr1 + window_length, itrial))) < 0.5)
                        BadEpoch = [BadEpoch, itrial]; 
                    end
                end
            end
            BadEpoch = unique(BadEpoch);
            BadEpochs_Mask(BadEpoch) = 1;
            Log.Artefact3 = BadEpoch;
            EEG = pop_select(EEG, 'trial', find(~BadEpochs_Mask));                                   
               
            %% 14. FFT and save Files                
            nfft = 2^nextpow2(EEG.pnts);
            NrFreqs = nfft/2+1;
            FreqData = zeros(length({EEG.chanlocs.labels}), nfft, EEG.trials);
            Data = EEG.data();
            Data = Data.*[hamming(EEG.pnts)]';
            power = NaN(EEG.nbchan, NrFreqs, EEG.trials);
            
            for n = 1:EEG.nbchan
                for itrial = 1:EEG.trials
                    FreqData(n, :, itrial) = fft(Data(n, :, itrial),nfft)/EEG.pnts;
                end
            end

            power(:,1:NrFreqs,1:EEG.trials) = abs(FreqData(:,1:NrFreqs, :)).^2;
            Average = horzcat({EEG.chanlocs.labels}', num2cell(mean(power(:,1:60,:),3)));
    
            if i_run == 1
                writecell(Average, fullfile(DatafolderPreEyesOpen, ['rest_' num2str(i_participant) '_pre_Average_EO.csv'])); 
            else
                writecell(Average, fullfile(DatafolderPostEyesOpen, ['rest_' num2str(i_participant) '_post_Average_EO.csv']));
            end

            %% Error Management       
            savetofile(Log, fullfile(LogFolder, ['rest' num2str(i_run) '_' num2str(i_participant) '_Log.mat']));

        catch error
            ErrorMessage = string(error.message);
            for ierrors = 1:length(error.stack)
                ErrorMessage = strcat(ErrorMessage, '//', num2str(error.stack(ierrors).name), ', Line: ', num2str(error.stack(ierrors).line));
            end
            writematrix(ErrorMessage, fullfile(ErrorFolder, ['rest' num2str(i_run) '_' num2str(i_participant) '.csv']))
            savetofile(Log, fullfile(LogFolder, ['rest' num2str(i_run) '_' num2str(i_participant) '_Log.mat']))
        end        
    end % for run
end % for parfor
