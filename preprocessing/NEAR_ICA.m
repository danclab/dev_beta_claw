function NEAR_ICA(study_info)

addpath('NoiseTools');
addpath('/home/bonaiuto/eeglab2021.1');

ext='.set';

% Channel location file
channel_locations = '/home/bonaiuto/dev_beta_umd/data/GSN-HydroCel-64.sfp';

% Initialize the filters
% High-pass frequency
highpass = 1;
% Low-pass frequency
lowpass  = 100;

% event/condition markers
event_markers = {'OBM','EBM','FTGO','FTGE'};

% epoch length in seconds
epoch_length = [-1.5 1.5];

% lower and upper voltage threshold (in mV)
volt_threshold = [-150 150];

% list of frontal channels to check for epoch-level channel interpolation
% (see manuscript for detail)
frontal_channels = {'E1', 'E5', 'E10', 'E17'}; 

% Parameters for NEAR - Bad Channels Detection
% flat channels
isFlat  = 1; % flag variable to enable or disable Flat-lines detection method (default: 1)
flatWin = 5; % tolerance level in s(default: 5)

% LOF (density-based)
isLOF       = 1;  % flag variable to enable or disable LOF method (default: 1)
dist_metric = 'seuclidean'; % Distance metric to compute k-distance
thresh_lof  = 2.5; % Threshold cut-off for outlier detection on LOF scores
isAdapt = 10; % The threshold will be incremented by a factor of 1 if the given threshold detects more than xx %
%of total channels (eg., 10); if this variable left empty [], no adaptive thresholding is enabled.

% Periodogram (frequency based)
isPeriodogram = 0; % flag variable to enable or disable periodogram method (default: 0)
frange        = [1 20]; % Frequency Range in Hz
winsize       = 1; % window length in s
winov         = 0.66; % 66% overlap factor
pthresh       = 4.5; % Threshold Factor to predict outliers on the computed energy

% Parameters for NEAR- Bad Segments Correction/Rejection using ASR %
rej_cutoff = 13;   % A lower value implies severe removal (Recommended value range: 20 to 30)
rej_mode   = 'off'; % Set to 'off' for ASR Correction and 'on for ASR Removal (default: 'on')
add_reject = 'off'; % Set to 'on' for additional rejection of bad segments if any after ASR processing (default: 'off')

% Parameter for interpolation %
interp_type = 'spherical'; % other values can be 'v4'. Please refer to pop_interp.m for more details.

%% Initialize output variables
lof_flat_channels={};
lof_channels={};
lof_periodo_channels={};
% Bad channels identified using LOF
lof_bad_channels={};
% number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
ica_preparation_bad_channels=[];
% length of data (in second) fed into ICA decomposition
length_ica_data=[];
% total independent components (ICs)
total_ICs=[];
% number of artifacted ICs
ICs_removed=[];
% number of epochs before artifact rejection
total_epochs_before_artifact_rejection=[];
% number of epochs after artifact rejection
total_epochs_after_artifact_rejection=[];
% total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
total_channels_interpolated=[];
asr_tot_samples_modified=[];
asr_change_in_RMS=[];

% behavior
exe_df=readtable(fullfile(study_info.data_dir,'behavior_exe.csv'));

%% Loop over all data files
for s_idx=1:size(study_info.participant_info,1)
    % Get subject ID from study info
    subject=study_info.participant_info.participant_id{s_idx};
    
    % Where original raw data is located
    subject_raw_data_dir=fullfile(study_info.data_dir, 'data',subject, 'eeg');
    
    % Where to put processed (derived) data
    subject_output_data_dir=fullfile(study_info.data_dir, 'derivatives', 'NEARICA_behav', subject);
    
    if exist([subject_output_data_dir filesep '01_filtered_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '01_filtered_data'])
    end

    if exist([subject_output_data_dir filesep '02_near_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '02_near_data'])
    end

    if exist([subject_output_data_dir filesep '03_ica_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '03_ica_data'])
    end

    if exist([subject_output_data_dir filesep '04_rereferenced_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '04_rereferenced_data'])
    end
    
    if exist([subject_output_data_dir filesep 'processed_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep 'processed_data'])
    end
    
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subject);
    
    %% Import data
    data_file_name=sprintf('%s_task-%s_eeg.set',subject, study_info.task);
    EEG=pop_loadset('filename', data_file_name, 'filepath', subject_raw_data_dir);
    EEG = eeg_checkset(EEG);
    origEEG=EEG;
    
    % Remove events after last LBBS event
    lbbs_evts=find(strcmp('LBBS',{EEG.event.type}));
    EEG.event=EEG.event(1:lbbs_evts(end));
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    %% Label the task
    % Add observation baseline events
    base_obs_latencies=[EEG.event(find(strcmp('LBOB',{EEG.event.type}))).latency];
    for i=1:length(base_obs_latencies)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='OBM';
        EEG.event(n_events+1).latency=(base_obs_latencies(i)-3*EEG.srate)-1;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Add execution baseline events
    exe_obs_latencies=[EEG.event(find(strcmp('LBEX',{EEG.event.type}))).latency];
    for i=1:length(exe_obs_latencies)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='EBM';
        EEG.event(n_events+1).latency=(exe_obs_latencies(i)-3*EEG.srate)-1;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    
    % Add grasp completion and trial end events
    start_evts=find(strcmp('LEXT',{EEG.event.type}));
    end_evts=find(strcmp('LBBS',{EEG.event.type}));
    end_times=[EEG.event(end_evts).latency];
    
    exe_subj_rows=find(strcmp(subject,exe_df.Subject));
    for i=1:length(start_evts)
        start_time=EEG.event(start_evts(i)).latency;
        end_evt_time=end_times(find(end_times>start_time,1));
        
        if exe_df.GC(exe_subj_rows(i))>0
            gc_time=exe_df.Gcfinal(exe_subj_rows(i));
            n_events=length(EEG.event);
            EEG.event(n_events+1).type='EXGC';
            EEG.event(n_events+1).latency=(gc_time/1000)*EEG.srate;
            EEG.event(n_events+1).urevent=n_events+1;
        end
        
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='EXEND';
        EEG.event(n_events+1).latency=end_evt_time;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    
    %% Delete discontinuous data from the raw data file
    % remove data after last task event
    end_flags=find(strcmp({EEG.event.type},'LBBS'));
    latency=EEG.event(end_flags(end)).latency;
    % remove everything 1.5 seconds after the last event
    EEG = eeg_eegrej( EEG, [(latency+(1.5*EEG.srate)) EEG.pnts] );
    EEG = eeg_checkset( EEG );
    
    lbse_flags = find(strcmp({EEG.event.type},'LBSE'));
    latency=EEG.event(lbse_flags(1)).latency;
    % remove everything until 1.5 seconds before the first event
    EEG = eeg_eegrej( EEG, [1 (latency-(1.5*EEG.srate))] );
    EEG = eeg_checkset( EEG );
    
    %% Import channel locations
    EEG=pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    
    % Check whether the channel locations were properly imported. The EEG
    % signals and channel numbers should be same.
    if size(EEG.data, 1) ~= length(EEG.chanlocs)
        error('The size of the data does not match with channel numbers.');
    end
    
    % Plot channel layout
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'01-initial_ch_locations.png'));
        
    %% Adjust anti-aliasing and task related time offset
    % adjust anti-aliasing filter time offset
    filter_timeoffset = 18;
    for aafto=1:length(EEG.event)
        if ~strcmp(EEG.event(aafto).type,'FTGO') && ~strcmp(EEG.event(aafto).type,'FTGE') && ~strcmp(EEG.event(aafto).type,'EXGC') && ~strcmp(EEG.event(aafto).type,'EXEND')
            EEG.event(aafto).latency=EEG.event(aafto).latency+(filter_timeoffset/1000)*EEG.srate;
        end
    end
    
    %% Filter data
    % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
    % df = transition band width, dF = normalized transition width, fs = sampling rate
    % dF is specific for the window type. Hamming window dF = 3.3
    
    high_transband = highpass; % high pass transition band
    low_transband = 10; % low pass transition band
    
    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband / EEG.srate);
    
    % Round filter order to next higher even integer. Filter order is always even integer.
    if mod(floor(hp_fl_order),2) == 0
        hp_fl_order=floor(hp_fl_order);
    elseif mod(floor(hp_fl_order),2) == 1
        hp_fl_order=floor(hp_fl_order)+1;
    end
    
    if mod(floor(lp_fl_order),2) == 0
        lp_fl_order=floor(lp_fl_order)+2;
    elseif mod(floor(lp_fl_order),2) == 1
        lp_fl_order=floor(lp_fl_order)+1;
    end
    
    % Calculate cutoff frequency
    high_cutoff = highpass/2;
    low_cutoff = lowpass + (low_transband/2);
    
    % Performing high pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass',...
        'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % Performing low pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass',...
        'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % Plot PSD
    fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'05-filtered_psd.png'));
    
    %% NEAR Bad Channel Detection        
    [EEG, flat_ch, lof_ch, periodo_ch, LOF_vec, thresh_lof_update] = NEAR_getBadChannels(EEG, isFlat, flatWin, isLOF, thresh_lof, dist_metric, isAdapt, ...
        isPeriodogram, frange, winsize, winov, pthresh, 0);
    save(fullfile(subject_output_data_dir, 'LOF_Values.mat'), 'LOF_vec'); % save .mat format
    disp('Bad Channel Detection is performed successfully');
    badChans = sort(unique(union(union(flat_ch, lof_ch),periodo_ch)));

    if(~isempty(badChans))
        if(size(badChans,1) ~= 1)
            badChans = badChans';
        end
    end

    EEG = pop_select(EEG, 'nochannel', badChans);

    lof_flat_channels{s_idx}='';
    if numel(flat_ch)>0
        lof_flat_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(flat_ch,3), 'UniformOutput', false)',',');
    end
    lof_channels{s_idx}='';
    if numel(lof_ch)>0
        lof_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(lof_ch,3), 'UniformOutput', false)',',');
    end
    lof_periodo_channels{s_idx}='';
    if numel(periodo_ch)>0
        lof_periodo_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(periodo_ch,3), 'UniformOutput', false)',',');
    end
    lof_bad_channels{s_idx}='';
    if numel(badChans)>0
        lof_bad_channels(s_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(badChans,3), 'UniformOutput', false)',',');
    end
    
    %% Save data after running filter and LOF function
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', strrep(data_file_name, ext, '_filtered_data'));
    EEG = pop_saveset( EEG,'filename',strrep(data_file_name, ext, '_filtered_data.set'),...
        'filepath', [subject_output_data_dir filesep '01_filtered_data' filesep]); % save .set format
    
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'06-lof_removed.png'));
    
    %% Bad epochs correction/removal using ASR
    EEG_copy = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off', ...
        'Highpass','off','BurstCriterion',rej_cutoff,'WindowCriterion',add_reject,'BurstRejection',rej_mode,'Distance','Euclidian');

    if(strcmp(rej_mode, 'on'))
        modified_mask = ~EEG.etc.clean_sample_mask;
    else
        modified_mask = sum(abs(EEG_copy.data-EEG.data),1) > 1e-10;
    end

    tot_samples_modified = (length(find(modified_mask)) * 100) / EEG_copy.pnts;
    tot_samples_modified = round(tot_samples_modified * 100) / 100;
    asr_tot_samples_modified(s_idx)=tot_samples_modified;
    change_in_RMS = -(mean(rms(EEG.data,2)) - mean(rms(EEG_copy.data,2))*100)/mean(rms(EEG_copy.data,2)); % in percentage
    change_in_RMS = round(change_in_RMS * 100) / 100;
    asr_change_in_RMS(s_idx) =change_in_RMS;
    fprintf('\nArtifacted epochs are corrected by ASR algorithm\n');
    
    %% Save data after running ASR function
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', strrep(data_file_name, ext, '_asr_data'));
    EEG = pop_saveset( EEG,'filename',strrep(data_file_name, ext, '_asr_data.set'),...
        'filepath', [subject_output_data_dir filesep '02_near_data' filesep]); % save .set format

    fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'07-asr_psd.png'));
    
    %% STEP 8: Prepare data for ICA
    EEG_copy=EEG; % make a copy of the dataset
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Perform 1Hz high pass filter on copied dataset
    transband = 1;
    fl_cutoff = transband/2;
    fl_order = 3.3 / (transband / EEG.srate);
    
    if mod(floor(fl_order),2) == 0
        fl_order=floor(fl_order);
    elseif mod(floor(fl_order),2) == 1
        fl_order=floor(fl_order)+1;
    end
    
    EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff,...
        'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order,...
        'minphase', 0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Create 1 second epoch
    % insert temporary marker 1 second apart and create epochs
    EEG_copy=eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1],...
        'rmbase', [NaN], 'eventtype', '999'); 
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Find bad epochs and delete them from dataset
    % [lower upper] threshold limit(s) in mV.
    vol_thrs = [-1000 1000]; 
    
    % Find channel/s with xx% of artifacted 1-second epochs and delete them
    chanCounter = 1; ica_prep_badChans = [];
    numEpochs =EEG_copy.trials; % find the number of epochs
    all_bad_channels=0;
    
    for ch=1:EEG_copy.nbchan
        % Find artifaceted epochs by detecting outlier voltage
        EEG_copy = pop_eegthresh(EEG_copy,1, ch, vol_thrs(1), vol_thrs(2),...
            EEG_copy.xmin, EEG_copy.xmax, 0, 0);
        EEG_copy = eeg_checkset( EEG_copy );
        
        % Find number of artifacted epochs
        EEG_copy = eeg_checkset( EEG_copy );
        EEG_copy = eeg_rejsuperpose( EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        artifacted_epochs=EEG_copy.reject.rejglobal;
        
        % Find bad channel / channel with more than 20% artifacted epochs
        if sum(artifacted_epochs) > (numEpochs*20/100)
            ica_prep_badChans(chanCounter) = ch;
            chanCounter=chanCounter+1;
        end
    end
    
    % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
    if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan
        all_bad_channels=1;
        warning(['No usable data for datafile', data_file_name]);        
    else
        % Reject bad channel - channel with more than xx% artifacted epochs
        EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
        EEG_copy = eeg_checkset(EEG_copy);
    end
    
    if numel(ica_prep_badChans)==0
        ica_preparation_bad_channels{s_idx}='0';
    else
        ica_preparation_bad_channels{s_idx}=num2str(ica_prep_badChans);
    end
    
    if all_bad_channels == 1
        length_ica_data(s_idx)=0;
        total_ICs(s_idx)=0;
        ICs_removed{s_idx}='0';
        total_epochs_before_artifact_rejection(s_idx)=0;
        total_epochs_after_artifact_rejection(s_idx)=0;
        total_channels_interpolated(s_idx)=0;
        continue % ignore rest of the processing and go to next datafile
    end
    
    % Find the artifacted epochs across all channels and reject them before doing ICA.
    EEG_copy = pop_eegthresh(EEG_copy,1, 1:EEG_copy.nbchan, vol_thrs(1),...
        vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax,0,0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Find the number of artifacted epochs and reject them
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
    reject_artifacted_epochs=EEG_copy.reject.rejglobal;
    EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);
    
    fig=compute_and_plot_psd(EEG_copy, 1:EEG_copy.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'08-ica_copy_epochs_psd.png'));    
    
    %% Run ICA
    length_ica_data(s_idx)=EEG_copy.trials; % length of data (in second) fed into ICA
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1,...
        'stop', 1E-7, 'interupt','off');
    
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_editset(EEG_copy, 'setname',  strrep(data_file_name, ext, '_ica'));
    EEG_copy = pop_saveset(EEG_copy, 'filename', strrep(data_file_name, ext, '_ica.set'),...
        'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format
    
    % Find the ICA weights that would be transferred to the original dataset
    ICA_WINV=EEG_copy.icawinv;
    ICA_SPHERE=EEG_copy.icasphere;
    ICA_WEIGHTS=EEG_copy.icaweights;
    ICA_CHANSIND=EEG_copy.icachansind;
    
    % If channels were removed from copied dataset during preparation of ica, then remove
    % those channels from original dataset as well before transferring ica weights.
    EEG = eeg_checkset(EEG);
    EEG = pop_select(EEG,'nochannel', ica_prep_badChans);
    
    % Transfer the ICA weights of the copied dataset to the original dataset
    EEG.icawinv=ICA_WINV;
    EEG.icasphere=ICA_SPHERE;
    EEG.icaweights=ICA_WEIGHTS;
    EEG.icachansind=ICA_CHANSIND;
    EEG = eeg_checkset(EEG);
    
    %% Run adjusted-adjust to find artifacted ICA components
    badICs=[];    
    if size(EEG_copy.icaweights,1) == size(EEG_copy.icaweights,2)
        figure()
        badICs = adjusted_ADJUST(EEG_copy, [[subject_output_data_dir filesep '03_ica_data' filesep] strrep(data_file_name, ext, '_adjust_report')]);        
        close all;
    else % if rank is less than the number of electrodes, throw a warning message
        warning('The rank is less than the number of electrodes. ADJUST will be skipped. Artefacted ICs will have to be manually rejected for this participant');
    end
    
    % Mark the bad ICs found by ADJUST
    for ic=1:length(badICs)
        EEG.reject.gcompreject(1, badICs(ic))=1;
        EEG = eeg_checkset(EEG);
    end
    total_ICs(s_idx)=size(EEG.icasphere, 1);
    if numel(badICs)==0
        ICs_removed{s_idx}='0';
    else
        ICs_removed{s_idx}=num2str(double(badICs));
    end    
    
    %% Save dataset after ICA
    EEG = eeg_checkset(EEG);
    EEG = pop_editset(EEG, 'setname',  strrep(data_file_name, ext, '_ica_data'));
    EEG = pop_saveset(EEG, 'filename', strrep(data_file_name, ext, '_ica_data.set'),...
        'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format
    
    %% Remove artifacted ICA components from data
    all_bad_ICs=0;
    ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove
    
    % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
    if numel(ICs2remove)==total_ICs(s_idx)
        all_bad_ICs=1;
        warning(['No usable data for datafile', data_file_name]);        
    else
        EEG = eeg_checkset( EEG );
        EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
    end
    
    if all_bad_ICs==1
        total_epochs_before_artifact_rejection(s_idx)=0;
        total_epochs_after_artifact_rejection(s_idx)=0;
        total_channels_interpolated(s_idx)=0;
        continue % ignore rest of the processing and go to next datafile
    end
    
    fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'09-ica_art_rej_psd.png'));
    
    
    %% Segment data into fixed length epochs
    EEG = eeg_checkset(EEG);
    EEG = pop_epoch(EEG, event_markers, epoch_length, 'epochinfo', 'yes');
    
    fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'10-epoch_psd.png'));
    
    total_epochs_before_artifact_rejection(s_idx)=EEG.trials;
    
    %% Artifact rejection
    all_bad_epochs=0;
    chans=[]; chansidx=[];chans_labels2=[];
    chans_labels2=cell(1,EEG.nbchan);
    for i=1:EEG.nbchan
        chans_labels2{i}= EEG.chanlocs(i).labels;
    end
    [chans,chansidx] = ismember(frontal_channels, chans_labels2);
    frontal_channels_idx = chansidx(chansidx ~= 0);
    badChans = zeros(EEG.nbchan, EEG.trials);
    badepoch=zeros(1, EEG.trials);
    
    % find artifaceted epochs by detecting outlier voltage in the specified channels list and remove epoch if artifacted in those channels
    for ch =1:length(frontal_channels_idx)
        EEG = pop_eegthresh(EEG,1, frontal_channels_idx(ch), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        badChans(ch,:) = EEG.reject.rejglobal;
    end
    for ii=1:size(badChans, 2)
        badepoch(ii)=sum(badChans(:,ii));
    end
    badepoch=logical(badepoch);
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
        all_bad_epochs=1;
        warning(['No usable data for datafile', data_file_name]);                
    else
        EEG = pop_rejepoch( EEG, badepoch, 0);
        EEG = eeg_checkset(EEG);
    end

    if all_bad_epochs==0
        % Interpolate artifacted data for all reaming channels
        badChans = zeros(EEG.nbchan, EEG.trials);
        % Find artifacted epochs by detecting outlier voltage but don't remove
        for ch=1:EEG.nbchan
            EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
            EEG = eeg_checkset(EEG);
            EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            badChans(ch,:) = EEG.reject.rejglobal;
        end
        tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
        for e = 1:EEG.trials
            % Initialize variables EEGe and EEGe_interp;
            EEGe = []; EEGe_interp = []; badChanNum = [];
            % Select only this epoch (e)
            EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
            badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
            if length(badChanNum) < round((10/100)*EEG.nbchan)% check if more than 10% are bad
                EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
            end
        end
        EEG.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

        % If more than 10% of channels in an epoch were interpolated, reject that epoch
        badepoch=zeros(1, EEG.trials);
        for ei=1:EEG.trials
            NumbadChan = badChans(:,ei); % find how many channels are bad in an epoch
            if sum(NumbadChan) > round((10/100)*EEG.nbchan)% check if more than 10% are bad
                badepoch (ei)= sum(NumbadChan);
            end
        end
        badepoch=logical(badepoch);
    end
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
        all_bad_epochs=1;
        warning(['No usable data for datafile', data_file_name]);                
    else
        EEG = pop_rejepoch(EEG, badepoch, 0);
        EEG = eeg_checkset(EEG);
    end
    
    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
    if sum(EEG.reject.rejthresh)==EEG.trials || sum(EEG.reject.rejthresh)+1==EEG.trials
        all_bad_epochs=1;
        warning(['No usable data for datafile', data_file_name]);            
    else
        EEG = pop_rejepoch(EEG,(EEG.reject.rejthresh), 0);
        EEG = eeg_checkset(EEG);
    end
    
    % if all epochs are found bad during artifact rejection
    if all_bad_epochs==1
        total_epochs_after_artifact_rejection(s_idx)=0;
        total_channels_interpolated(s_idx)=0;
        continue % ignore rest of the processing and go to next datafile
    else
        total_epochs_after_artifact_rejection(s_idx)=EEG.trials;        
    end
    
    %% Interpolation  
    total_channels_interpolated(s_idx)=length(origEEG.chanlocs)-length(EEG.chanlocs);
    EEG = pop_interp(EEG, origEEG.chanlocs, interp_type);
    fprintf('\nMissed channels are spherically interpolated\n');
    
    %% Re-referencing
    EEG = pop_reref( EEG, []);
    
    fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'11-art_rej_reref_psd.png'));    
    
    %% Save processed data
    EEG = eeg_checkset(EEG);
    EEG = pop_editset(EEG, 'setname',  strrep(data_file_name, ext, '_rereferenced_data'));
    EEG = pop_saveset(EEG, 'filename', strrep(data_file_name, ext, '_rereferenced_data.set'),...
        'filepath', [subject_output_data_dir filesep '04_rereferenced_data' filesep ]); % save .set format
    
end

%% Create the report table for all the data files with relevant preprocessing outputs.
report_table=table(study_info.participant_info.participant_id,...
    lof_flat_channels', lof_channels', lof_periodo_channels', lof_bad_channels',...
    asr_tot_samples_modified', asr_change_in_RMS', ica_preparation_bad_channels',...
    length_ica_data', total_ICs', ICs_removed', total_epochs_before_artifact_rejection',...
    total_epochs_after_artifact_rejection',total_channels_interpolated');

report_table.Properties.VariableNames={'subject','lof_flat_channels', 'lof_channels', ...
    'lof_periodo_channels', 'lof_bad_channels', 'asr_tot_samples_modified', 'asr_change_in_RMS',...
    'ica_preparation_bad_channels', 'length_ica_data', 'total_ICs', 'ICs_removed', 'total_epochs_before_artifact_rejection', ...
    'total_epochs_after_artifact_rejection', 'total_channels_interpolated'};
writetable(report_table, fullfile(study_info.data_dir, 'derivatives', 'NEARICA_behav', ['NEARICA_preprocessing_report_', datestr(now,'dd-mm-yyyy'),'.csv']));
end