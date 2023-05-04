function analysis_05_extract_bursts(study_info, varargin)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

addpath('/home/bonaiuto/burst_detection/matlab');


pipeline='NEARICA_behav';
exe_epoch_types={'EBM','LEXT','FTGE','EXGC','EXEND'};
obs_epoch_types={'OBM','LOBS','FTGO','OBGC','OBEND'};

% behavior
exe_df=readtable(fullfile(study_info.data_dir,'behavior_exe.csv'));

subj_id=study_info.participant_info.participant_id{1};
subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);    
EEG=pop_loadset('filepath', subject_data_dir, 'filename', fname); 

load(fullfile('/home/bonaiuto/dev_beta_umd/data',study_info.age,'derivatives',pipeline,'processed_psd.mat'));
beta_idx=max(find(foi_ranges(:,2)<=30));
band_lims=foi_ranges(beta_idx,:);

[lags,foi,lagged_coh]=load_lagged_coherence(study_info, pipeline);        
lag_part=find(lags<=4.5);

c3_idx=find(strcmp(study_info.clusters,'C3'));
c3_chan_idx=cellfun(@(x) find(strcmp({EEG.chanlocs.labels},x)),...
        study_info.cluster_channels{c3_idx});
c4_idx=find(strcmp(study_info.clusters,'C4'));
c4_chan_idx=cellfun(@(x) find(strcmp({EEG.chanlocs.labels},x)),...
        study_info.cluster_channels{c4_idx});

band_idx=knnsearch(foi',band_lims');

band_lc=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,[c3_chan_idx c4_chan_idx],band_idx(1):band_idx(2),lag_part),3),2),1));
band_lc=band_lc-min(band_lc);
band_lc=band_lc./max(band_lc);
 
[pk_val,pk_idx]=max(band_lc);
fwhm_idx=min(find(band_lc<=.5*pk_val));
fwhm_lag=lags(fwhm_idx);
 
band_center=mean(band_lims);
cycle_dur=1/band_center;
win_size=2*cycle_dur*fwhm_lag;

search_lims=[band_lims(1)-5 band_lims(2)+5];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    disp(subj_id);
    
    % Subject rows in behav data
    exe_subj_rows=find(strcmp(subj_id,exe_df.Subject));
        
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
    c3_tf_fname='processed_superlet_tf_C3.mat';
    c4_tf_fname='processed_superlet_tf_C4.mat';
    if exist(fullfile(subject_data_dir,fname),'file')==2 && exist(fullfile(subject_data_dir,c3_tf_fname),'file')==2 && exist(fullfile(subject_data_dir,c4_tf_fname),'file')==2
        
        % Load data
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);        
        time=EEG.times;
        
        bursts=[];
        bursts.trial=[];
        bursts.waveform=[];
        bursts.peak_freq=[];
        bursts.peak_amp_iter=[];
        bursts.peak_amp_base=[];
        bursts.peak_time=[];
        bursts.peak_adjustment=[];
        bursts.fwhm_freq=[];
        bursts.fwhm_time=[];
        bursts.polarity=[];
        bursts.waveform_times=[];
        bursts.cluster={};
        bursts.chan=[];
        bursts.epoch={};
        bursts.condition={};

        load(fullfile(subject_data_dir,c3_tf_fname));
        c3_trial_tf=trial_tf;
        load(fullfile(subject_data_dir,c4_tf_fname));
        c4_trial_tf=trial_tf;

        %% Compute FOOOF threshold
        f_idx=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
        base_trials=find(cellfun(@length,[cellfun(@(x) find((strcmp(x,'EBM')) | (strcmp(x,'OBM'))), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        
        c3_base_tf=squeeze(c3_trial_tf(base_trials,:,c3_chan_idx,:));
        % Compute PSD by averaging over time, and then trials
        c3_base_psd=mean(mean(c3_base_tf,3),1);
        % Rough 1/f fit
        c3_oof=log10(1./foi(f_idx));
        c3_lm_psd=fitlm(c3_oof,log10(c3_base_psd(f_idx)),'RobustOpts','on');
        % FOOOF threshold (need to transform back from log10)
        c3_fooof_thresh=10.^predict(c3_lm_psd,log10(1./foi)');
        
        c4_base_tf=squeeze(c4_trial_tf(base_trials,:,c4_chan_idx,:));
        % Compute PSD by averaging over time, and then trials
        c4_base_psd=mean(mean(c4_base_tf,3),1);
        % Rough 1/f fit
        c4_oof=log10(1./foi(f_idx));
        c4_lm_psd=fitlm(c4_oof,log10(c4_base_psd(f_idx)),'RobustOpts','on');
        % FOOOF threshold (need to transform back from log10)
        c4_fooof_thresh=10.^predict(c4_lm_psd,log10(1./foi)');

        for i=1:length(exe_epoch_types)
        
            epoch_type=exe_epoch_types{i};            
            epoch_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,epoch_type)), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
            
            if length(epoch_trials)>=5
                
                for j=1:length(c3_chan_idx)
                    chan_tf=squeeze(c3_trial_tf(epoch_trials,j,:,:));
            
                    % Raw data for each trial in this condition
                    raw_trials=squeeze(EEG.data(c3_chan_idx(j),:,epoch_trials))';

                    % Range to look for beta bursts in
                    search_range=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                    search_freq=foi(search_range);
                    % Sampling rate
                    sfreq=EEG.srate;

                    % Extract bursts
                    ch_bursts=extract_bursts(raw_trials, chan_tf(:,search_range,:),...
                        time, search_freq, band_lims, c3_fooof_thresh(search_range),...
                        sfreq, 'win_size', win_size);
                    n_bursts=length(ch_bursts.trial);
                    ch_bursts.cluster=repmat({''},1,n_bursts);
                    for t=1:length(epoch_trials)
                        trial=EEG.epoch(epoch_trials(t)).eventtrial;
                        trial_row=find(exe_df.Trial(exe_subj_rows)==trial);
                        trial_bursts=find(ch_bursts.trial==trial);
                        if strcmp(exe_df.FTGhand(exe_subj_rows(trial_row)),'L')                            
                            ch_bursts.cluster(trial_bursts)=repmat({'ipsi'},1,length(trial_bursts));
                        elseif strcmp(exe_df.FTGhand(exe_subj_rows(trial_row)),'R')                            
                            ch_bursts.cluster(trial_bursts)=repmat({'contra'},1,length(trial_bursts));
                        elseif ~strcmp(epoch_type,'EBM')
                            disp('error');
                        end
                    end
                    ch_bursts.chan=ones(1,n_bursts).*c3_chan_idx(j);
                    ch_bursts.epoch=repmat({epoch_type},1,n_bursts);
                    
                    % Add to list of all bursts
                    bursts.trial(end+1:end+n_bursts)=ch_bursts.trial;
                    bursts.waveform(end+1:end+n_bursts,:)=ch_bursts.waveform;
                    bursts.peak_freq(end+1:end+n_bursts)=ch_bursts.peak_freq;
                    bursts.peak_amp_iter(end+1:end+n_bursts)=ch_bursts.peak_amp_iter;
                    bursts.peak_amp_base(end+1:end+n_bursts)=ch_bursts.peak_amp_base;
                    bursts.peak_time(end+1:end+n_bursts)=ch_bursts.peak_time;
                    bursts.peak_adjustment(end+1:end+n_bursts)=ch_bursts.peak_adjustment;
                    bursts.fwhm_freq(end+1:end+n_bursts)=ch_bursts.fwhm_freq;
                    bursts.fwhm_time(end+1:end+n_bursts)=ch_bursts.fwhm_time;
                    bursts.polarity(end+1:end+n_bursts)=ch_bursts.polarity;
                    bursts.cluster(end+1:end+n_bursts)=ch_bursts.cluster;
                    bursts.chan(end+1:end+n_bursts)=ch_bursts.chan;
                    bursts.epoch(end+1:end+n_bursts)=ch_bursts.epoch;
                    if ~isempty(ch_bursts.waveform_times)
                        bursts.waveform_times=ch_bursts.waveform_times;
                    end
                end

                for j=1:length(c4_chan_idx)
                    chan_tf=squeeze(c4_trial_tf(epoch_trials,j,:,:));
            
                    % Raw data for each trial in this condition
                    raw_trials=squeeze(EEG.data(c4_chan_idx(j),:,epoch_trials))';

                    % Range to look for beta bursts in
                    search_range=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                    search_freq=foi(search_range);
                    % Sampling rate
                    sfreq=EEG.srate;

                    % Extract bursts
                    ch_bursts=extract_bursts(raw_trials, chan_tf(:,search_range,:),...
                        time, search_freq, band_lims, c4_fooof_thresh(search_range),...
                        sfreq, 'win_size', win_size);
                    n_bursts=length(ch_bursts.trial);
                    ch_bursts.cluster=repmat({''},1,n_bursts);
                    for t=1:length(epoch_trials)
                        trial=EEG.epoch(epoch_trials(t)).eventtrial;
                        trial_row=find(exe_df.Trial(exe_subj_rows)==trial);
                        trial_bursts=find(ch_bursts.trial==trial);
                        if strcmp(exe_df.FTGhand(exe_subj_rows(trial_row)),'L')                            
                            ch_bursts.cluster(trial_bursts)=repmat({'contra'},1,length(trial_bursts));
                        elseif strcmp(exe_df.FTGhand(exe_subj_rows(trial_row)),'R')                            
                            ch_bursts.cluster(trial_bursts)=repmat({'ipsi'},1,length(trial_bursts));
                        elseif ~strcmp(epoch_type,'EBM')
                            disp('error');
                        end
                    end
                    ch_bursts.chan=ones(1,n_bursts).*c4_chan_idx(j);
                    ch_bursts.epoch=repmat({epoch_type},1,n_bursts);
                    
                    % Add to list of all bursts
                    bursts.trial(end+1:end+n_bursts)=ch_bursts.trial;
                    bursts.waveform(end+1:end+n_bursts,:)=ch_bursts.waveform;
                    bursts.peak_freq(end+1:end+n_bursts)=ch_bursts.peak_freq;
                    bursts.peak_amp_iter(end+1:end+n_bursts)=ch_bursts.peak_amp_iter;
                    bursts.peak_amp_base(end+1:end+n_bursts)=ch_bursts.peak_amp_base;
                    bursts.peak_time(end+1:end+n_bursts)=ch_bursts.peak_time;
                    bursts.peak_adjustment(end+1:end+n_bursts)=ch_bursts.peak_adjustment;
                    bursts.fwhm_freq(end+1:end+n_bursts)=ch_bursts.fwhm_freq;
                    bursts.fwhm_time(end+1:end+n_bursts)=ch_bursts.fwhm_time;
                    bursts.polarity(end+1:end+n_bursts)=ch_bursts.polarity;
                    bursts.cluster(end+1:end+n_bursts)=ch_bursts.cluster;
                    bursts.chan(end+1:end+n_bursts)=ch_bursts.chan;
                    bursts.epoch(end+1:end+n_bursts)=ch_bursts.epoch;
                    if ~isempty(ch_bursts.waveform_times)
                        bursts.waveform_times=ch_bursts.waveform_times;
                    end
                end
            end
        end
        
        for i=1:length(obs_epoch_types)
        
            epoch_type=obs_epoch_types{i};            
            epoch_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,epoch_type)), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
            
            if length(epoch_trials)>=5
                
                for j=1:length(c3_chan_idx)
                    chan_tf=squeeze(c3_trial_tf(epoch_trials,j,:,:));
            
                    % Raw data for each trial in this condition
                    raw_trials=squeeze(EEG.data(c3_chan_idx(j),:,epoch_trials))';

                    % Range to look for beta bursts in
                    search_range=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                    search_freq=foi(search_range);
                    % Sampling rate
                    sfreq=EEG.srate;

                    % Extract bursts
                    ch_bursts=extract_bursts(raw_trials, chan_tf(:,search_range,:),...
                        time, search_freq, band_lims, c3_fooof_thresh(search_range),...
                        sfreq, 'win_size', win_size);
                    n_bursts=length(ch_bursts.trial);
                    ch_bursts.cluster=repmat({'C3'},1,n_bursts);
                    ch_bursts.chan=ones(1,n_bursts).*c3_chan_idx(j);
                    ch_bursts.epoch=repmat({epoch_type},1,n_bursts);
                    
                    % Add to list of all bursts
                    bursts.trial(end+1:end+n_bursts)=ch_bursts.trial;
                    bursts.waveform(end+1:end+n_bursts,:)=ch_bursts.waveform;
                    bursts.peak_freq(end+1:end+n_bursts)=ch_bursts.peak_freq;
                    bursts.peak_amp_iter(end+1:end+n_bursts)=ch_bursts.peak_amp_iter;
                    bursts.peak_amp_base(end+1:end+n_bursts)=ch_bursts.peak_amp_base;
                    bursts.peak_time(end+1:end+n_bursts)=ch_bursts.peak_time;
                    bursts.peak_adjustment(end+1:end+n_bursts)=ch_bursts.peak_adjustment;
                    bursts.fwhm_freq(end+1:end+n_bursts)=ch_bursts.fwhm_freq;
                    bursts.fwhm_time(end+1:end+n_bursts)=ch_bursts.fwhm_time;
                    bursts.polarity(end+1:end+n_bursts)=ch_bursts.polarity;
                    bursts.cluster(end+1:end+n_bursts)=ch_bursts.cluster;
                    bursts.chan(end+1:end+n_bursts)=ch_bursts.chan;
                    bursts.epoch(end+1:end+n_bursts)=ch_bursts.epoch;
                    if ~isempty(ch_bursts.waveform_times)
                        bursts.waveform_times=ch_bursts.waveform_times;
                    end
                end

                for j=1:length(c4_chan_idx)
                    chan_tf=squeeze(c4_trial_tf(epoch_trials,j,:,:));
            
                    % Raw data for each trial in this condition
                    raw_trials=squeeze(EEG.data(c4_chan_idx(j),:,epoch_trials))';

                    % Range to look for beta bursts in
                    search_range=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                    search_freq=foi(search_range);
                    % Sampling rate
                    sfreq=EEG.srate;

                    % Extract bursts
                    ch_bursts=extract_bursts(raw_trials, chan_tf(:,search_range,:),...
                        time, search_freq, band_lims, c4_fooof_thresh(search_range),...
                        sfreq, 'win_size', win_size);
                    n_bursts=length(ch_bursts.trial);
                    ch_bursts.cluster=repmat({'C4'},1,n_bursts);
                    ch_bursts.chan=ones(1,n_bursts).*c4_chan_idx(j);
                    ch_bursts.epoch=repmat({epoch_type},1,n_bursts);
                    
                    % Add to list of all bursts
                    bursts.trial(end+1:end+n_bursts)=ch_bursts.trial;
                    bursts.waveform(end+1:end+n_bursts,:)=ch_bursts.waveform;
                    bursts.peak_freq(end+1:end+n_bursts)=ch_bursts.peak_freq;
                    bursts.peak_amp_iter(end+1:end+n_bursts)=ch_bursts.peak_amp_iter;
                    bursts.peak_amp_base(end+1:end+n_bursts)=ch_bursts.peak_amp_base;
                    bursts.peak_time(end+1:end+n_bursts)=ch_bursts.peak_time;
                    bursts.peak_adjustment(end+1:end+n_bursts)=ch_bursts.peak_adjustment;
                    bursts.fwhm_freq(end+1:end+n_bursts)=ch_bursts.fwhm_freq;
                    bursts.fwhm_time(end+1:end+n_bursts)=ch_bursts.fwhm_time;
                    bursts.polarity(end+1:end+n_bursts)=ch_bursts.polarity;
                    bursts.cluster(end+1:end+n_bursts)=ch_bursts.cluster;
                    bursts.chan(end+1:end+n_bursts)=ch_bursts.chan;
                    bursts.epoch(end+1:end+n_bursts)=ch_bursts.epoch;
                    if ~isempty(ch_bursts.waveform_times)
                        bursts.waveform_times=ch_bursts.waveform_times;
                    end
                end
            end
        end
        save(fullfile(study_info.data_dir,'derivatives', pipeline, subj_id, 'processed_data','processed_bursts.mat'),'bursts','-v7.3');
    end
end
disp('done');
