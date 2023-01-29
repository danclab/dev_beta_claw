function pipeline_06_extract_bursts(study_info, clusters, clus_name, varargin)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

addpath('/home/bonaiuto/burst_detection/matlab');

pipeline='NEARICA';

subj_id=study_info.participant_info.participant_id{1};
subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);    
EEG=pop_loadset('filepath', subject_data_dir, 'filename', fname); 

load(fullfile(study_info.data_dir, 'derivatives', pipeline, sprintf('processed_%s_fois.mat',clus_name)));
ranges=cell2mat({fois.range}');
[~,sorted_idx]=sort(ranges(:,1));
ranges=ranges(sorted_idx,:);

if strcmp(clus_name,'C')
    beta_idx=max(find(ranges(:,2)<30));
    band_lims=ranges(beta_idx,:);
elseif strcmp(clus_name,'P')
    beta_idx=max(find(ranges(:,2)<30))+1;
    band_lims=ranges(beta_idx,:);
end

[lags,foi,lagged_coh]=load_lagged_coherence(study_info, pipeline);        
lag_part=find(lags<=4.5);

clus_chan_idx=[];
for c=1:length(clusters)
    cluster=clusters{c};
    cluster_idx=find(strcmp(study_info.clusters,cluster));
    chan_idx=cellfun(@(x) find(strcmp({EEG.chanlocs.labels},x)),...
        study_info.cluster_channels{cluster_idx});
    clus_chan_idx(end+1:end+length(chan_idx))=chan_idx;
end

band_idx=knnsearch(foi',band_lims');

band_lc=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,clus_chan_idx,band_idx(1):band_idx(2),lag_part),3),2),1));
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
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
    tf_fname='processed_superlet_tf.mat';
    if exist(fullfile(subject_data_dir,fname),'file')==2 && exist(fullfile(subject_data_dir,tf_fname),'file')==2
        
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

        exe_base_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'EBM')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        exe_exp_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'FTGE')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        
        obs_base_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'OBM')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        obs_exp_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'FTGO')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        if length(exe_base_trials)>=5 && length(exe_exp_trials)>=5 && length(obs_base_trials)>=5 && length(obs_exp_trials)>=5
            conditions={'exe','exe','obs','obs'};
            epochs={'base','exp','base','exp'};
            cond_evts={'EBM','FTGE','OBM','FTGO'};
                        
            load(fullfile(subject_data_dir,tf_fname));

            for clus_idx=1:length(clusters)
                cluster=clusters{clus_idx};
                cluster_idx=find(strcmp(study_info.clusters, cluster));

                cluster_chans=study_info.cluster_channels{cluster_idx};
                for c_idx=1:length(cluster_chans)
                    chan=cluster_chans{c_idx};
                    chan_idx=find(strcmp({EEG.chanlocs.labels},chan));

                    %% Compute FOOOF threshold
                    base_trials=find(cellfun(@length,[cellfun(@(x) find((strcmp(x,'EBM')) | (strcmp(x,'OBM'))), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
                    base_tf=squeeze(trial_tf(base_trials,:,chan_idx,:));
                    % Compute PSD by averaging over time, and then trials
                    base_psd=mean(mean(base_tf,3),1);
                    f_idx=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                    % Rough 1/f fit
                    oof=log10(1./foi(f_idx));
                    lm_psd=fitlm(oof,log10(base_psd(f_idx)),'RobustOpts','on');
                    % FOOOF threshold (need to transform back from log10)
                    fooof_thresh=10.^predict(lm_psd,log10(1./foi)');
                    
                    for cond_idx=1:length(conditions)
                        condition=conditions{cond_idx};
                        epoch=epochs{cond_idx};
                        cond_evt=cond_evts{cond_idx};
                        trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,cond_evt)), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
                
                        chan_tf=squeeze(trial_tf(trials,:,chan_idx,:));
            
                        % Raw data for each trial in this condition
                        raw_trials=squeeze(EEG.data(chan_idx,:,trials))';

                        % Range to look for beta bursts in
                        search_range=find((foi>=search_lims(1)) & (foi<=search_lims(2)));
                        search_freq=foi(search_range);
                        % Sampling rate
                        sfreq=EEG.srate;

                        % Extract bursts
                        ch_bursts=extract_bursts(raw_trials, chan_tf(:,search_range,:),...
                            time, search_freq, band_lims, fooof_thresh(search_range),...
                            sfreq, 'win_size', win_size);
                        n_bursts=length(ch_bursts.trial);
                        ch_bursts.cluster=repmat({cluster},1,n_bursts);
                        ch_bursts.chan=ones(1,n_bursts).*chan_idx;
                        ch_bursts.epoch=repmat({epoch},1,n_bursts);
                        ch_bursts.condition=repmat({condition},1,n_bursts);

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
                        bursts.condition(end+1:end+n_bursts)=ch_bursts.condition;
                        if ~isempty(ch_bursts.waveform_times)
                            bursts.waveform_times=ch_bursts.waveform_times;
                        end
                        
                    end
                end
            end
            save(fullfile(study_info.data_dir,'derivatives', pipeline, subj_id, 'processed_data',sprintf('processed_%s_bursts.mat',clus_name)),'bursts','-v7.3');
        end

    end
end
disp('done');
