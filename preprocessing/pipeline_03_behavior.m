function pipeline_03_behavior(study_info)

pipeline='NEARICA_behav';

% behavior
exe_df=readtable(fullfile(study_info.data_dir,'behavior_exe.csv'));
obs_df=readtable(fullfile(study_info.data_dir,'behavior_obs.csv'));

% Number of subjects
n_subjects=size(study_info.participant_info,1);

epoch_types={'OBM','LOBS','FTGO','OBGC','OBEND','EBM','LEXT','FTGE','EXGC','EXEND'};
n_epochs=zeros(n_subjects,length(epoch_types));

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id);
    subject_data_dir=fullfile(subject_dir, '04_rereferenced_data');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_rereferenced_data.set',subj_id);
    
    if exist(fullfile(subject_data_dir,fname),'file')==2
        
        % Load data
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);   
        
        load(fullfile(subject_dir,'05_zapped_data','zapline.mat'));
        data=data.*1e6;
        EEG.data=permute(data,[2 3 1]);
        EEG = pop_editset(EEG, 'setname', sprintf('%s_task-tool_obs_exe_eeg_zapped_data',subj_id));
        pop_saveset(EEG, 'filepath', fullfile(subject_dir, '05_zapped_data'),...
            'filename', sprintf('%s_task-tool_obs_exe_eeg_zapped_data.set',subj_id));
        
        fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        saveas(fig, fullfile(subject_dir,'12-zapped_psd.png'));
        
        % Subject rows in behav data
        exe_subj_rows=find(strcmp(subj_id,exe_df.Subject));
        obs_subj_rows=find(strcmp(subj_id,obs_df.Subject));
        
        % Remove trials with no grasp~= or both hand grasp - only for LEXT,
        % FTGE, EXGC, and EXEND epochs
        other_trials=find(~strcmp(exe_df.FTGhand(exe_subj_rows),'L') & ~strcmp(exe_df.FTGhand(exe_subj_rows),'R'));
        for i=1:length(other_trials)
            t=other_trials(i);
            epochs_to_remove=(([EEG.epoch.eventtrial]==t) & (strcmp({EEG.epoch.eventtype},'LEXT') | strcmp({EEG.epoch.eventtype},'FTGE') | strcmp({EEG.epoch.eventtype},'EXGC') | strcmp({EEG.epoch.eventtype},'EXEND')));
            EEG = pop_rejepoch(EEG,epochs_to_remove, 0);
        end
        
        % Remove trials with two touches - only for FTGE, EXGC, and EXEND
        ft_to_ftg=exe_df.FTG(exe_subj_rows)-exe_df.FT(exe_subj_rows);
        two_touch_trials=find(ft_to_ftg>0);
        for i=1:length(two_touch_trials)
            t=two_touch_trials(i);
            epochs_to_remove=(([EEG.epoch.eventtrial]==t) & (strcmp({EEG.epoch.eventtype},'FTGE') | strcmp({EEG.epoch.eventtype},'EXGC') | strcmp({EEG.epoch.eventtype},'EXEND')));
            EEG = pop_rejepoch(EEG,epochs_to_remove, 0);
        end
        
        % Remove trials with too long grasp duration - only for EXGC and
        % EXEND
        ftg_to_gc=exe_df.GC(exe_subj_rows)-exe_df.FTG(exe_subj_rows);
        long_grasp_trials=find(ftg_to_gc/1000>1.6);
        for i=1:length(long_grasp_trials)
            epochs_to_remove=(([EEG.epoch.eventtrial]==t) & (strcmp({EEG.epoch.eventtype},'EXGC') | strcmp({EEG.epoch.eventtype},'EXEND')));
            EEG = pop_rejepoch(EEG,epochs_to_remove, 0);
        end
        
        EEG = pop_editset(EEG, 'setname', sprintf('%s_task-tool_obs_exe_eeg_processed_data',subj_id));
        pop_saveset(EEG, 'filepath', fullfile(subject_dir, 'processed_data'),...
            'filename', sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id));
        
        fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        saveas(fig, fullfile(subject_dir,'13-final_psd.png'));
        
        for i=1:length(epoch_types)
            epoch_type=epoch_types{i};
            trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,epoch_type)), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
            n_epochs(s,i)=length(trials);
        end
    end
end

disp(sprintf('Total subjects: %d',n_subjects));
for i=1:length(epoch_types)
    epoch_type=epoch_types{i};
    above_thresh=find(n_epochs(:,i)>=5);
    n_above_thresh=length(above_thresh);
    min_n=min(n_epochs(above_thresh,i));
    max_n=max(n_epochs(above_thresh,i));
    mean_n=mean(n_epochs(above_thresh,i));
    sd_n=std(n_epochs(above_thresh,i));
    disp(sprintf('%s: %d subjects, %d-%d, M=%.2f, SD=%.2f', epoch_type, n_above_thresh, min_n, max_n, mean_n, sd_n));
end