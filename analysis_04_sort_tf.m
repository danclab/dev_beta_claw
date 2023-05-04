function analysis_04_sort_tf(study_info, varargin)

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

pipeline='NEARICA_behav';
exe_epoch_types={'EBM','LEXT','FTGE','EXGC','EXEND'};
obs_epoch_types={'OBM','LOBS','FTGO','OBGC','OBEND'};

% behavior
exe_df=readtable(fullfile(study_info.data_dir,'behavior_exe.csv'));

% Number of subjects
n_subjects=size(study_info.participant_info,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    disp(subj_id);
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
    if exist(fullfile(subject_data_dir,fname),'file')==2
        
        % Load data
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);        
        time=EEG.times;
        
        % Subject rows in behav data
        exe_subj_rows=strcmp(subj_id,exe_df.Subject);
        
        for i=1:length(exe_epoch_types)
        
            epoch_type=exe_epoch_types{i};            
            
            % Epochs for this epoch type
            eeg_epochs=find(strcmp({EEG.epoch.eventtype},epoch_type));

            % TF matrices after sorting into ipsi and contralateral
            c_ipsi_tf=[];
            c_contra_tf=[];
            
            % Load C3 and C4 TFs
            load(fullfile(subject_data_dir,sprintf("processed_superlet_tf_C3.mat")));
            c3_trial_tf=trial_tf;
            load(fullfile(subject_data_dir,sprintf("processed_superlet_tf_C4.mat")));
            c4_trial_tf=trial_tf;
            
            % For all epochs of this type
            for t=1:length(eeg_epochs)
                % Index of this epoch in EEG
                eeg_idx=eeg_epochs(t);
                
                % Trial number for this epoch
                trial=EEG.epoch(eeg_idx).eventtrial;
                
                % For for this trial in behavior
                trial_row=exe_subj_rows & (exe_df.Trial==trial);
                
                % Hand used
                hand_used=exe_df.FTGhand(trial_row);
                
                if strcmp(hand_used,'L')
                    c_ipsi_tf(end+1,:,:,:)=c3_trial_tf(eeg_idx,:,:,:);
                    c_contra_tf(end+1,:,:,:)=c4_trial_tf(eeg_idx,:,:,:);
                elseif strcmp(hand_used,'R')
                    c_ipsi_tf(end+1,:,:,:)=c4_trial_tf(eeg_idx,:,:,:);
                    c_contra_tf(end+1,:,:,:)=c3_trial_tf(eeg_idx,:,:,:);
                elseif ~strcmp(epoch_type,'EBM')
                    disp('error');
                end
            end
            if size(c_ipsi_tf,1)>=5
                c_ipsi_tf=squeeze(mean(mean(c_ipsi_tf,2),1));
                c_contra_tf=squeeze(mean(mean(c_contra_tf,2),1));
                save(fullfile(subject_data_dir,sprintf('processed_superlet_tf_%s.mat',epoch_type)),...
                    'foi','c_ipsi_tf','c_contra_tf','time','-v7.3');
            end
        
        end

        for i=1:length(obs_epoch_types)
        
            epoch_type=obs_epoch_types{i};            
            
            % Epochs for this epoch type
            eeg_epochs=find(strcmp({EEG.epoch.eventtype},epoch_type));

            % TF matrices but ipsi=C3 and contra=C4
            c_ipsi_tf=[];
            c_contra_tf=[];
            
            % Load C3 and C4 TFs
            load(fullfile(subject_data_dir,sprintf("processed_superlet_tf_C3.mat")));
            c3_trial_tf=trial_tf;
            load(fullfile(subject_data_dir,sprintf("processed_superlet_tf_C4.mat")));
            c4_trial_tf=trial_tf;
            
            % For all epochs of this type
            for t=1:length(eeg_epochs)
                % Index of this epoch in EEG
                eeg_idx=eeg_epochs(t);
                
                % Trial number for this epoch
                trial=EEG.epoch(eeg_idx).eventtrial;
                
                c_ipsi_tf(end+1,:,:,:)=c3_trial_tf(eeg_idx,:,:,:);
                c_contra_tf(end+1,:,:,:)=c4_trial_tf(eeg_idx,:,:,:);
            end
            if size(c_ipsi_tf,1)>=5
                c_ipsi_tf=squeeze(mean(mean(c_ipsi_tf,2),1));
                c_contra_tf=squeeze(mean(mean(c_contra_tf,2),1));
                save(fullfile(subject_data_dir,sprintf('processed_superlet_tf_%s.mat',epoch_type)),...
                    'foi','c_ipsi_tf','c_contra_tf','time','-v7.3');
            end
        
        end
    end
end
