function pipeline_05_compute_mean_tf(study_info, varargin)

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

pipeline='NEARICA';

% Number of subjects
n_subjects=size(study_info.participant_info,1);

foi=linspace(1,100,200);

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
        
        exe_base_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'EBM')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        exe_exp_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'FTGE')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        if length(exe_base_trials)>=5 && length(exe_exp_trials)>=5
            load(fullfile(subject_data_dir,tf_fname));

            base_tf=squeeze(trial_tf(exe_base_trials,:,:,:));
            exp_tf=squeeze(trial_tf(exe_exp_trials,:,:,:));

            mean_base=squeeze(mean(mean(base_tf,4),1));
            base_tf=(squeeze(mean(base_tf,1)));
            exp_tf=(squeeze(mean(exp_tf,1)));

            save(fullfile(subject_data_dir,'processed_superlet_exe_mean_tf.mat'),...
                'foi','base_tf','exp_tf','mean_base','time','-v7.3');
        end
        
        obs_base_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'OBM')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        obs_exp_trials=find(cellfun(@length,[cellfun(@(x) find(strcmp(x,'FTGO')), {EEG.epoch.eventtype},'UniformOutput',false)])>0);
        if length(obs_base_trials)>=5 && length(obs_exp_trials)>=5
            load(fullfile(subject_data_dir,tf_fname));

            base_tf=squeeze(trial_tf(obs_base_trials,:,:,:));
            exp_tf=squeeze(trial_tf(obs_exp_trials,:,:,:));

            mean_base=squeeze(mean(mean(base_tf,4),1));
            base_tf=(squeeze(mean(base_tf,1)));
            exp_tf=(squeeze(mean(exp_tf,1)));

            save(fullfile(subject_data_dir,'processed_superlet_obs_mean_tf.mat'),...
                'foi','base_tf','exp_tf','mean_base','time','-v7.3');
        end
    end
end
