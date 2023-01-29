function compute_num_epochs_final(study_info)

pipeline='NEARICA';

% Number of subjects
n_subjects=size(study_info.participant_info,1);

subj_exe_trials=zeros(n_subjects,1);
subj_obs_trials=zeros(n_subjects,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    disp(subj_id);
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
    
    % Load data
    EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', fname);        

    n_exe_trials=length(find(strcmp({EEG.event.type},'FTGE')));
    n_obs_trials=length(find(strcmp({EEG.event.type},'FTGO')));
    subj_exe_trials(s)=n_exe_trials;
    subj_obs_trials(s)=n_obs_trials;
end
disp('Execution');
disp(sprintf('%d-%d, M = %.2f, SD = %.2f', min(subj_exe_trials), ...
    max(subj_exe_trials), mean(subj_exe_trials), std(subj_exe_trials)));
disp('Observation');
disp(sprintf('%d-%d, M = %.2f, SD = %.2f', min(subj_obs_trials), ...
    max(subj_obs_trials), mean(subj_obs_trials), std(subj_obs_trials)));