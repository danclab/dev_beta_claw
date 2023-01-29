function compute_num_trials_initial(study_info)

% Number of subjects
n_subjects=size(study_info.participant_info,1);

subj_trials=zeros(n_subjects,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    disp(subj_id);
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'data', subj_id, 'eeg');
    
    fname=sprintf('%s_task-tool_obs_exe_eeg.set',subj_id);
    
    % Load data
    EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', fname);        

    n_obs_trials=length(find(strcmp({EEG.event.type},'FTGO')));
    subj_trials(s)=n_obs_trials;
end
disp(sprintf('%d-%d, M = %.2f, SD = %.2f', min(subj_trials), ...
    max(subj_trials), mean(subj_trials), std(subj_trials)));