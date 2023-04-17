function compute_num_epochs_final(study_info)

pipeline='NEARICA_behav';

% Number of subjects
n_subjects=size(study_info.participant_info,1);

subj_lext_trials=zeros(n_subjects,1);
subj_ftge_trials=zeros(n_subjects,1);
subj_exgc_trials=zeros(n_subjects,1);
subj_exend_trials=zeros(n_subjects,1);

subj_lobs_trials=zeros(n_subjects,1);
subj_ftgo_trials=zeros(n_subjects,1);
subj_obgc_trials=zeros(n_subjects,1);
subj_obend_trials=zeros(n_subjects,1);

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

    n_lext_trials=length(find(strcmp({EEG.event.type},'LEXT')));
    n_ftge_trials=length(find(strcmp({EEG.event.type},'FTGE')));
    n_exgc_trials=length(find(strcmp({EEG.event.type},'EXGC')));
    n_exend_trials=length(find(strcmp({EEG.event.type},'EXEND')));
    subj_lext_trials(s)=n_lext_trials;
    subj_ftge_trials(s)=n_ftge_trials;
    subj_exgc_trials(s)=n_exgc_trials;
    subj_exend_trials(s)=n_exend_trials;
    
    n_lobs_trials=length(find(strcmp({EEG.event.type},'LOBS')));
    n_ftgo_trials=length(find(strcmp({EEG.event.type},'FTGO')));
    n_obgc_trials=length(find(strcmp({EEG.event.type},'OBGC')));
    n_obend_trials=length(find(strcmp({EEG.event.type},'OBEND')));
    subj_lobs_trials(s)=n_lobs_trials;
    subj_ftgo_trials(s)=n_ftgo_trials;
    subj_obgc_trials(s)=n_obgc_trials;
    subj_obend_trials(s)=n_obend_trials;
end
disp('Execution');
idx=find(subj_lext_trials>=5);
subj_lext_trials=subj_lext_trials(idx);
disp(sprintf('LEXT: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_lext_trials), ...
    max(subj_lext_trials), mean(subj_lext_trials), std(subj_lext_trials)));
idx=find(subj_ftge_trials>=5);
subj_ftge_trials=subj_ftge_trials(idx);
disp(sprintf('FTGE: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_ftge_trials), ...
    max(subj_ftge_trials), mean(subj_ftge_trials), std(subj_ftge_trials)));
idx=find(subj_exgc_trials>=5);
subj_exgc_trials=subj_exgc_trials(idx);
disp(sprintf('EXGC: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_exgc_trials), ...
    max(subj_exgc_trials), mean(subj_exgc_trials), std(subj_exgc_trials)));
idx=find(subj_exend_trials>=5);
subj_exend_trials=subj_exend_trials(idx);
disp(sprintf('EXEND: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_exend_trials), ...
    max(subj_exend_trials), mean(subj_exend_trials), std(subj_exend_trials)));

disp('Observation');
idx=find(subj_lobs_trials>=5);
subj_lobs_trials=subj_lobs_trials(idx);
disp(sprintf('LOBS: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_lobs_trials), ...
    max(subj_lobs_trials), mean(subj_lobs_trials), std(subj_lobs_trials)));
idx=find(subj_ftgo_trials>=5);
subj_ftgo_trials=subj_ftgo_trials(idx);
disp(sprintf('FTGO: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_ftgo_trials), ...
    max(subj_ftgo_trials), mean(subj_ftgo_trials), std(subj_ftgo_trials)));
idx=find(subj_obgc_trials>=5);
subj_obgc_trials=subj_obgc_trials(idx);
disp(sprintf('OBGC: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_obgc_trials), ...
    max(subj_obgc_trials), mean(subj_obgc_trials), std(subj_obgc_trials)));
idx=find(subj_obend_trials>=5);
subj_obend_trials=subj_obend_trials(idx);
disp(sprintf('OBEND: %d subjects, %d-%d, M = %.2f, SD = %.2f', length(idx), min(subj_obend_trials), ...
    max(subj_obend_trials), mean(subj_obend_trials), std(subj_obend_trials)));