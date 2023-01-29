function apply_zapline(study_info, pipeline, step)

addpath('NoiseTools');

% Number of subjects
n_subjects=size(study_info.participant_info,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id);
    subject_data_dir=fullfile(subject_dir, step);
    
    fname=sprintf('%s_task-tool_obs_exe_eeg_rereferenced_data.set',subj_id);
    
    if exist(fullfile(subject_data_dir,fname),'file')==2
        
        % Load data
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);   
        
        load(fullfile(subject_dir,'zapline.mat'));
        data=data.*1e6;
        EEG.data=permute(data,[2 3 1]);
        EEG = pop_editset(EEG, 'setname', sprintf('%s_task-tool_obs_exe_eeg_processed_NEAR_data',subj_id));
        pop_saveset(EEG, 'filepath', fullfile(subject_dir, 'processed_data'),...
            'filename', sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id));
        
        fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        saveas(fig, fullfile(subject_dir,'final-zapped_psd.png'));
    end
end


