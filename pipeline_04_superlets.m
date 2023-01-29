function pipeline_04_superlets(study_info, varargin)

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
    
    if exist(fullfile(subject_data_dir,fname),'file')==2
        
        % Load data
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);        
        
        % If min number of trials in baseline and experimental epochs
        trial_tf=zeros(EEG.trials, length(foi), EEG.nbchan, EEG.pnts);
        for c=1:EEG.nbchan
            for t=1:EEG.trials                
                trial_tf(t,:,c,:)=faslt(double(squeeze(EEG.data(c,:,t))), EEG.srate, foi([1 end]), 200, 4, [1 40], 1, 1);                
            end
        end

        save(fullfile(subject_data_dir,'derivatives',pipeline,'processed_superlet_tf.mat'),'foi','trial_tf','-v7.3');
    end
end