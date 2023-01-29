function check_event_timing(study_info)

all_obs_base=[];
all_obs_reach=[];
all_obs_grasp=[];
all_exe_base=[];
all_exe_reach=[];
all_exe_grasp=[];

for s_idx=1:size(study_info.participant_info,1)
    
    % Get subject ID from study info
    subject=study_info.participant_info.participant_id{s_idx};
    
    % Where original raw data is located
    subject_raw_data_dir=fullfile(study_info.data_dir, 'data',subject, 'eeg');
    data_file_name=sprintf('%s_task-%s_eeg.set',subject, study_info.task);
    EEG=pop_loadset('filename', data_file_name, 'filepath', subject_raw_data_dir);
    
    obs_base_start=find(strcmp({EEG.event.type},'LBOB'))-1;
    obs_base_end=find(strcmp({EEG.event.type},'LBOB'));
    obs_base_duration=[EEG.event(obs_base_end).latency]./EEG.srate-[EEG.event(obs_base_start).latency]./EEG.srate;
    all_obs_base(end+1:end+length(obs_base_duration))=obs_base_duration;
    
    obs_start=find(strcmp({EEG.event.type},'FTGO'))-1;
    touch_start=find(strcmp({EEG.event.type},'FTGO'));
    obs_reach_duration=[EEG.event(touch_start).latency]./EEG.srate-[EEG.event(obs_start).latency]./EEG.srate;
    all_obs_reach(end+1:end+length(obs_reach_duration))=obs_reach_duration;
    
    obs_end=find(strcmp({EEG.event.type},'FTGO'))+1;
    obs_grasp_duration=[EEG.event(obs_end).latency]./EEG.srate-[EEG.event(touch_start).latency]./EEG.srate;
    all_obs_grasp(end+1:end+length(obs_grasp_duration))=obs_grasp_duration;
    
    exe_base_start=find(strcmp({EEG.event.type},'LBEX'))-1;
    exe_base_end=find(strcmp({EEG.event.type},'LBEX'));
    exe_base_duration=[EEG.event(exe_base_end).latency]./EEG.srate-[EEG.event(exe_base_start).latency]./EEG.srate;
    all_exe_base(end+1:end+length(exe_base_duration))=exe_base_duration;
    
    exe_start=find(strcmp({EEG.event.type},'FTGE'))-1;
    touch_start=find(strcmp({EEG.event.type},'FTGE'));
    exe_reach_duration=[EEG.event(touch_start).latency]./EEG.srate-[EEG.event(exe_start).latency]./EEG.srate;
    all_exe_reach(end+1:end+length(exe_reach_duration))=exe_reach_duration;
    
    exe_end=find(strcmp({EEG.event.type},'FTGE'))+1;
    exe_grasp_duration=[EEG.event(exe_end).latency]./EEG.srate-[EEG.event(touch_start).latency]./EEG.srate;
    all_exe_grasp(end+1:end+length(exe_grasp_duration))=exe_grasp_duration;
end
    
print('');
