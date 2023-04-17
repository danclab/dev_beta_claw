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
    
    base_starts=find(strcmp({EEG.event.type},'LBSE'));
    for i=1:length(base_starts)
        base_start=EEG.event(base_starts(i)).latency/EEG.srate;
        next_event=base_starts(i)+1;
        if next_event<=length(EEG.event)
            if strcmp(EEG.event(next_event).type,'LBOB')
                obs_base_end=EEG.event(next_event).latency/EEG.srate;
                obs_base_duration=obs_base_end-base_start;
                all_obs_base(end+1)=obs_base_duration;
            elseif strcmp(EEG.event(next_event).type,'LBEX')
                exe_base_end=EEG.event(next_event).latency/EEG.srate;
                exe_base_duration=exe_base_end-base_start;
                all_exe_base(end+1)=exe_base_duration;
            end            
        end
    end
    disp(sprintf('Baseline epochs: %d', length(base_starts)));
    
    obs_starts=find(strcmp({EEG.event.type},'LOBS'));
    for i=1:length(obs_starts)
        obs_start=EEG.event(obs_starts(i)).latency/EEG.srate;
        next_event=obs_starts(i)+1;
        if next_event<=length(EEG.event)
            if strcmp(EEG.event(next_event).type,'FTGO')
                obs_reach_end=EEG.event(next_event).latency/EEG.srate;
                obs_reach_duration=obs_reach_end-obs_start;
                all_obs_reach(end+1)=obs_reach_duration;
                
                next_next_event=next_event+1;
                obs_grasp_end=EEG.event(next_next_event).latency/EEG.srate;
                obs_grasp_duration=obs_grasp_end-obs_reach_end;
                all_obs_grasp(end+1)=obs_grasp_duration;                
            else
                all_obs_reach(end+1)=NaN;
                all_obs_grasp(end+1)=NaN;
            end
        else
            all_obs_reach(end+1)=NaN;
            all_obs_grasp(end+1)=NaN;
        end
    end
    obs_without_contact=length(obs_starts)-length(find(strcmp({EEG.event.type},'FTGO')));
    disp(sprintf('Observation epochs: %d', length(obs_starts)));
    disp(sprintf('Observation epochs without contact: %d', obs_without_contact));
    
    exe_starts=find(strcmp({EEG.event.type},'LEXT'));
    for i=1:length(exe_starts)
        exe_start=EEG.event(exe_starts(i)).latency/EEG.srate;
        next_event=exe_starts(i)+1;
        if next_event<length(EEG.event)
            if strcmp(EEG.event(next_event).type,'FTGE')
                exe_reach_end=EEG.event(next_event).latency/EEG.srate;
                exe_reach_duration=exe_reach_end-exe_start;
                all_exe_reach(end+1)=exe_reach_duration;
                
                next_next_event=next_event+1;
                exe_grasp_end=EEG.event(next_next_event).latency/EEG.srate;
                exe_grasp_duration=exe_grasp_end-exe_reach_end;
                all_exe_grasp(end+1)=exe_grasp_duration;                
            else
                all_exe_reach(end+1)=NaN;
                all_exe_grasp(end+1)=NaN;
            end            
        end
    end
    exe_without_contact=length(exe_starts)-length(find(strcmp({EEG.event.type},'FTGE')));
    disp(sprintf('Execution epochs: %d', length(exe_starts)));
end
    
figure();hist(all_exe_base,20);
figure();hist(all_obs_base,20);
figure();hist(all_exe_reach,20);
figure();hist(all_exe_grasp,20);
figure();hist(all_obs_reach,20);
figure();hist(all_obs_grasp,20);
