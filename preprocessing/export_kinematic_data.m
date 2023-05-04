function export_kinematic_data(study_info)

subject={};
condition={};
trial=[];
reach_dur=[];
grasp_dur=[];
end_dur=[];

% behavior
exe_df=readtable(fullfile(study_info.data_dir,'behavior_exe.csv'));
obs_df=readtable(fullfile(study_info.data_dir,'behavior_obs.csv'));

%% Loop over all data files
for s_idx=1:size(study_info.participant_info,1)
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s_idx};
    
    % Where original raw data is located
    subject_raw_data_dir=fullfile(study_info.data_dir, 'data',subj_id, 'eeg');
    
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subj_id);
    
    %% Import data
    data_file_name=sprintf('%s_task-%s_eeg.set',subj_id, study_info.task);
    EEG=pop_loadset('filename', data_file_name, 'filepath', subject_raw_data_dir);
    EEG = eeg_checkset(EEG);
    
    % Remove events after last LBBS event
    lbbs_evts=find(strcmp('LBBS',{EEG.event.type}));
    EEG.event=EEG.event(1:lbbs_evts(end));
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Remove extra events
    evts_to_keep=find(strcmp('LBBS',{EEG.event.type}) | strcmp('LBSE',{EEG.event.type}) | strcmp('LBOB',{EEG.event.type}) | strcmp('LOBS',{EEG.event.type}) | strcmp('FTGO',{EEG.event.type}) | strcmp('LBEX',{EEG.event.type}) | strcmp('LEXT',{EEG.event.type}) | strcmp('FTGE',{EEG.event.type}));
    EEG.event=EEG.event(evts_to_keep);
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    %% Label the task
    % Add observation baseline events
    base_obs_latencies=[EEG.event(find(strcmp('LBOB',{EEG.event.type}))).latency];
    for i=1:length(base_obs_latencies)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='OBM';
        EEG.event(n_events+1).latency=(base_obs_latencies(i)-1.5*EEG.srate)-1;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Add execution baseline events
    exe_obs_latencies=[EEG.event(find(strcmp('LBEX',{EEG.event.type}))).latency];
    for i=1:length(exe_obs_latencies)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='EBM';
        EEG.event(n_events+1).latency=(exe_obs_latencies(i)-1.5*EEG.srate)-1;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    
    exe_subj_rows=find(strcmp(subj_id,exe_df.Subject));
    
    % Add grasp completion and trial end events, and fake FTGE and GC
    % events if no grasp (will be removed later)
    start_evts=find(strcmp('LEXT',{EEG.event.type}));
    start_times=[EEG.event(start_evts).latency];
    ende_evts=find(strcmp('LBBS',{EEG.event.type}));
    end_times=[EEG.event(ende_evts).latency];
    
    for i=1:length(start_evts)
        start_time=start_times(i);
        end_evt_time=end_times(find(end_times>start_time,1));
        
        % Add fake FTG event
        if exe_df.FTG(exe_subj_rows(i))==0
            % put 25% between start and end
            ftg_time=start_time+.25*(end_evt_time-start_time);
            n_events=length(EEG.event);
            EEG.event(n_events+1).type='FTGE';
            EEG.event(n_events+1).latency=ftg_time;
            EEG.event(n_events+1).urevent=n_events+1;
        end
        
        % Add GC event
        if exe_df.GC(exe_subj_rows(i))>0
            gc_time=exe_df.Gcfinal(exe_subj_rows(i));
            n_events=length(EEG.event);
            EEG.event(n_events+1).type='EXGC';
            EEG.event(n_events+1).latency=(gc_time/1000)*EEG.srate;
            EEG.event(n_events+1).urevent=n_events+1;
        % Add fake GC event
        else
            % put 75% between start and end
            gc_time=start_time+.75*(end_evt_time-start_time);
            n_events=length(EEG.event);
            EEG.event(n_events+1).type='EXGC';
            EEG.event(n_events+1).latency=gc_time;
            EEG.event(n_events+1).urevent=n_events+1;
        end
        
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='EXEND';
        EEG.event(n_events+1).latency=end_evt_time;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Set L/R hand for each event
    n_events=length(EEG.event);
    trial_vals=zeros(1, n_events);
    EEG=pop_editeventfield(EEG, 'trial',trial_vals);
    base_evts=find(strcmp('EBM',{EEG.event.type}));
    start_evts=find(strcmp('LEXT',{EEG.event.type}));
    ftge_evts=find(strcmp('FTGE',{EEG.event.type}));
    gc_evts=find(strcmp('EXGC',{EEG.event.type}));
    ende_evts=find(strcmp('EXEND',{EEG.event.type}));    
    for i=1:length(start_evts)
        hand=exe_df.FTGhand{exe_subj_rows(i)};
        EEG.event(base_evts(i)).value=hand;
        EEG.event(base_evts(i)).trial=i;
        EEG.event(start_evts(i)).value=hand;
        EEG.event(start_evts(i)).trial=i;
        EEG.event(ftge_evts(i)).value=hand;
        EEG.event(ftge_evts(i)).trial=i;
        EEG.event(gc_evts(i)).value=hand;
        EEG.event(gc_evts(i)).trial=i;
        EEG.event(ende_evts(i)).value=hand;
        EEG.event(ende_evts(i)).trial=i;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    obs_subj_rows=find(strcmp(subj_id,obs_df.Subject));
    
    % Add grasp completion and trial end events
    start_evts=find(strcmp('LOBS',{EEG.event.type}));
    start_times=[EEG.event(start_evts).latency];
    ende_evts=find(strcmp('LBBS',{EEG.event.type}));
    end_times=[EEG.event(ende_evts).latency];
    
    for i=1:length(start_evts)
        start_time=start_times(i);
        end_evt_time=end_times(find(end_times>start_time,1));
        
        % Add GC event
        gc_time=obs_df.Gcfinal(obs_subj_rows(i));
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='OBGC';
        EEG.event(n_events+1).latency=(gc_time/1000)*EEG.srate;
        EEG.event(n_events+1).urevent=n_events+1;
        
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='OBEND';
        EEG.event(n_events+1).latency=end_evt_time;
        EEG.event(n_events+1).urevent=n_events+1;
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Set L/R hand for each event
    base_evts=find(strcmp('OBM',{EEG.event.type}));
    start_evts=find(strcmp('LOBS',{EEG.event.type}));
    ftgo_evts=find(strcmp('FTGO',{EEG.event.type}));
    gc_evts=find(strcmp('OBGC',{EEG.event.type}));
    ende_evts=find(strcmp('OBEND',{EEG.event.type}));    
    for i=1:length(start_evts)
        if i<=length(base_evts)
            EEG.event(base_evts(i)).trial=i;
            EEG.event(base_evts(i)).value='R';
        end
        EEG.event(start_evts(i)).trial=i;
        EEG.event(start_evts(i)).value='R';
        EEG.event(ftgo_evts(i)).trial=i;
        EEG.event(ftgo_evts(i)).value='R';
        EEG.event(gc_evts(i)).trial=i;
        EEG.event(gc_evts(i)).value='R';
        EEG.event(ende_evts(i)).trial=i;
        EEG.event(ende_evts(i)).value='R';
    end
    %check for consistency and reorder the events chronologically...
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    %% Delete discontinuous data from the raw data file
    % remove data after last task event
    end_flags=find(strcmp({EEG.event.type},'LBBS'));
    latency=EEG.event(end_flags(end)).latency;
    % remove everything 1.5 seconds after the last event
    EEG = eeg_eegrej( EEG, [(latency+(1.5*EEG.srate)) EEG.pnts] );
    EEG = eeg_checkset( EEG );
    
    lbse_flags = find(strcmp({EEG.event.type},'LBSE'));
    latency=EEG.event(lbse_flags(1)).latency;
    % remove everything until 1.5 seconds before the first event
    EEG = eeg_eegrej( EEG, [1 (latency-(1.5*EEG.srate))] );
    EEG = eeg_checkset( EEG );
    
    % Remove extra events
    evts_to_keep=find(strcmp('OBM',{EEG.event.type}) | strcmp('EBM',{EEG.event.type}) | strcmp('LOBS',{EEG.event.type}) | strcmp('FTGO',{EEG.event.type}) | strcmp('OBGC',{EEG.event.type}) | strcmp('OBEND',{EEG.event.type}) | strcmp('LEXT',{EEG.event.type}) | strcmp('FTGE',{EEG.event.type}) | strcmp('EXGC',{EEG.event.type}) | strcmp('EXEND',{EEG.event.type}));
    EEG.event=EEG.event(evts_to_keep);
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    evts_to_keep=find(strcmp('L',{EEG.event.value}) | strcmp('R',{EEG.event.value}));
    EEG.event=EEG.event(evts_to_keep);
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    ex_start_evts=find(strcmp('LEXT',{EEG.event.type}));
    ex_start_latency=[EEG.event(ex_start_evts).latency];    
    ftge_evts=find(strcmp('FTGE',{EEG.event.type}));
    ftge_latency=[EEG.event(ftge_evts).latency];    
    exe_reach_dur=(ftge_latency-ex_start_latency)/EEG.srate;
    
    grspe_evts=find(strcmp('EXGC',{EEG.event.type}));
    grspe_latency=[EEG.event(grspe_evts).latency];    
    exe_grsp_dur=(grspe_latency-ftge_latency)/EEG.srate;
    
    ende_evts=find(strcmp('EXEND',{EEG.event.type}));
    ende_latency=[EEG.event(ende_evts).latency];    
    exe_end_dur=(ende_latency-grspe_latency)/EEG.srate;
    
    exe_trials=[EEG.event(ex_start_evts).trial];
    for j=1:length(exe_trials)
        if exe_reach_dur(j)>=0 && exe_grsp_dur(j)>=0 && exe_end_dur(j)>=0 && exe_grsp_dur(j)<=1.6
            subject{end+1}=subj_id;
            condition{end+1}='exe';
            trial(end+1)=exe_trials(j);
            reach_dur(end+1)=exe_reach_dur(j);
            grasp_dur(end+1)=exe_grsp_dur(j);
            end_dur(end+1)=exe_end_dur(j);
        end
    end
    
    ob_start_evts=find(strcmp('LOBS',{EEG.event.type}));
    ob_start_latency=[EEG.event(ob_start_evts).latency];    
    ftgo_evts=find(strcmp('FTGO',{EEG.event.type}));
    ftgo_latency=[EEG.event(ftgo_evts).latency];    
    obs_reach_dur=(ftgo_latency-ob_start_latency)/EEG.srate;
    
    grspo_evts=find(strcmp('OBGC',{EEG.event.type}));
    grspo_latency=[EEG.event(grspo_evts).latency];    
    obs_grsp_dur=(grspo_latency-ftgo_latency)/EEG.srate;
    
    endo_evts=find(strcmp('OBEND',{EEG.event.type}));
    endo_latency=[EEG.event(endo_evts).latency];    
    obs_end_dur=(endo_latency-grspo_latency)/EEG.srate;
    
    obs_trials=[EEG.event(ob_start_evts).trial];
    for j=1:length(obs_trials)
        if obs_reach_dur(j)>=0 && obs_grsp_dur(j)>=0 && obs_end_dur(j)>=0 && obs_grsp_dur(j)<=1.6
            subject{end+1}=subj_id;
            condition{end+1}='obs';
            trial(end+1)=obs_trials(j);
            reach_dur(end+1)=obs_reach_dur(j);
            grasp_dur(end+1)=obs_grsp_dur(j);
            end_dur(end+1)=obs_end_dur(j);
        end
    end
end
df=table(subject', condition', trial', reach_dur', grasp_dur', end_dur', 'VariableNames', {'Subject','Condition','Trial','ReachDur','GraspDur','EndDur'});
writetable(df, fullfile(study_info.data_dir, 'derivatives', 'NEARICA_behav', 'processed_kinematics.csv'));
