function analysis_01_lagged_coherence(study_info, varargin)
% Computes lagged coherence across a range of frequencies and lags

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

pipeline='NEARICA_behav';
n_chans=64;

% Frequencies to run on. We start at 5Hz because the trials aren't long
% enough to look at lower frequencies over long lags
foi=[5:1:100];

% Lags to run on
lags=[2:.1:4.5];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Lagged coherence for each subject in each cluster over all frequencies
% and lags
lagged_coh=[];

s_idx=1;

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
        
        if EEG.trials>=5
            % Create fieldtrip data structure
            data=create_ft_data({EEG.chanlocs.labels}, EEG.data, EEG.times, EEG.srate);                        
            fsample=data.fsample;

            % Run lagged coherence
            l_lagged_coh=zeros(length(data.label),length(foi));

            for l_idx=1:length(lags)
                lag=lags(l_idx);                

                parfor f_idx = 1:length(foi)

                    % Configuration for frequency analysis
                    cfg_F             = [];
                    cfg_F.method      = 'mtmconvol';
                    cfg_F.taper       = 'hanning';
                    cfg_F.output      = 'fourier';
                    cfg_F.keeptrials  = 'yes';
                    cfg_F.pad         = 'nextpow2';

                    % Configuration for lagged coherence
                    cfg_LC            = [];            
                    cfg_LC.method     = 'laggedcoherence';
                    cfg_LC.trialsets  = 'all';
                    fs                = fsample;            
                    cfg_LC.lag        = lag;
                    cfg_F.width       = lag;

                    % Set freq range
                    cfg_F.foi     = foi(f_idx);
                    cfg_LC.foi    = foi(f_idx);

                    % width of time windows for phase comparison (seconds)
                    width         = cfg_F.width/cfg_F.foi;                        
                    cfg_F.t_ftimwin = width;

                    % half width of time window (seconds)
                    halfwidth     = ceil(fs*width/2)/fs;

                    % Go from half window width after trial start to half window
                    % width before trial end
                    toi_start     = data.time{1}(1) + halfwidth;
                    toi_stop      = data.time{1}(end) - halfwidth;

                    % Step size
                    step          = ceil(fs*cfg_LC.lag/cfg_F.foi)/fs;
                    cfg_F.toi     = toi_start:step:toi_stop;

                    % Run frequency analysis
                    freqout       = ft_freqanalysis(cfg_F,data);

                    % Compute lagged coherence
                    lcoh = ft_connectivity_laggedcoherence(cfg_LC,freqout);

                    l_lagged_coh(:,f_idx)=lcoh.laggedcoh;
                end  
                lagged_coh(s_idx,:,:,l_idx)=l_lagged_coh;
                s_idx=s_idx+1;
            end    
            save(fullfile(study_info.data_dir,'derivatives', pipeline, 'processed_lagged_coherence.mat'),'lags','foi','lagged_coh');
        end
    end
end

save(fullfile(study_info.data_dir,'derivatives', pipeline, 'processed_lagged_coherence.mat'),'lags','foi','lagged_coh');

