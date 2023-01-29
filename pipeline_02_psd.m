function pipeline_02_psd(study_info, varargin)
% Compute power spectrial density in all electrodes

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

pipeline='NEARICA';

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Power spectra for each subject in each cluster
spectra=[];
periodic=[];
aperiodic=[];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

subjects={};

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
    
    % Load data
    fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
    
    if exist(fullfile(subject_data_dir,fname),'file')==2
        EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', fname);        
                
        ch_space=[[EEG.chanlocs.X];[EEG.chanlocs.Y];[EEG.chanlocs.Z]]';
        for i=1:3
            ch_space(:,i)=ch_space(:,i)-min(ch_space(:,i));
            ch_space(:,i)=ch_space(:,i)./max(ch_space(:,i));
        end
        % Use a window size of 1s with a 50% overlap
        winsize=EEG.srate;
        overlap=round(winsize/2);

        % Compute power spectral density for each cluster using Welch's
        % method
        [subj_spectra,frex,~,~,~] = spectopo(EEG.data,...
                EEG.pnts, EEG.srate, 'winsize', winsize,...
                'overlap', overlap, 'plot', 'off', 'freqfac',10);
        freq_idx=find((frex>=2) & (frex<=100));
        frex=frex(freq_idx);
        subj_spectra=subj_spectra(:,freq_idx);

        subj_spectra=subj_spectra./10;
    
        ch_aperiodic=zeros(size(subj_spectra));
        ch_resids=zeros(size(subj_spectra));
        for i=1:size(subj_spectra,1)
            % Channel spectrum
            chan_psd=subj_spectra(i,:);
            % 1/f
            oof=log10(1./frex);
            % Fit 1/f to spectrum
            lm_psd=fitlm(oof,chan_psd,'RobustOpts','on');
            % Get residuals
            ch_resids(i,:)=lm_psd.Residuals.Raw;%-min(lm_psd.Residuals.Raw);                
            ch_aperiodic(i,:)=lm_psd.Fitted;
        end
            
        figure();
        hold all
        for ch=1:length(EEG.chanlocs)
            plot(frex,ch_resids(ch,:),'color',ch_space(ch,:));
        end
        
        subjects{end+1}=subj_id;
        spectra(end+1,:,:)=subj_spectra;    
        periodic(end+1,:,:)=ch_resids;
        aperiodic(end+1,:,:)=ch_aperiodic;
    end
end

save(fullfile(study_info.data_dir,'derivatives', pipeline, 'processed_psd.mat'),'subjects','frex','spectra','periodic','aperiodic');

