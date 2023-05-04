function [lags,foi,lagged_coh]=load_lagged_coherence(study_info, pipeline, varargin)

% Parse optional arguments
defaults=struct('foi',[5 80]);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Load lagged coherence
load(fullfile(study_info.data_dir,'derivatives',pipeline,'processed_lagged_coherence.mat'));
freq_idx=find((foi>=params.foi(1)) & (foi<=params.foi(2)));
foi=foi(freq_idx);
lagged_coh=lagged_coh(:,:,freq_idx,:);

for s=1:size(lagged_coh,1)
    subj_lagged_coh=squeeze(lagged_coh(s,:,:,:));
    lagged_coh(s,:,:,:)=subj_lagged_coh./max(subj_lagged_coh(:));
end
