function fois=pipeline_03_determine_fois(study_info, clusters, clus_name, varargin)

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

pipeline='NEARICA';

subj_id=study_info.participant_info.participant_id{1};
subject_data_dir=fullfile(study_info.data_dir, 'derivatives', pipeline, subj_id, 'processed_data');
fname=sprintf('%s_task-tool_obs_exe_eeg_processed_data.set',subj_id);
EEG=pop_loadset('filepath', subject_data_dir, 'filename', fname);

% Load power spectral densities
load(fullfile(study_info.data_dir,'derivatives', pipeline, 'processed_psd.mat'));

% Look at residuals from 5-40Hz because we don't have lagged coherence at
% lower frequencies (trials aren't long enough)
freq_idx=find(frex>=1);
frex=frex(freq_idx);
periodic=periodic(:,:,freq_idx);

fwhm_min=2;

fois=[];

clus_chan_idx=[];
for c=1:length(clusters)
    cluster=clusters{c};
    cluster_idx=find(strcmp(study_info.clusters,cluster));
    chan_idx=cellfun(@(x) find(strcmp({EEG.chanlocs.labels},x)),...
        study_info.cluster_channels{cluster_idx});
    clus_chan_idx(end+1:end+length(chan_idx))=chan_idx;
end

figure();

mean_resids=squeeze(nanmean(nanmean(periodic(:,clus_chan_idx,:),2),1));
mean_resids=mean_resids-min(mean_resids);
plot(frex,mean_resids);
hold all

iter=1;
while true
    thresh=std(mean_resids);
    % Find peaks and sort in descending order
    [pks,locs]=peak_finder(mean_resids,frex);
    [sorted_pks,sort_idx]=sort(pks,'descend');
    sorted_locs=locs(sort_idx);

    % Try first peak
    pk_freq=sorted_locs(1);
    pk_idx=find(frex==pk_freq);
    pk_pow=mean_resids(pk_idx);

    if pk_pow<thresh
        break
    end

    % Find fwhm
    l_idx=find(mean_resids(1:pk_idx)<pk_pow*.5);
    r_idx=find(mean_resids(pk_idx+1:end)<pk_pow*.5);
    if length(l_idx) && length(r_idx)
        l_freq=frex(l_idx(end));
        r_freq=frex(pk_idx-1+r_idx(1));
        r_side=(r_freq-pk_freq);
        l_side=(pk_freq-l_freq);
        fwhm=2*min([r_side l_side]);
    elseif length(l_idx)
        l_freq=frex(l_idx(end));
        fwhm=2*(pk_freq-l_freq);
    elseif length(r_idx)
        r_freq=frex(pk_idx-1+r_idx(1));
        fwhm=2*(r_freq-pk_freq);
    end
    % Band limits
    l_freq=pk_freq-fwhm*.5;
    r_freq=pk_freq+fwhm*.5;
    sd=fwhm/(2*sqrt(2*log(2)));
    A=pk_pow*exp(-.5*((frex-pk_freq)/sd).^2);

    % If within range and not obscuring another peak
    if l_freq>=frex(1) && r_freq<=frex(end)
        inside=false;
        for f=1:iter-1
            if (l_freq>=fois(f).range(1) && r_freq<=fois(f).range(2))
               inside=true;
               break
            end
        end
        if ~inside && fwhm>=fwhm_min
            plot([pk_freq pk_freq],[mean_resids(pk_idx) mean_resids(pk_idx)],'or');
            plot(frex,A);

            fois(iter).range=[l_freq r_freq];
            fois(iter).pk_pow=pk_pow;
            fois(iter).pk_freq=pk_freq;
            fois(iter).fwhm=fwhm;
            disp(sprintf('peak pow=%.4f, peak freq=%.2fHz, fwhm=%.2fHz, range=%.2f-%.2fHz', pk_pow, pk_freq, fwhm, l_freq, r_freq));
            iter=iter+1;
        end
    end

    % Subtract Gaussian from residuals
    mean_resids=mean_resids-A;
    plot(frex,mean_resids,'k--');
end
save(fullfile(study_info.data_dir, 'derivatives', pipeline, sprintf('processed_%s_fois.mat',clus_name)), 'fois');