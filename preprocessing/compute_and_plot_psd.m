function fig=compute_and_plot_psd(EEG,ch_idx)
    % Use a window size of twice the sampling rate with a 50%
    % overlap
%     winsize=EEG.srate;%*2;
%     overlap=round(winsize/2);
% 
%     % Compute power spectral density for each cluster using Welch's
%     % method
%     [subj_spectra,frex,~,~,~] = spectopo(EEG.data,...
%             EEG.pnts, EEG.srate, 'winsize', winsize,...
%             'overlap', overlap, 'plot', 'off', 'nfft',2*EEG.pnts,...
%             'percent',perc);
    ch_space=[[EEG.chanlocs.X];[EEG.chanlocs.Y];[EEG.chanlocs.Z]]';
    for i=1:3
        ch_space(:,i)=ch_space(:,i)-min(ch_space(:,i));
        ch_space(:,i)=ch_space(:,i)./max(ch_space(:,i));
    end    

    fig=figure();
    hold all
    for ch=1:length(ch_idx)
        [pxx,f]=nt_spect_plot(EEG.data(ch_idx(ch),:)'/sqrt(mean(EEG.data(ch_idx(ch),:)'.^2)),1024,[],[],500);
        freq_idx=(f>=1) & (f<=100);
        plot(f(freq_idx),log(abs(pxx(freq_idx))),'color',ch_space(ch_idx(ch),:),'LineWidth',1);
    end
    xlabel('Frequency (Hz)');
    ylabel('log(power)');
    axes('Position',[.65 0.65 .2 .2]);
    plot_sensors(zeros(1,75),EEG.chanlocs,'electrodes','on',...
        'style','blank','emarker','o','electcolor',ch_space,'whitebk','on');
end