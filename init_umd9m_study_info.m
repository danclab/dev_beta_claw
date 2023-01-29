%% Individualized_Info
function study_info=init_umd9m_study_info()

study_info=[];

study_info.age='9m';

study_info.task='tool_obs_exe';

%% Specify Paths
% Locate data on hard drive
study_info.data_dir = '/home/bonaiuto/dev_beta_umd/data/9m';

%% Subject list
study_info.participant_info=readtable([fullfile(study_info.data_dir, 'data', 'participants.tsv')],...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{''});

study_info.clusters={'C3','C4','P3','P4'};
study_info.cluster_channels={
    {'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'},...
    {'E26', 'E27', 'E28', 'E31'},...
    {'E40', 'E42', 'E45', 'E46'}};
