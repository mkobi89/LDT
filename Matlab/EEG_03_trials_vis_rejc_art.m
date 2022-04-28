%% EEG read in data from matlab preprocessed with automagic (mara) %%

%% Preparation:
clear

%% Change directory to the current m. file
% by hand or by code:
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Path to server depending on OS
% if IsOSX == true
%     dataPath = '/Volumes/CLINT/LDT_segmented_preprocessed/'; % All Data Path
%     savePath = '/Volumes/CLINT/LDT_segmented_preprocessed/'; % LDT Path
%     addpath('/Volumes/CLINT/fieldtrip-20200607/'); % add fieldtrip path
%     addpath('/Volumes/CLINT/eeglab2020_0/'); % add eeglab path
%
% else
%     dataPath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT_segmented_preprocessed/'; % All Data Path
%     savePath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT_segmented_preprocessed/'; % LDT Path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/fieldtrip-20200607/'); % add fieldtrip path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/eeglab2020_0/'); % add eeglab path
% end

dataPath = 'F:\LDT_final_2_preprocessed\/'; % All Data Path
savePath = 'F:\LDT_final_2_preprocessed\/'; % LDT Path

% Add EEGLAB and Fieldtrip
addpath('F:\CLINT backup_15.02.2022\fieldtrip-20200607/') % add fieldtrip path
% addpath('//130.60.235.123/users/neuro/Desktop/CLINT/eeglab2020_0/'); % add eeglab path

%% Starting EEGLAB and fieldtrip


diary([dataPath 'settings_EEG_03_erp_ga.txt'])
diary on

ft_defaults;

diary off

%% load files for headmodeling

load('headmodeling/eeglabchans');
load('headmodeling/elec_aligned');
load('headmodeling/labels105');
load('headmodeling/lay129_head');


%% change CD to data and get list of all subjects

cd(dataPath);

vplist = dir('C*');
vpNames = {vplist.name};

cd(fileparts(tmp.Filename));

%% use trialfun_LDT_automagic to define trials and calculate erp


for zz = 1:length(vpNames)
    %     if (exist([savePath vpNames{zz} '/'  vpNames{zz} '_LDT_mara_preprocessed.mat'])==2 && ~(exist([vpNames{zz} '_LDT_EEG_fieldtrip.mat'])==2) && ~(exist([vpNames{zz} '_LDT_EEG_preprocessed.mat'])==2))
    
    load([dataPath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_18.mat'])
    
    % define trials
    cfg = [];
    cfg.dataset = [dataPath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_18.mat'];
    cfg.trialfun  = 'ft_trialfun_LDT_automagic';
    cfg.trialdef.prestim    = 0.5;
    cfg.trialdef.poststim   = 1;
    cfg = ft_definetrial(cfg);
    
    dataorig = ft_redefinetrial(cfg,data);
    
    %% average reference
    cfg = [];
    cfg.channel = 'all'; % this is the default
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    cfg.refmethod = 'avg';
    cfg.implicitref = 'Cz';
    
    datapreprocessed = ft_preprocessing(cfg,dataorig);
    
    
    %% eeg detrend
    cfg.detrend = 'yes';
    datapreprocessed = ft_preprocessing(cfg,datapreprocessed);
    
    
     %% Select 105 eeg channels
    cfg = [];
    cfg.channel = {'E7', 'E106', 'E80', 'E55', 'E31', 'Cz'};
    datapreprocessed = ft_selectdata(cfg,datapreprocessed);
    
    
    %% Visual artefact rejection
    disp(zz)

    cfg          = [];
    cfg.method   = 'summary';
    datapreprocessed_all   = ft_rejectvisual(cfg,datapreprocessed);
    
    save([savePath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_vis_rej_all.mat'], 'datapreprocessed_all', '-v7.3');
    
    
    cfg          = [];
    cfg.method   = 'summary';
    cfg.channel = {'E7', 'E106', 'E80', 'E55', 'E31', 'Cz'};
    cfg.keepchannel = 'yes';
    datapreprocessed_red   = ft_rejectvisual(cfg,datapreprocessed);


    save([savePath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_vis_rej_red.mat'], 'datapreprocessed_red', '-v7.3');
   
end
