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

dataPath = 'C:\Users\neuro\Desktop\CLINT\LDT_preprocessed/'; % All Data Path
savePath = 'C:\Users\neuro\Desktop\CLINT/LDT_preprocessed/'; % LDT Path

% Add EEGLAB and Fieldtrip
addpath('C:\Users\neuro\Desktop\CLINT/fieldtrip-20200607/') % add fieldtrip path
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

%% calculate erp


jj = 1;
for zz = 1:length(vpNames)
    %     if (exist([savePath vpNames{zz} '/'  vpNames{zz} '_LDT_mara_preprocessed.mat'])==2 && ~(exist([vpNames{zz} '_LDT_EEG_fieldtrip.mat'])==2) && ~(exist([vpNames{zz} '_LDT_EEG_preprocessed.mat'])==2))
    
    disp(jj)
    
    load([dataPath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_vis_rej_red.mat'])
    
    datapreprocessed = datapreprocessed_red;
    
    if ~(vpNames{zz}== "C50") && ~(vpNames{zz}== "C55") &&  ~(vpNames{zz}== "CE5") && ~(vpNames{zz}== "CM7")
        
        
        %% compute ERP
        % all_correct
        cfg = [];
        
        ind_all_correct = find(datapreprocessed.trialinfo(:,2)== 150 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 151 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 152 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)== 250 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 251 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 252 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)==301 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==302 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==304 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==305 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==307 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==308 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)==303 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==306 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==309 & datapreprocessed.trialinfo(:,4)==1);
        
        if length(ind_all_correct) ~= 0
            
            cfg = [];
            cfg.trials = ind_all_correct;
            erp_all_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_all_correct = ft_selectdata(cfg,erp_all_correct);
            
            all_erp_all_correct{1,jj} = vpNames{zz};
            all_erp_all_correct{2,jj} = erp_all_correct;
            all_erp_all_correct{3,jj} = length(ind_all_correct);
            
        end
        % Word vs Nonword all tasks
        
        cfg = [];
        
        ind_word_all_task = find(datapreprocessed.trialinfo(:,2)== 150 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 151 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)== 250 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 251 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)==301 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==302 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==304 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==305 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==307 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==308 & datapreprocessed.trialinfo(:,4)==1);
        
        if length(ind_word_all_task) ~= 0
            cfg = [];
            cfg.trials = ind_word_all_task;
            erp_word_all_task = ft_timelockanalysis(cfg,datapreprocessed);
            
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_word_all_task = ft_selectdata(cfg,erp_word_all_task);
            
            all_erp_word_all_task{1,jj} = vpNames{zz};
            all_erp_word_all_task{2,jj} = erp_word_all_task;
            all_erp_word_all_task{3,jj} = length(ind_word_all_task);
            
        end
        
        
        
        ind_nonword_all_task = find(datapreprocessed.trialinfo(:,2)== 152 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)== 252 & datapreprocessed.trialinfo(:,4)==1 | ...
            datapreprocessed.trialinfo(:,2)==303 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==306 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==309 & datapreprocessed.trialinfo(:,4)==1);
        
        if length(ind_nonword_all_task) ~= 0
            cfg = [];
            cfg.trials = ind_nonword_all_task;
            erp_nonword_all_task = ft_timelockanalysis(cfg,datapreprocessed);
            
            
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_nonword_all_task = ft_selectdata(cfg,erp_nonword_all_task);
            
            
            
            all_erp_nonword_all_task{1,jj} = vpNames{zz};
            all_erp_nonword_all_task{2,jj} = erp_nonword_all_task;
            all_erp_nonword_all_task{3,jj} = length(ind_nonword_all_task);
            
        end
        
        % T1_T2_correct
        cfg = [];
        
        ind_T1_2_correct = find(datapreprocessed.trialinfo(:,2)== 150 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 151 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 152 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 250 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 251 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 252 & datapreprocessed.trialinfo(:,4)==1);
        
        if length(ind_T1_2_correct) ~= 0
            
            cfg = [];
            cfg.trials = ind_T1_2_correct;
            erp_T1_2_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_2_correct = ft_selectdata(cfg,erp_T1_2_correct);
            
            all_erp_T1_2_correct{1,jj} = vpNames{zz};
            all_erp_T1_2_correct{2,jj} = erp_T1_2_correct;
            all_erp_T1_2_correct{3,jj} = length(ind_T1_2_correct);
            
        end
        
        % T1_correct, T2_correct, T1_word, T1_nonword, T2_word, T2_nonword
        cfg = [];
        
        ind_T1_correct = find(datapreprocessed.trialinfo(:,2)== 150 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 151 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 152 & datapreprocessed.trialinfo(:,4)==1);
        ind_T2_correct = find(datapreprocessed.trialinfo(:,2)== 250 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 251 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 252 & datapreprocessed.trialinfo(:,4)==1);
        
        ind_T1_word_correct = find(datapreprocessed.trialinfo(:,2)== 150 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 151 & datapreprocessed.trialinfo(:,4)==1);
        ind_T1_nonword_correct = find(datapreprocessed.trialinfo(:,2)== 152 & datapreprocessed.trialinfo(:,4)==1);
        ind_T2_word_correct = find(datapreprocessed.trialinfo(:,2)== 250 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)== 251 & datapreprocessed.trialinfo(:,4)==1);
        ind_T2_nonword_correct = find(datapreprocessed.trialinfo(:,2)== 252 & datapreprocessed.trialinfo(:,4)==1);
        
        
        if length(ind_T1_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T1_correct;
            erp_T1_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_correct = ft_selectdata(cfg,erp_T1_correct);
            
            all_erp_T1_correct{1,jj} = vpNames{zz};
            all_erp_T1_correct{2,jj} = erp_T1_correct;
            all_erp_T1_correct{3,jj} = length(ind_T1_correct);
        end
        
        if length(ind_T2_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T2_correct;
            erp_T2_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T2_correct = ft_selectdata(cfg,erp_T2_correct);
            
            all_erp_T2_correct{1,jj} = vpNames{zz};
            all_erp_T2_correct{2,jj} = erp_T2_correct;
            all_erp_T2_correct{3,jj} = length(ind_T2_correct);
        end
        
        if length(ind_T1_word_correct) ~= 0
            
            cfg = [];
            cfg.trials = ind_T1_word_correct;
            erp_T1_word_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_word_correct = ft_selectdata(cfg,erp_T1_word_correct);
            
            all_erp_T1_word{1,jj} = vpNames{zz};
            all_erp_T1_word{2,jj} = erp_T1_word_correct;
            all_erp_T1_word{3,jj} = length(ind_T1_word_correct);
            
        end
        
        if length(ind_T1_nonword_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T1_nonword_correct;
            erp_T1_nonword_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_nonword_correct = ft_selectdata(cfg,erp_T1_nonword_correct);
            
            all_erp_T1_nonword{1,jj} = vpNames{zz};
            all_erp_T1_nonword{2,jj} = erp_T1_nonword_correct;
            all_erp_T1_nonword{3,jj} = length(ind_T1_nonword_correct);
        end
        
        if length(ind_T2_word_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T2_word_correct;
            erp_T2_word_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T2_word_correct = ft_selectdata(cfg,erp_T2_word_correct);
            
            all_erp_T2_word{1,jj} = vpNames{zz};
            all_erp_T2_word{2,jj} = erp_T2_word_correct;
            all_erp_T2_word{3,jj} = length(ind_T2_word_correct);
        end
        
        
        if length(ind_T2_nonword_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T2_nonword_correct;
            erp_T2_nonword_correct = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T2_nonword_correct = ft_selectdata(cfg,erp_T2_nonword_correct);
            
            
            all_erp_T2_nonword{1,jj} = vpNames{zz};
            all_erp_T2_nonword{2,jj} = erp_T2_nonword_correct;
            all_erp_T2_nonword{3,jj} = length(ind_T2_nonword_correct);
        end
        
        % T1_HF, T1_LF, T2_HF, T2_LF
        
        
        cfg = [];
        
        ind_T1_LF_correct = find(datapreprocessed.trialinfo(:,2)==150 & datapreprocessed.trialinfo(:,4)==1);
        ind_T1_HF_correct = find(datapreprocessed.trialinfo(:,2)==151 & datapreprocessed.trialinfo(:,4)==1);
        ind_T2_LF_correct = find(datapreprocessed.trialinfo(:,2)==250 & datapreprocessed.trialinfo(:,4)==1);
        ind_T2_HF_correct = find(datapreprocessed.trialinfo(:,2)==251 & datapreprocessed.trialinfo(:,4)==1);
        
        if length(ind_T1_LF_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T1_LF_correct;
            erp_T1_LF_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_LF = ft_selectdata(cfg,erp_T1_LF_correct);
            
            all_erp_T1_LF{1,jj} = vpNames{zz};
            all_erp_T1_LF{2,jj} = erp_T1_LF;
            all_erp_T1_LF{3,jj} = length(ind_T1_LF_correct);
        end
        
        
        
        if length(ind_T1_HF_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T1_HF_correct;
            erp_T1_HF_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T1_HF = ft_selectdata(cfg,erp_T1_HF_correct);
            
            all_erp_T1_HF{1,jj} = vpNames{zz};
            all_erp_T1_HF{2,jj} = erp_T1_HF;
            all_erp_T1_HF{3,jj} = length(ind_T1_HF_correct);
        end
        
        if length(ind_T2_LF_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T2_LF_correct;
            erp_T2_LF_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T2_LF = ft_selectdata(cfg,erp_T2_LF_correct);
            
            all_erp_T2_LF{1,jj} = vpNames{zz};
            all_erp_T2_LF{2,jj} = erp_T2_LF;
            all_erp_T2_LF{3,jj} = length(ind_T2_LF_correct);
        end
        
        if length(ind_T2_HF_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T2_HF_correct;
            erp_T2_HF_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T2_HF = ft_selectdata(cfg,erp_T2_HF_correct);
            
            all_erp_T2_HF{1,jj} = vpNames{zz};
            all_erp_T2_HF{2,jj} = erp_T2_HF;
            all_erp_T2_HF{3,jj} = length(ind_T2_HF_correct);
        end
        
        
        
        
        
        
        
        
        % Switch task
        
        cfg = [];
        
        ind_T3_word_correct = find(datapreprocessed.trialinfo(:,2)==301 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==302 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==304 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==305 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==307 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==308 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_nonword_correct = find(datapreprocessed.trialinfo(:,2)==303 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==306 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==309 & datapreprocessed.trialinfo(:,4)==1);
        
        ind_T3_G_correct = find(datapreprocessed.trialinfo(:,2)==301 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==304 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==307 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_E_correct = find(datapreprocessed.trialinfo(:,2)==302 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==305 & datapreprocessed.trialinfo(:,4)==1 | datapreprocessed.trialinfo(:,2)==308 & datapreprocessed.trialinfo(:,4)==1);
        
        ind_T3_G_G_correct = find(datapreprocessed.trialinfo(:,2)==301 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_G_E_correct = find(datapreprocessed.trialinfo(:,2)==302 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_G_NW_correct = find(datapreprocessed.trialinfo(:,2)==303 & datapreprocessed.trialinfo(:,4)==1);
        
        ind_T3_E_G_correct = find(datapreprocessed.trialinfo(:,2)==304 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_E_E_correct = find(datapreprocessed.trialinfo(:,2)==305 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_E_NW_correct = find(datapreprocessed.trialinfo(:,2)==306 & datapreprocessed.trialinfo(:,4)==1);
        
        ind_T3_NW_G_correct = find(datapreprocessed.trialinfo(:,2)==307 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_NW_E_correct = find(datapreprocessed.trialinfo(:,2)==308 & datapreprocessed.trialinfo(:,4)==1);
        ind_T3_NW_NW_correct = find(datapreprocessed.trialinfo(:,2)==309 & datapreprocessed.trialinfo(:,4)==1);
        
        
        if length(ind_T3_word_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_word_correct;
            erp_T3_word_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_word = ft_selectdata(cfg,erp_T3_word_correct);
            
            all_erp_T3_word{1,jj} = vpNames{zz};
            all_erp_T3_word{2,jj} = erp_T3_word;
            all_erp_T3_word{3,jj} = length(ind_T3_word_correct);
        end
        
        
        if length(ind_T3_nonword_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_nonword_correct;
            erp_T3_nonword_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_nonword = ft_selectdata(cfg,erp_T3_nonword_correct);
            
            all_erp_T3_nonword{1,jj} = vpNames{zz};
            all_erp_T3_nonword{2,jj} = erp_T3_nonword;
            all_erp_T3_nonword{3,jj} = length(ind_T3_nonword_correct);
        end
        
        
        if length(ind_T3_G_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_G_correct;
            erp_T3_G_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_G = ft_selectdata(cfg,erp_T3_G_correct);
            all_erp_T3_G{1,jj} = vpNames{zz};
            all_erp_T3_G{2,jj} = erp_T3_G;
            all_erp_T3_G{3,jj} = length(ind_T3_G_correct);
        end
        
        if length(ind_T3_E_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_E_correct;
            erp_T3_E_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_E = ft_selectdata(cfg,erp_T3_E_correct);
            
            all_erp_T3_E{1,jj} = vpNames{zz};
            all_erp_T3_E{2,jj} = erp_T3_E;
            all_erp_T3_E{3,jj} = length(ind_T3_E_correct);
        end
        
        if length(ind_T3_G_G_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_G_G_correct;
            erp_T3_G_G_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_G_G = ft_selectdata(cfg,erp_T3_G_G_correct);
            
            all_erp_T3_G_G{1,jj} = vpNames{zz};
            all_erp_T3_G_G{2,jj} = erp_T3_G_G;
            all_erp_T3_G_G{3,jj} = length(ind_T3_G_G_correct);
        end
        
        if length(ind_T3_G_E_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_G_E_correct;
            erp_T3_G_E_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_G_E = ft_selectdata(cfg,erp_T3_G_E_correct);
            
            all_erp_T3_G_E{1,jj} = vpNames{zz};
            all_erp_T3_G_E{2,jj} = erp_T3_G_E;
            all_erp_T3_G_E{3,jj} = length(ind_T3_G_E_correct);
        end
        
        if length(ind_T3_G_NW_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_G_NW_correct;
            erp_T3_G_NW_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_G_NW = ft_selectdata(cfg,erp_T3_G_NW_correct);
            
            all_erp_T3_G_NW{1,jj} = vpNames{zz};
            all_erp_T3_G_NW{2,jj} = erp_T3_G_NW;
            all_erp_T3_G_NW{3,jj} = length(ind_T3_G_NW_correct);
        end
        
        
        if length(ind_T3_E_G_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_E_G_correct;
            erp_T3_E_G_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_E_G = ft_selectdata(cfg,erp_T3_E_G_correct);
            
            all_erp_T3_E_G{1,jj} = vpNames{zz};
            all_erp_T3_E_G{2,jj} = erp_T3_E_G;
            all_erp_T3_E_G{3,jj} = length(ind_T3_E_G_correct);
        end
        
        if length(ind_T3_E_E_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_E_E_correct;
            erp_T3_E_E_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_E_E = ft_selectdata(cfg,erp_T3_E_E_correct);
            
            
            all_erp_T3_E_E{1,jj} = vpNames{zz};
            all_erp_T3_E_E{2,jj} = erp_T3_E_E;
            all_erp_T3_E_E{3,jj} = length(ind_T3_E_E_correct);
        end
        
        if length(ind_T3_E_NW_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_E_NW_correct;
            erp_T3_E_NW_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_E_NW = ft_selectdata(cfg,erp_T3_E_NW_correct);
            
            all_erp_T3_E_NW{1,jj} = vpNames{zz};
            all_erp_T3_E_NW{2,jj} = erp_T3_E_NW;
            all_erp_T3_E_NW{3,jj} = length(ind_T3_E_NW_correct);
        end
        
        
        if length(ind_T3_NW_G_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_NW_G_correct;
            erp_T3_NW_G_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_NW_G = ft_selectdata(cfg,erp_T3_NW_G_correct);
            
            all_erp_T3_NW_G{1,jj} = vpNames{zz};
            all_erp_T3_NW_G{2,jj} = erp_T3_NW_G;
            all_erp_T3_NW_G{3,jj} = length(ind_T3_NW_G_correct);
        end
        
        if length(ind_T3_NW_E_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_NW_E_correct;
            erp_T3_NW_E_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_NW_E = ft_selectdata(cfg,erp_T3_NW_E_correct);
            
            all_erp_T3_NW_E{1,jj} = vpNames{zz};
            all_erp_T3_NW_E{2,jj} = erp_T3_NW_E;
            all_erp_T3_NW_E{3,jj} = length(ind_T3_NW_E_correct);
        end
        
        
        
        if length(ind_T3_NW_NW_correct) ~= 0
            cfg = [];
            cfg.trials = ind_T3_NW_NW_correct;
            erp_T3_NW_NW_correct   = ft_timelockanalysis(cfg,datapreprocessed);
            
            cfg = [];
            cfg.latency = [-.1 1];
            erp_T3_NW_NW = ft_selectdata(cfg,erp_T3_NW_NW_correct);
            
            
            all_erp_T3_NW_NW{1,jj} = vpNames{zz};
            all_erp_T3_NW_NW{2,jj} = erp_T3_NW_NW;
            all_erp_T3_NW_NW{3,jj} = length(ind_T3_NW_NW_correct);
        end
        
        jj = jj+1;
        
    end
end


save([dataPath 'LDT_erp_visrej_ndm_nbl_18.mat'], ...
    'all_erp_all_correct', 'all_erp_word_all_task', 'all_erp_nonword_all_task', ...
    'all_erp_T1_2_correct','all_erp_T1_correct','all_erp_T2_correct',...
    ...
    'all_erp_T1_word', 'all_erp_T1_nonword', 'all_erp_T2_word' , 'all_erp_T2_nonword',...
    ...
    'all_erp_T1_HF', 'all_erp_T1_LF',...
    'all_erp_T2_HF', 'all_erp_T2_LF',...
    ...
    'all_erp_T3_word','all_erp_T3_nonword','all_erp_T3_G','all_erp_T3_E',...
    'all_erp_T3_G_G', 'all_erp_T3_G_E', 'all_erp_T3_G_NW',...
    'all_erp_T3_E_G', 'all_erp_T3_E_E', 'all_erp_T3_E_NW',...
    'all_erp_T3_NW_G', 'all_erp_T3_NW_E', 'all_erp_T3_NW_NW')

cd(fileparts(tmp.Filename));






