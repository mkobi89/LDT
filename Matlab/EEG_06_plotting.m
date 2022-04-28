% exract mean amplitudes with fieldtrip

%% Preparation:
clear

%% Change directory to the current m. file
% by hand or by code:
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Path to server depending on OS

% if IsOSX == true
%     dataPath = '/Volumes/CLINT/LDT_preprocessed/'; 
%     addpath('/Volumes/CLINT/fieldtrip-20200607/'); % add fieldtrip path
% else
%     dataPath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT_final_preprocessed/'; % LDT Path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/fieldtrip-20200607/'); % add fieldtrip path
% end

dataPath = 'F:\LDT_final_2_preprocessed/'; % All Data Path
addpath('F:\CLINT backup_15.02.2022\fieldtrip-20200607/'); % add fieldtrip path

%% Starting fieldtrip

ft_defaults


%% load files for headmodeling

load('headmodeling/eeglabchans');
load('headmodeling/elec_aligned');
load('headmodeling/labels105');
load('headmodeling/lay129_head');

%% Layout plot
figure
ft_plot_layout(lay129_head, 'box', 'no' )
figure
ft_plot_layout(lay129_head, 'box', 'no', 'chanindx', find(ismember(lay129_head.label,labels105)) )


%% load data ROI vis rej version
cd(dataPath);
load('LDT_erp_visrej_ndm_nbl_18.mat');
vpNames = all_erp_T1_HF(1,:);

%% grand averages over all tasks and trials

% remove unwanted subjects
ind_subjects = logical(zeros(1,length(vpNames)));
ind_subjects(1:end) = true; 
deselect = [20,31,42,49,57,72,75,85];

for i = 1:length(deselect)
    ind_subjects(deselect(i)) = false;
end


% Prepare erp for GA all trials and words vs nonwords

% all trials
all_erp_all_correct_reduced = all_erp_all_correct(:,ind_subjects);
cfg = [];
ga_erp_all_tasks = ft_timelockgrandaverage(cfg, all_erp_all_correct_reduced{2,:});

% word vs nonword all tasks
all_erp_word_all_task_reduced = all_erp_word_all_task(:,ind_subjects);
all_erp_nonword_all_task_reduced = all_erp_nonword_all_task(:,ind_subjects);
cfg = [];
ga_erp_word_all_task = ft_timelockgrandaverage(cfg, all_erp_word_all_task_reduced{2,:});
ga_erp_nonword_all_task = ft_timelockgrandaverage(cfg, all_erp_nonword_all_task_reduced{2,:});


lay129_head.height = lay129_head.height*2;
lay129_head.width = lay129_head.width*2;
%subjects = all_erp_T1_2_correct_reduced(1,:);

%% Select Cz + 5 electrodes around: E7, E106, E80, E55, E31, ECz
% defining for highlighting in the topoplot
chans = {'E7', 'E106', 'E80', 'E55', 'E31', 'Cz'};

%% plot topoplots to select electrodes, grand average across all conditions
cfg = [];
cfg.layout = lay129_head;
% cfg.height = 1;
% cfg.width = 1;
cfg.channel = {'all'};
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq   = 35;
cfg.highlight        = 'on';
cfg.highlightsymbol  = '*';
cfg.highlightchannel = chans;
cfg.higlightfontsize = 15;
cfg.showlabels = 'yes';
cfg.xlim = 0.2:0.05:0.5;
cfg.zlim = [-1 1];
% plot every 50ms
cfg.comment = 'xlim';
cfg.commentpos = 'leftbottom';
% cfg.colorbar = 'yes';
% cfg.colorbartext     = 'µV';
% cfg.gridscale = 500;
figure; ft_topoplotTFR(cfg, ga_erp_all_tasks);


%% ERP GA plots for condition
all_erp_T1_HF_reduced = all_erp_T1_HF(:,ind_subjects);
all_erp_T1_LF_reduced = all_erp_T1_LF(:,ind_subjects);
all_erp_T1_nonword_reduced = all_erp_T1_nonword(:,ind_subjects);

all_erp_T2_HF_reduced = all_erp_T2_HF(:,ind_subjects);
all_erp_T2_LF_reduced = all_erp_T2_LF(:,ind_subjects);
all_erp_T2_nonword_reduced = all_erp_T2_nonword(:,ind_subjects);

all_erp_T3_G_G_reduced = all_erp_T3_G_G(:,ind_subjects);
all_erp_T3_G_E_reduced = all_erp_T3_G_E(:,ind_subjects);
all_erp_T3_E_G_reduced = all_erp_T3_E_G(:,ind_subjects);
all_erp_T3_E_E_reduced = all_erp_T3_E_E(:,ind_subjects);

clearvars all_erp_T1_HF all_erp_T1_LF all_erp_T1_nonword all_erp_T2_HF all_erp_T2_LF all_erp_T2_nonword all_erp_T3_G_G all_erp_T3_G_E all_erp_T3_E_G all_erp_T3_E_E

cfg = [];
ga_erp_T1_HF_correct = ft_timelockgrandaverage(cfg, all_erp_T1_HF_reduced{2,:});
ga_erp_T1_LF_correct = ft_timelockgrandaverage(cfg, all_erp_T1_LF_reduced{2,:});
ga_erp_T1_nonword = ft_timelockgrandaverage(cfg, all_erp_T1_nonword_reduced{2,:});

ga_erp_T2_HF_correct = ft_timelockgrandaverage(cfg, all_erp_T2_HF_reduced{2,:});
ga_erp_T2_LF_correct = ft_timelockgrandaverage(cfg, all_erp_T2_LF_reduced{2,:});
ga_erp_T2_nonword = ft_timelockgrandaverage(cfg, all_erp_T2_nonword_reduced{2,:});

ga_erp_T3_G_G_correct = ft_timelockgrandaverage(cfg, all_erp_T3_G_G_reduced{2,:});
ga_erp_T3_G_E_correct = ft_timelockgrandaverage(cfg, all_erp_T3_G_E_reduced{2,:});
ga_erp_T3_E_G_correct = ft_timelockgrandaverage(cfg, all_erp_T3_E_G_reduced{2,:});
ga_erp_T3_E_E_correct = ft_timelockgrandaverage(cfg, all_erp_T3_E_E_reduced{2,:});

clearvars all_erp_T1_HF_reduced all_erp_T1_LF_reduced all_erp_T1_nonword_reduced all_erp_T2_HF_reduced all_erp_T2_LF_reduced all_erp_T2_nonword_reduced all_erp_T3_G_G_reduced all_erp_T3_G_E_reduced all_erp_T3_E_G_reduced all_erp_T3_E_E_reduced

% German high-, low-frequency and pseudowords
cfg = [];
cfg.layout = lay129_head;
cfg.channel = chans;
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq   = 35;
%cfg.showlabels = 'yes';
cfg.xlim = [-0.1 0.6];
cfg.ylim = [-1.2 1.4];
figure; ft_singleplotER(cfg, ga_erp_T1_HF_correct, ga_erp_T1_LF_correct, ga_erp_T1_nonword)
legend('German HF','German LF', 'German PW');

% English high-, low-frequency and pseudowords
cfg = [];
cfg.layout = lay129_head;
cfg.channel = chans;
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq   = 35;
cfg.showlabels = 'yes';
cfg.xlim = [-0.1 0.6];
cfg.ylim = [-1.2 1.4];
figure; ft_singleplotER(cfg, ga_erp_T2_HF_correct, ga_erp_T2_LF_correct, ga_erp_T2_nonword)
legend('English HF','English LF', 'English PW');

% Switch conditions
cfg = [];
cfg.layout = lay129_head;
cfg.channel = chans;
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq   = 35;
cfg.showlabels = 'yes';
cfg.xlim = [-0.1 0.6];
cfg.ylim = [-1.2 1.4];
figure; ft_singleplotER(cfg, ga_erp_T3_G_G_correct, ga_erp_T3_E_G_correct, ga_erp_T3_E_E_correct,ga_erp_T3_G_E_correct)
legend('GG nonswich','EG switch', 'EE nonswitch', 'GE switch');



%% topoplots with selected channels highlighted for words and nonwords
cfg = [];
cfg.layout = lay129_head;
cfg.xlim = [0.25 0.45];
% plot every 50ms
cfg.highlight        = 'on';
cfg.highlightsymbol  = '*';
cfg.highlightchannel = chans;
cfg.higlightfontsize = 15;
cfg.colorbar         = 'eastoutside';
cfg.colorbartext     = 'µV';
cfg.showlabels       = 'yes';
cfg.zlim             =  [-1 1];
%cfg.title('GA word all task ROI vis rej 0.25 - 0.45');
% figure; ft_topoplotTFR(cfg, ga_erp_all_tasks);
% legend('GA all cond all task ROI vis rej 0.25 - 0.45');
figure; ft_topoplotTFR(cfg, ga_erp_word_all_task);

figure; ft_topoplotTFR(cfg, ga_erp_nonword_all_task);
%legend('GA nonword all task ROI vis rej 0.25 - 0.45');


