%% exract mean amplitudes with fieldtrip

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
%     liesefeldPath = ('/Users/ninaraduner/Documents/Studium/A_Studium/Masterarbeit/LDT/Matlab/latency-master/');
% else
%     dataPath = '\\130.60.235.123\Users\neuro\Desktop\CLINT\LDT_preprocessed/'; % LDT Path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/fieldtrip-20200607/'); % add fieldtrip path
%     liesefeldPath = ('C:\Users\matth\Documents\LDT\Matlab\latency-master\'); 
% end


dataPath = 'F:\LDT_final_2_preprocessed/'; % All Data Path
addpath('F:\CLINT backup_15.02.2022\fieldtrip-20200607/'); % add fieldtrip path
liesefeldPath = ('C:\Users\matth\Documents\LDT\Matlab\latency-master\');

%% Starting fieldtrip

ft_defaults

%% load files for headmodeling

load('headmodeling/eeglabchans');
load('headmodeling/elec_aligned');
load('headmodeling/labels105');
load('headmodeling/lay129_head');


%% load data
cd(dataPath);
load('LDT_erp_visrej_ndm_nbl_18.mat');
vpNames = all_erp_T1_HF(1,:);


%% Liesefeld approach
cd(liesefeldPath)

% select channels
% Cz + 5 Elektroden drum herum: E7, E106, E80, E55, E31, ECz
chans = {'E7', 'E106', 'E80', 'E55', 'E31', 'Cz'};

% bring data into right format
nSubs  = length(all_erp_T1_HF);
nElecs = size(all_erp_T1_HF{2,1}.avg,1);
nTime  = size(all_erp_T1_HF{2,1}.avg,2);

T1_G_HF = nan(nSubs,nElecs,nTime);
T1_G_LF = nan(nSubs,nElecs,nTime);
T1_G_NW = nan(nSubs,nElecs,nTime);

T2_E_HF = nan(nSubs,nElecs,nTime);
T2_E_LF = nan(nSubs,nElecs,nTime);
T2_E_NW = nan(nSubs,nElecs,nTime);

T3_G_G = nan(nSubs,nElecs,nTime);
T3_G_E = nan(nSubs,nElecs,nTime);
T3_G_NW = nan(nSubs,nElecs,nTime);

T3_E_G = nan(nSubs,nElecs,nTime);
T3_E_E = nan(nSubs,nElecs,nTime);
T3_E_NW = nan(nSubs,nElecs,nTime);

T3_NW_G = nan(nSubs,nElecs,nTime);
T3_NW_E = nan(nSubs,nElecs,nTime);
T3_NW_NW = nan(nSubs,nElecs,nTime);

for subi= 1:nSubs
    T1_G_HF(subi,:,:)  = all_erp_T1_HF{2,subi}.avg;
    T1_G_LF(subi,:,:)  = all_erp_T1_LF{2,subi}.avg;
    T1_G_NW(subi,:,:)  = all_erp_T1_nonword{2,subi}.avg;

    T2_E_HF(subi,:,:)  = all_erp_T2_HF{2,subi}.avg;
    T2_E_LF(subi,:,:)  = all_erp_T2_LF{2,subi}.avg;
    T2_E_NW(subi,:,:)  = all_erp_T2_nonword{2,subi}.avg;

    T3_G_G(subi,:,:)   = all_erp_T3_G_G{2,subi}.avg;
    T3_G_E(subi,:,:)   = all_erp_T3_G_E{2,subi}.avg;
    T3_G_NW(subi,:,:)  = all_erp_T3_G_NW{2,subi}.avg;
    
    T3_E_G(subi,:,:)   = all_erp_T3_E_G{2,subi}.avg;
    T3_E_E(subi,:,:)   = all_erp_T3_E_E{2,subi}.avg;
    T3_E_NW(subi,:,:)  = all_erp_T3_E_NW{2,subi}.avg;
    
    T3_NW_G(subi,:,:)  = all_erp_T3_NW_G{2,subi}.avg;
    T3_NW_E(subi,:,:)  = all_erp_T3_NW_E{2,subi}.avg;
    T3_NW_NW(subi,:,:) = all_erp_T3_NW_NW{2,subi}.avg;
end

meanERF_T1=(T1_G_HF + T1_G_LF + T1_G_NW) / 3;
meanERF_T2=(T2_E_HF + T2_E_LF + T2_E_NW) / 3;
meanERF_T3=(T3_G_G + T3_G_E + T3_G_NW + T3_E_G + T3_E_E + T3_E_NW + T3_NW_G + T3_NW_E + T3_NW_NW) / 9;

% baseline correction
%the baseline currently is at 0.5 x 10e-12, because of the combination of
%two gradiometers by root-mean-square; this is removed, because it would
%confound various latency measures and makes interpretation of some
%paramters more cumbersome
T1_G_HF = T1_G_HF - repmat(mean(T1_G_HF(:,:,1:51),3), [1,1,nTime]);
T1_G_LF = T1_G_LF - repmat(mean(T1_G_LF(:,:,1:51),3), [1,1,nTime]);
T1_G_NW = T1_G_NW - repmat(mean(T1_G_NW(:,:,1:51),3), [1,1,nTime]);

T2_E_HF = T2_E_HF - repmat(mean(T2_E_HF(:,:,1:51),3), [1,1,nTime]);
T2_E_LF = T2_E_LF - repmat(mean(T2_E_LF(:,:,1:51),3), [1,1,nTime]);
T2_E_NW = T2_E_NW - repmat(mean(T2_E_NW(:,:,1:51),3), [1,1,nTime]);

T3_G_G   = T3_G_G  - repmat(mean(T3_G_G(:,:,1:51),3), [1,1,nTime]);
T3_G_E   = T3_G_E  - repmat(mean(T3_G_E(:,:,1:51),3), [1,1,nTime]);
T3_G_NW  = T3_G_NW - repmat(mean(T3_G_NW(:,:,1:51),3), [1,1,nTime]);

T3_E_G   = T3_E_G  - repmat(mean(T3_E_G(:,:,1:51),3), [1,1,nTime]);
T3_E_E   = T3_E_E  - repmat(mean(T3_E_E(:,:,1:51),3), [1,1,nTime]);
T3_E_NW  = T3_E_NW - repmat(mean(T3_E_NW(:,:,1:51),3), [1,1,nTime]);

T3_NW_G  = T3_NW_G  - repmat(mean(T3_NW_G(:,:,1:51),3), [1,1,nTime]);
T3_NW_E  = T3_NW_E  - repmat(mean(T3_NW_E(:,:,1:51),3), [1,1,nTime]);
T3_NW_NW = T3_NW_NW - repmat(mean(T3_NW_NW(:,:,1:51),3), [1,1,nTime]);

%extract information on time and channels
time      = all_erp_T1_HF{2,1}.time;
chanName  = all_erp_T1_HF{2,1}.label;
ind_chans = find(ismember(chanName,chans));

%% Liesefeld with parameters for our data
cfg = [];
cfg.chanName  = chanName;       %tell the function about the channel names
cfg.times     = time;           %tell the function about the time axis
cfg.sign      = -1;             %search for a negative component
cfg.chans     = chans;          %set channel to analyze
cfg.areaWin   = 'ampLat';       %area is confined by the on- and offsets
cfg.percAmp   = 0.3;            %percentage of amplitude for on- and offsets
cfg.percArea  = 0.5;            %defaults to 50% anyway
cfg.peakWin   = [0.25 0.45];      %set window for searching the peak
%cfg.cWinWidth = -0.15;           %set window for searching the adjacent peak (relative to cfg.peakWin)
%cfg.cWinStart = 0.1;            %where to start search for previouse peak
cfg.peakWidth = 0.010;              %set window for percent area latency
cfg.fig       = true;           %request figures for individual averages
cfg.areaBase  = 'zero';      %individualize baseline

cfg.extract   = {'peakLat', 'onset', 'offset', 'areaLat', 'mean', 'peakAmp', 'area', 'width', 'baseline'};        %'peak2peak' left, need to define cWinWidth or cWinStart

% for i = 1:83
 cfg.subs  = [47];
latency_T1_G_HF = latency(cfg,T1_G_HF);
latency_T1_G_LF = latency(cfg,T1_G_LF);
latency_T1_G_NW = latency(cfg,T1_G_NW);


latency_T2_E_HF = latency(cfg,T2_E_HF);
latency_T2_E_LF = latency(cfg,T2_E_LF);
latency_T2_E_NW = latency(cfg,T2_E_NW);

latency_T3_G_G  = latency(cfg,T3_G_G);
latency_T3_E_E  = latency(cfg,T3_E_E);
latency_T3_E_G  = latency(cfg,T3_E_G);
latency_T3_G_E  = latency(cfg,T3_G_E);

% pause(25)
% close all
% end

latency_T3_G_NW = latency(cfg,T3_G_NW);
latency_T3_E_NW = latency(cfg,T3_E_NW);
latency_T3_NW_G  = latency(cfg,T3_NW_G);
latency_T3_NW_E  = latency(cfg,T3_NW_E);
latency_T3_NW_NW = latency(cfg,T3_NW_NW);

% extract onset and offset in time points
cfg = [];
cfg.chanName  = chanName;       %tell the function about the channel names
%cfg.times     = time;           %tell the function about the time axis
cfg.sign      = -1;             %search for a negative component
cfg.chans     = chans;          %set channel to analyze
cfg.areaWin   = 'ampLat';       %area is confined by the on- and offsets
cfg.percAmp   = 0.3;            %percentage of amplitude for on- and offsets
cfg.percArea  = 0.5;            %defaults to 50% anyway
cfg.peakWin   = [176 276];      %set window for searching the peak
%cfg.cWinWidth = -0.15;           %set window for searching the adjacent peak (relative to cfg.peakWin)
%cfg.cWinStart = 0.1;            %where to start search for previouse peak
cfg.peakWidth = 5;              %set window for percent area latency
cfg.fig       = false;           %request figures for individual averages
cfg.areaBase  = 'zero';      %individualize baseline

cfg.extract   = {'onset', 'offset'};%'peak2peak' left, need to define cWinWidth or cWinStart

%cfg.subs = 61;
on_off_T1_G_HF = latency(cfg,T1_G_HF);
on_off_T1_G_LF = latency(cfg,T1_G_LF);
on_off_T1_G_NW = latency(cfg,T1_G_NW);

on_off_T2_E_HF = latency(cfg,T2_E_HF);
on_off_T2_E_LF = latency(cfg,T2_E_LF);
on_off_T2_E_NW = latency(cfg,T2_E_NW);

on_off_T3_G_G = latency(cfg,T3_G_G);
on_off_T3_G_E = latency(cfg,T3_G_E);
on_off_T3_G_NW = latency(cfg,T3_G_NW);

on_off_T3_E_G = latency(cfg,T3_E_G);
on_off_T3_E_E = latency(cfg,T3_E_E);
on_off_T3_E_NW = latency(cfg,T3_E_NW);

on_off_T3_NW_G = latency(cfg,T3_NW_G);
on_off_T3_NW_E = latency(cfg,T3_NW_E);
on_off_T3_NW_NW = latency(cfg,T3_NW_NW);


%% find cases, where there is no N400
ind_special_T1_HF = find(isnan(latency_T1_G_HF.areaLat) & latency_T1_G_HF.peakAmp >= 0);
ind_special_T1_LF = find(isnan(latency_T1_G_LF.areaLat) & latency_T1_G_LF.peakAmp >= 0);
ind_special_T1_NW = find(isnan(latency_T1_G_NW.areaLat) & latency_T1_G_NW.peakAmp >= 0);

ind_special_T2_HF = find(isnan(latency_T2_E_HF.areaLat) & latency_T2_E_HF.peakAmp >= 0);
ind_special_T2_LF = find(isnan(latency_T2_E_LF.areaLat) & latency_T2_E_LF.peakAmp >= 0);
ind_special_T2_NW = find(isnan(latency_T2_E_NW.areaLat) & latency_T2_E_NW.peakAmp >= 0);

ind_special_T3_G_G = find(isnan(latency_T3_G_G.areaLat) & latency_T3_G_G.peakAmp >= 0);
ind_special_T3_G_E = find(isnan(latency_T3_G_E.areaLat) & latency_T3_G_E.peakAmp >= 0);
% ind_special_T3_G_NW = find(isnan(latency_T3_G_NW.areaLat) & latency_T3_G_NW.peakAmp >= 0);

ind_special_T3_E_G = find(isnan(latency_T3_E_G.areaLat) & latency_T3_E_G.peakAmp >= 0);
ind_special_T3_E_E = find(isnan(latency_T3_E_E.areaLat) & latency_T3_E_E.peakAmp >= 0);
% ind_special_T3_E_NW = find(isnan(latency_T3_E_NW.areaLat) & latency_T3_E_NW.peakAmp >= 0);
% 
% ind_special_T3_NW_G = find(isnan(latency_T3_NW_G.areaLat) & latency_T3_NW_G.peakAmp >= 0);
% ind_special_T3_NW_E = find(isnan(latency_T3_NW_E.areaLat) & latency_T3_NW_E.peakAmp >= 0);
% ind_special_T3_NW_NW = find(isnan(latency_T3_NW_NW.areaLat) & latency_T3_NW_NW.peakAmp >= 0);

ind_special_T1 = union(ind_special_T1_HF, ind_special_T1_LF);
ind_special_T1 = union(ind_special_T1, ind_special_T1_NW);

ind_special_T2 = union(ind_special_T2_HF, ind_special_T2_LF);
ind_special_T2 = union(ind_special_T2, ind_special_T2_NW);

ind_special_T1_T2 = union(ind_special_T1, ind_special_T2);

ind_special_T3 = union(ind_special_T3_G_G, ind_special_T3_G_E);
% ind_special_T3 = union(ind_special_T3, ind_special_T3_G_NW);
ind_special_T3 = union(ind_special_T3, ind_special_T3_E_G);
ind_special_T3 = union(ind_special_T3, ind_special_T3_E_E);
% ind_special_T3 = union(ind_special_T3, ind_special_T3_E_NW);
% ind_special_T3 = union(ind_special_T3, ind_special_T3_NW_G);
% ind_special_T3 = union(ind_special_T3, ind_special_T3_NW_E);
% ind_special_T3 = union(ind_special_T3, ind_special_T3_NW_NW);

ind_special_all = union(ind_special_T1_T2, ind_special_T3);

%% Task 1
% column 1: HF; column 2: LF; column 3: NW

onsets_special_T1 = nan(length(ind_special_T1), 3);
offsets_special_T1 = nan(length(ind_special_T1), 3);

keep_onsets_special_T1 = ones(length(ind_special_T1), 3);
keep_offsets_special_T1 = ones(length(ind_special_T1), 3);

for i = 1:length(ind_special_T1)
    % create array for onsets
    onsets_special_T1(i,1) = on_off_T1_G_HF.onset(ind_special_T1(i));
    onsets_special_T1(i,2) = on_off_T1_G_LF.onset(ind_special_T1(i));
    onsets_special_T1(i,3) = on_off_T1_G_NW.onset(ind_special_T1(i));
  
    % create array for offsets
    offsets_special_T1(i,1) = on_off_T1_G_HF.offset(ind_special_T1(i));
    offsets_special_T1(i,2) = on_off_T1_G_LF.offset(ind_special_T1(i));
    offsets_special_T1(i,3) = on_off_T1_G_NW.offset(ind_special_T1(i));
    
    if ismember(ind_special_T1(i), ind_special_T1_HF)
        keep_onsets_special_T1(i,1) = 0;
        keep_offsets_special_T1(i,1) = 0;
    end
    if ismember(ind_special_T1(i), ind_special_T1_LF)
        keep_onsets_special_T1(i,2) = 0;
        keep_offsets_special_T1(i,2) = 0;
    end
    if ismember(ind_special_T1(i), ind_special_T1_NW)
        keep_onsets_special_T1(i,3) = 0;
        keep_offsets_special_T1(i,3) = 0;
    end
end

means_onset = floor(sum([onsets_special_T1.*keep_onsets_special_T1], 2)./sum(keep_onsets_special_T1,2));
means_offset = floor(sum([offsets_special_T1.*keep_offsets_special_T1], 2)./sum(keep_offsets_special_T1,2));

% s = tps/500-0.102
% tps = 500*s + 51
means_onset_s = means_onset/500 - 0.102;
means_offset_s = means_offset/500 - 0.102;

% replace onsets and offsets
on_off_T1_G_HF_new = on_off_T1_G_HF;
on_off_T1_G_LF_new = on_off_T1_G_LF;
on_off_T1_G_NW_new = on_off_T1_G_NW;

for i = 1:length(ind_special_T1)    
    if means_onset(i) > 0;
    if ismember(ind_special_T1(i), ind_special_T1_HF)
        %on_off_new with new onset for special cases
        on_off_T1_G_HF_new.onset(ind_special_T1(i)) =  means_onset(i,1);
        on_off_T1_G_HF_new.offset(ind_special_T1(i)) =  means_offset(i,1);
        
        % exchange onset/offset in latency struct for each special case
        latency_T1_G_HF.onset(ind_special_T1(i)) =  means_onset_s(i,1);
        latency_T1_G_HF.offset(ind_special_T1(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T1_G_HF, 2));
        erp_mean = erp_mean(ind_special_T1(i),on_off_T1_G_HF_new.onset(ind_special_T1(i)):on_off_T1_G_HF_new.offset(ind_special_T1(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T1_G_HF.areaLat(ind_special_T1(i)) = mean_latency_s;    
    end
    
    if ismember(ind_special_T1(i), ind_special_T1_LF)
        on_off_T1_G_LF_new.onset(ind_special_T1(i)) =  means_onset(i,1);
        on_off_T1_G_LF_new.offset(ind_special_T1(i)) =  means_offset(i,1);
        
        latency_T1_G_LF.onset(ind_special_T1(i)) =  means_onset_s(i,1);
        latency_T1_G_LF.offset(ind_special_T1(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T1_G_LF, 2));
        erp_mean = erp_mean(ind_special_T1(i),on_off_T1_G_LF_new.onset(ind_special_T1(i)):on_off_T1_G_LF_new.offset(ind_special_T1(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T1_G_LF.areaLat(ind_special_T1(i)) = mean_latency_s;
    end 
    
    if ismember(ind_special_T1(i), ind_special_T1_NW)
        on_off_T1_G_NW_new.onset(ind_special_T1(i)) =  means_onset(i,1);
        on_off_T1_G_NW_new.offset(ind_special_T1(i)) =  means_offset(i,1);
        
        latency_T1_G_NW.onset(ind_special_T1(i)) =  means_onset_s(i,1);
        latency_T1_G_NW.offset(ind_special_T1(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T1_G_NW, 2));
        erp_mean = erp_mean(ind_special_T1(i),on_off_T1_G_NW_new.onset(ind_special_T1(i)):on_off_T1_G_NW_new.offset(ind_special_T1(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T1_G_NW.areaLat(ind_special_T1(i)) = mean_latency_s;
    end
    end
end

%% Task 2
% column 1: HF; column 2: LF; column 3: NW

onsets_special_T2 = nan(length(ind_special_T2), 3);
offsets_special_T2 = nan(length(ind_special_T2), 3);

keep_onsets_special_T2 = ones(length(ind_special_T2), 3);
keep_offsets_special_T2 = ones(length(ind_special_T2), 3);

for i = 1:length(ind_special_T2)
    % create array for onsets
    onsets_special_T2(i,1) = on_off_T2_E_HF.onset(ind_special_T2(i));
    onsets_special_T2(i,2) = on_off_T2_E_LF.onset(ind_special_T2(i));
    onsets_special_T2(i,3) = on_off_T2_E_NW.onset(ind_special_T2(i));
  
    % create array for offsets
    offsets_special_T2(i,1) = on_off_T2_E_HF.offset(ind_special_T2(i));
    offsets_special_T2(i,2) = on_off_T2_E_LF.offset(ind_special_T2(i));
    offsets_special_T2(i,3) = on_off_T2_E_NW.offset(ind_special_T2(i));
    
    if ismember(ind_special_T2(i), ind_special_T2_HF)
        keep_onsets_special_T2(i,1) = 0;
        keep_offsets_special_T2(i,1) = 0;
    end
    if ismember(ind_special_T2(i), ind_special_T2_LF)
        keep_onsets_special_T2(i,2) = 0;
        keep_offsets_special_T2(i,2) = 0;
    end
    if ismember(ind_special_T2(i), ind_special_T2_NW)
        keep_onsets_special_T2(i,3) = 0;
        keep_offsets_special_T2(i,3) = 0;
    end
end

means_onset = floor(sum([onsets_special_T2.*keep_onsets_special_T2], 2)./sum(keep_onsets_special_T2,2));
means_offset = floor(sum([offsets_special_T2.*keep_offsets_special_T2], 2)./sum(keep_offsets_special_T2,2));

% s = tps/500-0.102
% tps = 500*s + 51
means_onset_s = means_onset/500 - 0.102;
means_offset_s = means_offset/500 - 0.102;

% replace onsets and offsets
on_off_T2_E_HF_new = on_off_T2_E_HF;
on_off_T2_E_LF_new = on_off_T2_E_LF;
on_off_T2_E_NW_new = on_off_T2_E_NW;

for i = 1:length(ind_special_T2)
        if means_onset(i) > 0;
    if ismember(ind_special_T2(i), ind_special_T2_HF)
        %on_off_new with new onset for special cases
        on_off_T2_E_HF_new.onset(ind_special_T2(i)) =  means_onset(i,1);
        on_off_T2_E_HF_new.offset(ind_special_T2(i)) =  means_offset(i,1);
        
        % exchange onset/offset in latency struct for each special case
        latency_T2_E_HF.onset(ind_special_T2(i)) =  means_onset_s(i,1);
        latency_T2_E_HF.offset(ind_special_T2(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T2_E_HF, 2));
        erp_mean = erp_mean(ind_special_T2(i),on_off_T2_E_HF_new.onset(ind_special_T2(i)):on_off_T2_E_HF_new.offset(ind_special_T2(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T2_E_HF.areaLat(ind_special_T2(i)) = mean_latency_s;    
    end
    
    if ismember(ind_special_T2(i), ind_special_T2_LF)
        on_off_T2_E_LF_new.onset(ind_special_T2(i)) =  means_onset(i,1);
        on_off_T2_E_LF_new.offset(ind_special_T2(i)) =  means_offset(i,1);
        
        latency_T2_E_LF.onset(ind_special_T2(i)) =  means_onset_s(i,1);
        latency_T2_E_LF.offset(ind_special_T2(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T2_E_LF, 2));
        erp_mean = erp_mean(ind_special_T2(i),on_off_T2_E_LF_new.onset(ind_special_T2(i)):on_off_T2_E_LF_new.offset(ind_special_T2(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T2_E_LF.areaLat(ind_special_T2(i)) = mean_latency_s;
    end 
    
    if ismember(ind_special_T2(i), ind_special_T2_NW)
        on_off_T2_E_NW_new.onset(ind_special_T2(i)) =  means_onset(i,1);
        on_off_T2_E_NW_new.offset(ind_special_T2(i)) =  means_offset(i,1);
        
        latency_T2_E_NW.onset(ind_special_T2(i)) =  means_onset_s(i,1);
        latency_T2_E_NW.offset(ind_special_T2(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T2_E_NW, 2));
        erp_mean = erp_mean(ind_special_T2(i),on_off_T2_E_NW_new.onset(ind_special_T2(i)):on_off_T2_E_NW_new.offset(ind_special_T2(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T2_E_NW.areaLat(ind_special_T2(i)) = mean_latency_s;
    end
        end
end

% %% task 3 with 9 conditions
% % column 1: G_G; column 2: G_E; column 3: G_NW
% % column 4: E_G; column 5: E_E; column 6: E_NW
% % column 7: NW_G; column 8: NW_E; column 9: NW_NW
% 
% onsets_special_T3 = nan(length(ind_special_T3), 9);
% offsets_special_T3 = nan(length(ind_special_T3), 9);
% 
% keep_onsets_special_T3 = ones(length(ind_special_T3), 9);
% keep_offsets_special_T3 = ones(length(ind_special_T3), 9);
% 
% for i = 1:length(ind_special_T3)
%     % create array for onsets
%     onsets_special_T3(i,1) = on_off_T3_G_G.onset(ind_special_T3(i));
%     onsets_special_T3(i,2) = on_off_T3_G_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,3) = on_off_T3_G_NW.onset(ind_special_T3(i));
%     
%     onsets_special_T3(i,4) = on_off_T3_E_G.onset(ind_special_T3(i));
%     onsets_special_T3(i,5) = on_off_T3_E_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,6) = on_off_T3_E_NW.onset(ind_special_T3(i));
%     
%     onsets_special_T3(i,7) = on_off_T3_NW_G.onset(ind_special_T3(i));
%     onsets_special_T3(i,8) = on_off_T3_NW_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,9) = on_off_T3_NW_NW.onset(ind_special_T3(i));
%     
%     % create array for offsets
%     offsets_special_T3(i,1) = on_off_T3_G_G.offset(ind_special_T3(i));
%     offsets_special_T3(i,2) = on_off_T3_G_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,3) = on_off_T3_G_NW.offset(ind_special_T3(i));
%     
%     offsets_special_T3(i,4) = on_off_T3_E_G.offset(ind_special_T3(i));
%     offsets_special_T3(i,5) = on_off_T3_E_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,6) = on_off_T3_E_NW.offset(ind_special_T3(i));
%     
%     offsets_special_T3(i,7) = on_off_T3_NW_G.offset(ind_special_T3(i));
%     offsets_special_T3(i,8) = on_off_T3_NW_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,9) = on_off_T3_NW_NW.offset(ind_special_T3(i));
%     
%     if ismember(ind_special_T3(i), ind_special_T3_G_G)
%         keep_onsets_special_T3(i,1) = 0;
%         keep_offsets_special_T3(i,1) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_G_E)
%         keep_onsets_special_T3(i,2) = 0;
%         keep_offsets_special_T3(i,2) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_G_NW)
%         keep_onsets_special_T3(i,3) = 0;
%         keep_offsets_special_T3(i,3) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_E_G)
%         keep_onsets_special_T3(i,4) = 0;
%         keep_offsets_special_T3(i,4) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_E_E)
%         keep_onsets_special_T3(i,5) = 0;
%         keep_offsets_special_T3(i,5) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_E_NW)
%         keep_onsets_special_T3(i,6) = 0;
%         keep_offsets_special_T3(i,6) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_G)
%         keep_onsets_special_T3(i,7) = 0;
%         keep_offsets_special_T3(i,7) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_E)
%         keep_onsets_special_T3(i,8) = 0;
%         keep_offsets_special_T3(i,8) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_NW)
%         keep_onsets_special_T3(i,9) = 0;
%         keep_offsets_special_T3(i,9) = 0;
%     end
% end
% 
% means_onset = floor(sum([onsets_special_T3.*keep_onsets_special_T3], 2)./sum(keep_onsets_special_T3,2));
% means_offset = floor(sum([offsets_special_T3.*keep_offsets_special_T3], 2)./sum(keep_offsets_special_T3,2));
% 
% % s = tps/500-0.102
% % tps = 500*s + 51
% means_onset_s = means_onset/500 - 0.102;
% means_offset_s = means_offset/500 - 0.102;
% 
% % replace onsets and offsets
% on_off_T3_G_G_new = on_off_T3_G_G;
% on_off_T3_G_E_new = on_off_T3_G_E;
% on_off_T3_G_NW_new = on_off_T3_G_NW;
% 
% on_off_T3_E_G_new = on_off_T3_E_G;
% on_off_T3_E_E_new = on_off_T3_E_E;
% on_off_T3_E_NW_new = on_off_T3_E_NW;
% 
% on_off_T3_NW_G_new = on_off_T3_NW_G;
% on_off_T3_NW_E_new = on_off_T3_NW_E;
% on_off_T3_NW_NW_new = on_off_T3_NW_NW;
% 
% for i = 1:length(ind_special_T3)    
%     if ismember(ind_special_T3(i), ind_special_T3_G_G)
%         %on_off_new with new onset for special cases
%         on_off_T3_G_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_G_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         % exchange onset/offset in latency struct for each special case
%         latency_T3_G_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_G_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_G_G, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_G_new.onset(ind_special_T3(i)):on_off_T3_G_G_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_G_G.areaLat(ind_special_T3(i)) = mean_latency_s;    
%     end
%     
%     if ismember(ind_special_T3(i), ind_special_T3_G_E)
%         on_off_T3_G_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_G_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_G_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_G_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_G_E, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_E_new.onset(ind_special_T3(i)):on_off_T3_G_E_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_G_E.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end 
%     
%     if ismember(ind_special_T3(i), ind_special_T3_G_NW)
%         on_off_T3_G_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_G_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_G_NW.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_G_NW.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_G_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_NW_new.onset(ind_special_T3(i)):on_off_T3_G_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_G_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end
%     
%     if ismember(ind_special_T3(i), ind_special_T3_E_G)
%         on_off_T3_E_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_E_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_E_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_E_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_E_G, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_G_new.onset(ind_special_T3(i)):on_off_T3_E_G_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_E_G.areaLat(ind_special_T3(i)) = mean_latency_s; 
%     end
%     
%     if ismember(ind_special_T3(i), ind_special_T3_E_E)
%         on_off_T3_E_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_E_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_E_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_E_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_E_E, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_E_new.onset(ind_special_T3(i)):on_off_T3_E_E_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_E_E.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end 
%     
%     if ismember(ind_special_T3(i), ind_special_T3_E_NW)
%         on_off_T3_E_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_E_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_E_NW.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_E_NW.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_E_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_NW_new.onset(ind_special_T3(i)):on_off_T3_E_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_E_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end 
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_G)
%         on_off_T3_NW_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_NW_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_NW_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_G, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_G_new.onset(ind_special_T3(i)):on_off_T3_NW_G_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_G.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end     
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_E)
%         on_off_T3_NW_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_NW_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_NW_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_E, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_E_new.onset(ind_special_T3(i)):on_off_T3_NW_E_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_E.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end    
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_NW)
%         on_off_T3_NW_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1); 
%         
%         latency_T3_NW_NW.onset(ind_special_T3(i)) =  means_onset(i,1);
%         latency_T3_NW_NW.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_NW_new.onset(ind_special_T3(i)):on_off_T3_NW_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end
% end


%% task 3 with 3 conditions
% column 1: G_G; column 2: G_E; column 3: E_G; column 4: E_E;

onsets_special_T3 = nan(length(ind_special_T3), 4);
offsets_special_T3 = nan(length(ind_special_T3), 4);

keep_onsets_special_T3 = ones(length(ind_special_T3), 4);
keep_offsets_special_T3 = ones(length(ind_special_T3), 4);

for i = 1:length(ind_special_T3)
    % create array for onsets
    onsets_special_T3(i,1) = on_off_T3_G_G.onset(ind_special_T3(i));
    onsets_special_T3(i,2) = on_off_T3_G_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,3) = on_off_T3_G_NW.onset(ind_special_T3(i));
    
    onsets_special_T3(i,3) = on_off_T3_E_G.onset(ind_special_T3(i));
    onsets_special_T3(i,4) = on_off_T3_E_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,6) = on_off_T3_E_NW.onset(ind_special_T3(i));
%     
%     onsets_special_T3(i,7) = on_off_T3_NW_G.onset(ind_special_T3(i));
%     onsets_special_T3(i,8) = on_off_T3_NW_E.onset(ind_special_T3(i));
%     onsets_special_T3(i,9) = on_off_T3_NW_NW.onset(ind_special_T3(i));
    
    % create array for offsets
    offsets_special_T3(i,1) = on_off_T3_G_G.offset(ind_special_T3(i));
    offsets_special_T3(i,2) = on_off_T3_G_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,3) = on_off_T3_G_NW.offset(ind_special_T3(i));
    
    offsets_special_T3(i,3) = on_off_T3_E_G.offset(ind_special_T3(i));
    offsets_special_T3(i,4) = on_off_T3_E_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,6) = on_off_T3_E_NW.offset(ind_special_T3(i));
%     
%     offsets_special_T3(i,7) = on_off_T3_NW_G.offset(ind_special_T3(i));
%     offsets_special_T3(i,8) = on_off_T3_NW_E.offset(ind_special_T3(i));
%     offsets_special_T3(i,9) = on_off_T3_NW_NW.offset(ind_special_T3(i));
    
    if ismember(ind_special_T3(i), ind_special_T3_G_G)
        keep_onsets_special_T3(i,1) = 0;
        keep_offsets_special_T3(i,1) = 0;
    end
    if ismember(ind_special_T3(i), ind_special_T3_G_E)
        keep_onsets_special_T3(i,2) = 0;
        keep_offsets_special_T3(i,2) = 0;
    end
%     if ismember(ind_special_T3(i), ind_special_T3_G_NW)
%         keep_onsets_special_T3(i,3) = 0;
%         keep_offsets_special_T3(i,3) = 0;
%     end
    if ismember(ind_special_T3(i), ind_special_T3_E_G)
        keep_onsets_special_T3(i,3) = 0;
        keep_offsets_special_T3(i,3) = 0;
    end
    if ismember(ind_special_T3(i), ind_special_T3_E_E)
        keep_onsets_special_T3(i,4) = 0;
        keep_offsets_special_T3(i,4) = 0;
    end
%     if ismember(ind_special_T3(i), ind_special_T3_E_NW)
%         keep_onsets_special_T3(i,6) = 0;
%         keep_offsets_special_T3(i,6) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_G)
%         keep_onsets_special_T3(i,7) = 0;
%         keep_offsets_special_T3(i,7) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_E)
%         keep_onsets_special_T3(i,8) = 0;
%         keep_offsets_special_T3(i,8) = 0;
%     end
%     if ismember(ind_special_T3(i), ind_special_T3_NW_NW)
%         keep_onsets_special_T3(i,9) = 0;
%         keep_offsets_special_T3(i,9) = 0;
%     end
end

means_onset = floor(sum([onsets_special_T3.*keep_onsets_special_T3], 2)./sum(keep_onsets_special_T3,2));
means_offset = floor(sum([offsets_special_T3.*keep_offsets_special_T3], 2)./sum(keep_offsets_special_T3,2));

% s = tps/500-0.102
% tps = 500*s + 51
means_onset_s = means_onset/500 - 0.102;
means_offset_s = means_offset/500 - 0.102;

% replace onsets and offsets
on_off_T3_G_G_new = on_off_T3_G_G;
on_off_T3_G_E_new = on_off_T3_G_E;
% on_off_T3_G_NW_new = on_off_T3_G_NW;

on_off_T3_E_G_new = on_off_T3_E_G;
on_off_T3_E_E_new = on_off_T3_E_E;
% on_off_T3_E_NW_new = on_off_T3_E_NW;
% 
% on_off_T3_NW_G_new = on_off_T3_NW_G;
% on_off_T3_NW_E_new = on_off_T3_NW_E;
% on_off_T3_NW_NW_new = on_off_T3_NW_NW;

for i = 1:length(ind_special_T3) 
        if means_onset(i) > 0;
    if ismember(ind_special_T3(i), ind_special_T3_G_G)
        %on_off_new with new onset for special cases
        on_off_T3_G_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
        on_off_T3_G_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
        
        % exchange onset/offset in latency struct for each special case
        latency_T3_G_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
        latency_T3_G_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T3_G_G, 2));
        erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_G_new.onset(ind_special_T3(i)):on_off_T3_G_G_new.offset(ind_special_T3(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T3_G_G.areaLat(ind_special_T3(i)) = mean_latency_s;    
    end
    
    if ismember(ind_special_T3(i), ind_special_T3_G_E)
        on_off_T3_G_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
        on_off_T3_G_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
        
        latency_T3_G_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
        latency_T3_G_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T3_G_E, 2));
        erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_E_new.onset(ind_special_T3(i)):on_off_T3_G_E_new.offset(ind_special_T3(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T3_G_E.areaLat(ind_special_T3(i)) = mean_latency_s;
    end 
    
%     if ismember(ind_special_T3(i), ind_special_T3_G_NW)
%         on_off_T3_G_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_G_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_G_NW.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_G_NW.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_G_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_G_NW_new.onset(ind_special_T3(i)):on_off_T3_G_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_G_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end
    
    if ismember(ind_special_T3(i), ind_special_T3_E_G)
        on_off_T3_E_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
        on_off_T3_E_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
        
        latency_T3_E_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
        latency_T3_E_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T3_E_G, 2));
        erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_G_new.onset(ind_special_T3(i)):on_off_T3_E_G_new.offset(ind_special_T3(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T3_E_G.areaLat(ind_special_T3(i)) = mean_latency_s; 
    end
    
    if ismember(ind_special_T3(i), ind_special_T3_E_E)
        on_off_T3_E_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
        on_off_T3_E_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
        
        latency_T3_E_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
        latency_T3_E_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
        
        % recalculate mean latency for special cases
        erp_mean = squeeze(mean(T3_E_E, 2));
        erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_E_new.onset(ind_special_T3(i)):on_off_T3_E_E_new.offset(ind_special_T3(i)));
        [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
        mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
        latency_T3_E_E.areaLat(ind_special_T3(i)) = mean_latency_s;
    end 
    
%     if ismember(ind_special_T3(i), ind_special_T3_E_NW)
%         on_off_T3_E_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_E_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_E_NW.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_E_NW.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_E_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_E_NW_new.onset(ind_special_T3(i)):on_off_T3_E_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_E_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end 
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_G)
%         on_off_T3_NW_G_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_G_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_NW_G.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_NW_G.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_G, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_G_new.onset(ind_special_T3(i)):on_off_T3_NW_G_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_G.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end     
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_E)
%         on_off_T3_NW_E_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_E_new.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         latency_T3_NW_E.onset(ind_special_T3(i)) =  means_onset_s(i,1);
%         latency_T3_NW_E.offset(ind_special_T3(i)) =  means_offset_s(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_E, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_E_new.onset(ind_special_T3(i)):on_off_T3_NW_E_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_E.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end    
%     
%     if ismember(ind_special_T3(i), ind_special_T3_NW_NW)
%         on_off_T3_NW_NW_new.onset(ind_special_T3(i)) =  means_onset(i,1);
%         on_off_T3_NW_NW_new.offset(ind_special_T3(i)) =  means_offset(i,1); 
%         
%         latency_T3_NW_NW.onset(ind_special_T3(i)) =  means_onset(i,1);
%         latency_T3_NW_NW.offset(ind_special_T3(i)) =  means_offset(i,1);
%         
%         % recalculate mean latency for special cases
%         erp_mean = squeeze(mean(T3_NW_NW, 2));
%         erp_mean = erp_mean(ind_special_T3(i),on_off_T3_NW_NW_new.onset(ind_special_T3(i)):on_off_T3_NW_NW_new.offset(ind_special_T3(i)));
%         [~,latency_index] = min(abs(cumsum(erp_mean)-sum(erp_mean)/2));
%         mean_latency_s = (means_onset_s(i,1) + (latency_index/500));
%         latency_T3_NW_NW.areaLat(ind_special_T3(i)) = mean_latency_s;
%     end
        end
end

%% mean amplitude von onset zu offset von Hand berechnen
for subi= 1:nSubs
    latency_T1_G_HF.meanAmp(subi,1) = mean(T1_G_HF(subi,ind_chans,on_off_T1_G_HF_new.onset(subi):on_off_T1_G_HF_new.offset(subi)),[3 2]);
    latency_T1_G_LF.meanAmp(subi,1) = mean(T1_G_LF(subi,ind_chans,on_off_T1_G_LF_new.onset(subi):on_off_T1_G_LF_new.offset(subi)),[3 2]);
    latency_T1_G_NW.meanAmp(subi,1) = mean(T1_G_NW(subi,ind_chans,on_off_T1_G_NW_new.onset(subi):on_off_T1_G_NW_new.offset(subi)),[3 2]);

    latency_T2_E_HF.meanAmp(subi,1) = mean(T2_E_HF(subi,ind_chans,on_off_T2_E_HF_new.onset(subi):on_off_T2_E_HF_new.offset(subi)),[3 2]);
    latency_T2_E_LF.meanAmp(subi,1) = mean(T2_E_LF(subi,ind_chans,on_off_T2_E_LF_new.onset(subi):on_off_T2_E_LF_new.offset(subi)),[3 2]);
    latency_T2_E_NW.meanAmp(subi,1) = mean(T2_E_NW(subi,ind_chans,on_off_T2_E_NW_new.onset(subi):on_off_T2_E_NW_new.offset(subi)),[3 2]);

    latency_T3_G_G.meanAmp(subi,1) = mean(T3_G_G(subi,ind_chans,on_off_T3_G_G_new.onset(subi):on_off_T3_G_G_new.offset(subi)),[3 2]);
    latency_T3_G_E.meanAmp(subi,1) = mean(T3_G_E(subi,ind_chans,on_off_T3_G_E_new.onset(subi):on_off_T3_G_E_new.offset(subi)),[3 2]);
%     latency_T3_G_NW.meanAmp(subi,1) = mean(T3_G_NW(subi,ind_chans,on_off_T3_G_NW_new.onset(subi):on_off_T3_G_NW_new.offset(subi)),[3 2]);
    
    latency_T3_E_G.meanAmp(subi,1) = mean(T3_E_G(subi,ind_chans,on_off_T3_E_G_new.onset(subi):on_off_T3_E_G_new.offset(subi)),[3 2]);
    latency_T3_E_E.meanAmp(subi,1) = mean(T3_E_E(subi,ind_chans,on_off_T3_E_E_new.onset(subi):on_off_T3_E_E_new.offset(subi)),[3 2]);
%     latency_T3_E_NW.meanAmp(subi,1) = mean(T3_E_NW(subi,ind_chans,on_off_T3_E_NW_new.onset(subi):on_off_T3_E_NW_new.offset(subi)),[3 2]);
% 
%     latency_T3_NW_G.meanAmp(subi,1) = mean(T3_NW_G(subi,ind_chans,on_off_T3_NW_G_new.onset(subi):on_off_T3_NW_G_new.offset(subi)),[3 2]);
%     latency_T3_NW_E.meanAmp(subi,1) = mean(T3_NW_E(subi,ind_chans,on_off_T3_NW_E_new.onset(subi):on_off_T3_NW_E_new.offset(subi)),[3 2]);
%     latency_T3_NW_NW.meanAmp(subi,1) = mean(T3_NW_NW(subi,ind_chans,on_off_T3_NW_NW_new.onset(subi):on_off_T3_NW_NW_new.offset(subi)),[3 2]);
end


%% save files
cd(dataPath);
filename = 'output_liesefeld.mat';
save(filename, 'latency_T1_G_HF', 'latency_T1_G_LF', 'latency_T1_G_NW', ...
    'latency_T2_E_HF', 'latency_T2_E_LF', 'latency_T2_E_NW', ...
    'latency_T3_G_G', 'latency_T3_G_E', 'latency_T3_G_NW', ...
    'latency_T3_E_G', 'latency_T3_E_E', 'latency_T3_E_NW', ...
    'latency_T3_NW_G', 'latency_T3_NW_E', 'latency_T3_NW_NW', ...
    '-v7.3')


%% conver to table
%load('output_liesefeld.mat');

% get vpNames for each task
name = 'vpNames';
vpNames = (all_erp_T1_HF(1,:))';

T_latency_T1_HF = struct2table(latency_T1_G_HF);T_latency_T1_HF.(name) = vpNames;T_latency_T1_HF = [T_latency_T1_HF(:,11), T_latency_T1_HF(:,1:10)];
T_latency_T1_LF = struct2table(latency_T1_G_LF);T_latency_T1_LF.(name) = vpNames;T_latency_T1_LF = [T_latency_T1_LF(:,11), T_latency_T1_LF(:,1:10)];
T_latency_T1_NW = struct2table(latency_T1_G_NW);T_latency_T1_NW.(name) = vpNames;T_latency_T1_NW = [T_latency_T1_NW(:,11), T_latency_T1_NW(:,1:10)];

T_latency_T2_HF = struct2table(latency_T2_E_HF);T_latency_T2_HF.(name) = vpNames;T_latency_T2_HF = [T_latency_T2_HF(:,11), T_latency_T2_HF(:,1:10)];
T_latency_T2_LF = struct2table(latency_T2_E_LF);T_latency_T2_LF.(name) = vpNames;T_latency_T2_LF = [T_latency_T2_LF(:,11), T_latency_T2_LF(:,1:10)];
T_latency_T2_NW = struct2table(latency_T2_E_NW);T_latency_T2_NW.(name) = vpNames;T_latency_T2_NW = [T_latency_T2_NW(:,11), T_latency_T2_NW(:,1:10)];

T_latency_T3_G_G  = struct2table(latency_T3_G_G);T_latency_T3_G_G.(name) = vpNames;T_latency_T3_G_G = [T_latency_T3_G_G(:,11), T_latency_T3_G_G(:,1:10)];
T_latency_T3_G_E  = struct2table(latency_T3_G_E);T_latency_T3_G_E.(name) = vpNames;T_latency_T3_G_E = [T_latency_T3_G_E(:,11), T_latency_T3_G_E(:,1:10)];
% T_latency_T3_G_NW = struct2table(latency_T3_G_NW);T_latency_T3_G_NW.(name) = vpNames;T_latency_T3_G_NW = [T_latency_T3_G_NW(:,11), T_latency_T3_G_NW(:,1:10)];

T_latency_T3_E_G  = struct2table(latency_T3_E_G);T_latency_T3_E_G.(name) = vpNames;T_latency_T3_E_G = [T_latency_T3_E_G(:,11), T_latency_T3_E_G(:,1:10)];
T_latency_T3_E_E  = struct2table(latency_T3_E_E);T_latency_T3_E_E.(name) = vpNames;T_latency_T3_E_E = [T_latency_T3_E_E(:,11), T_latency_T3_E_E(:,1:10)];
% T_latency_T3_E_NW = struct2table(latency_T3_E_NW);T_latency_T3_E_NW.(name) = vpNames;T_latency_T3_E_NW = [T_latency_T3_E_NW(:,11), T_latency_T3_E_NW(:,1:10)];
% 
% T_latency_T3_NW_G  = struct2table(latency_T3_NW_G);T_latency_T3_NW_G.(name) = vpNames;T_latency_T3_NW_G = [T_latency_T3_NW_G(:,11), T_latency_T3_NW_G(:,1:10)];
% T_latency_T3_NW_E  = struct2table(latency_T3_NW_E);T_latency_T3_NW_E.(name) = vpNames;T_latency_T3_NW_E = [T_latency_T3_NW_E(:,11), T_latency_T3_NW_E(:,1:10)];
% T_latency_T3_NW_NW = struct2table(latency_T3_NW_NW);T_latency_T3_NW_NW.(name) = vpNames;T_latency_T3_NW_NW = [T_latency_T3_NW_NW(:,11), T_latency_T3_NW_NW(:,1:10)];

% save individual tables as .txt
writetable(T_latency_T1_HF);
writetable(T_latency_T1_LF);
writetable(T_latency_T1_NW);

writetable(T_latency_T2_HF);
writetable(T_latency_T2_LF);
writetable(T_latency_T2_NW);

writetable(T_latency_T3_G_G);
writetable(T_latency_T3_G_E);
% writetable(T_latency_T3_G_NW);
writetable(T_latency_T3_E_G);
writetable(T_latency_T3_E_E);
% writetable(T_latency_T3_E_NW);
% writetable(T_latency_T3_NW_G);
% writetable(T_latency_T3_NW_E);
% writetable(T_latency_T3_NW_NW);

cd(fileparts(tmp.Filename));