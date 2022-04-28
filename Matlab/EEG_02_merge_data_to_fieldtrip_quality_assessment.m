%% Merge preprocessed LDT EEG data from each task into one file in fieldtrip structure
% get Automagic quality assessment
% save files to Server

%% Preparation:
clear

%% Change directory to the current m. file
% by hand or by code:
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Path to server depending on OS

% if IsOSX == true
%     dataPath = '/Volumes/CLINT/All_Data/'; % All Data Path
%     savePath = '/Volumes/CLINT/LDT_preprocessed/'; % LDT Path
%     addpath('/Volumes/CLINT/eeglab2020_0/'); % add eeglab path
%     addpath('/Volumes/CLINT/fieldtrip-20200607/'); % add fieldtrip path
% else
%     dataPath = '//130.60.235.123/users/neuro/Desktop/CLINT/All_Data/'; % All Data Path
%     savePath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT_preprocessed/'; % LDT Path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/eeglab2020_0/'); % add eeglab path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/fieldtrip-20200607/'); % add fieldtrip path
% end

% running on server
dataPath = 'C:\Users\neuro\Desktop\CLINT\LDT_results\';
savePath = 'C:\Users\neuro\Desktop\CLINT\LDT_preprocessed\';

% Add EEGLAB and Fieldtrip
addpath('C:\Users/neuro/Desktop/CLINT/fieldtrip-20200607/'); % add fieldtrip path
addpath('C:\Users\neuro/Desktop/CLINT/eeglab2021_0/'); % add eeglab path

%% Starting EEGLAB and Fieldtrip
eeglab
close()

ft_defaults;

%% Select desired subjects
cd(dataPath);

vplist = dir('C*');
vpNames = {vplist.name};

%% Check for target directory and create folders to save files

for zz = 1:length(vpNames)
    if ~(exist([savePath vpNames{zz}])==7)
        mkdir([savePath vpNames{zz}]);
    else
        continue
    end
end


%% Merge data

 for zz = 1:length(vpNames)
    if exist([savePath vpNames{zz} '/' vpNames{zz} '_LDT_preprocessed_18.mat'])==2
        continue
    else
     
    cd([dataPath vpNames{zz}]);

    
    file = dir('*p*_EEG_t1.mat');
    EEG_t1 = load(file.name);
    
    file = dir('*p*_EEG_t2.mat');
    EEG_t2 = load(file.name);
    
    file = dir('*p*_EEG_t3.mat');
    EEG_t3 = load(file.name);
    
    EEG_merged_1 = pop_mergeset(EEG_t1.EEG, EEG_t2.EEG);
    EEG_final = pop_mergeset(EEG_merged_1, EEG_t3.EEG);
    
    del = 0;
    for i = 1:length(EEG_final.event)
        if isequal(EEG_final.event(i-del).type, 'boundary')
            EEG_final.event(i-del) = [];
            del = del+1;
        end
    end
    
    EEG_final.setname = 'EGI file';
    EEG = EEG_final;
    
    
    
    
    %% convert to fieldtrip
    
    data = eeglab2fieldtrip(EEG, 'raw','none');
    data.hdr = ft_fetch_header(data);
    
    
    %% take elec --> unused
    %         ind=ismember(elec_aligned.label,data.label);
    %         elec.chanpos= elec_aligned.chanpos(find(ind==1),:);
    %         elec.chantype= elec_aligned.chantype(find(ind==1));
    %         elec.chanunit= elec_aligned.chanunit(find(ind==1));
    %         elec.elecpos= elec_aligned.elecpos(find(ind==1),:);
    %         elec.homogeneous= elec_aligned.homogeneous;
    %         elec.label= elec_aligned.label(find(ind==1));
    %         elec.type= elec_aligned.type;
    %         elec.unit= elec_aligned.unit;
    %
    
    % create hdr for ft data; dimord for EEGLAB should be chan x time x trial double
    %         hdr.Fs = data.fsample;
    %         hdr.chantype = elec.chantype;
    %         hdr.chanunit = elec.chanunit;
    %         hdr.label    = elec.label;
    %         hdr.nChans   = numel(elec.label);
    %         hdr.elec = elec;
    %         hdr.nTrials = numel(data.trial);
    %         hdr.nSamples = length(data.time{1});
    %         hdr.nSamplesPre = 2*data.fsample;
    
    %% create fieldtrip event structure
    for fn = 1:numel(EEG.event)
        data.event(fn).sample = EEG.event(fn).latency;
        data.event(fn).offset = [];
        data.event(fn).duration = 0;
        data.event(fn).type = 'trigger';
        data.event(fn).value = str2double(EEG.event(fn).type(1:3));
    end
    
    %% get quality rates
    
    quality_scores{zz,1} = vpNames{zz};
    quality_scores{zz,2} = EEG_t1.automagic.rate;
    quality_scores{zz,3} = EEG_t2.automagic.rate;
    quality_scores{zz,4} = EEG_t3.automagic.rate;
    
    %% add quality scores to fieltrip data structure
    
    data.automagic_rate.t1 = EEG_t1.automagic.rate;
    data.automagic_rate.t2 = EEG_t2.automagic.rate;
    data.automagic_rate.t3 = EEG_t3.automagic.rate;
    
    % save fieldtrip data

   save([savePath vpNames{zz} '/'  vpNames{zz} '_LDT_preprocessed_18.mat'], 'data', '-v7.3');
    
    
    clearvars EEG EEG_final EEG_t1 EEG_t2 EEG_t3 data
    end
 end

% arrange quality score as table and save it
quality_scores = cell2table(quality_scores, 'VariableNames',{'SubjectID' 'quality_t1' 'quality_t2' 'quality_t3'});
writetable(quality_scores,[savePath 'quality_scores.csv'],'FileType','text')

cd(fileparts(tmp.Filename));
