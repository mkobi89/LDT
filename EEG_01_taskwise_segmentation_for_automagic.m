%% Segment LDT EEG data into each task

%% Preparation:
clear
%% Change directory to the current m. file
% by hand or by code:
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Path to server depending on OS

% if IsOSX == true
%     dataPath = '/Volumes/CLINT/All_Data/'; % All Data Path
%     savePath = '/Volumes/CLINT/LDT_segmented/'; % LDT results Path
%     addpath('/Volumes/CLINT/eeglab2020_0/'); % add eeglab path
% else
%     dataPath = '//130.60.235.123/users/neuro/Desktop/CLINT/All_Data/'; % All Data Path
%     savePath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT_segmented/'; % LDT Path
%     addpath('//130.60.235.123/users/neuro/Desktop/CLINT/eeglab2020_0/'); % add eeglab path
% end

%% Corona Homeoffice 2.0 path definition
dataPath = '//130.60.235.123/users/neuro/Desktop/CLINT/All_Data/'; % All Data Path
savePath = '//130.60.235.123/users/neuro/Desktop/CLINT/LDT/';
%% Add EEGlab
addpath('//130.60.235.123/users/neuro/Desktop/CLINT/eeglab2021_0/'); % add eeglab path

%% Starting EEGLAB
eeglab
close()

%% Select desired subjects

cd(dataPath);

vplist = dir('C*');
vpNames = {vplist.name};

%% Check for target directory and create folders

for zz = 1:length(vpNames)
    if ~(exist([savePath vpNames{zz}])==7)
        mkdir([savePath vpNames{zz}]);
    else
        continue
    end
end

%% Segment data
%
for zz = 1:length(vpNames)

    if exist([savePath vpNames{zz} '/' vpNames{zz} '_LDT_EEG_t1.mat'])==2
        continue
    else
        cd([dataPath vpNames{zz}]);
        
        file = dir('*LDT_EEG.mat');
        file = file.name;
        
        load(file)
        
    %% adjust to trigger delay
    % trippers are 18ms late with N300 Amp at 500Hz sampling rate
    % see: https://sccn.ucsd.edu/pipermail/eeglablist/2015/010429.html for
    % more information 
    
    srate = EEG.srate;
    for index = 1:length(EEG.event)
        EEG.event(index).latency = EEG.event(index).latency-0.018*srate;
    end


        
        %% find latencies of start and end point of the 3 tasks
        for i = 1:length(EEG.event)
            if isequal(EEG.event(i).type, '190 ')
                t1_start = i; %  position of the event
                t1_start_latency = EEG.event(i).latency; % latenncy of the event
                
            elseif isequal(EEG.event(i).type, '290 ')
                t2_start = i;
                t2_start_latency = EEG.event(i).latency;
                
            elseif isequal(EEG.event(i).type, '390 ')
                t3_start = i;
                t3_start_latency = EEG.event(i).latency;
                
            elseif isequal(EEG.event(i).type, '195 ')
                t1_end = i;
                t1_end_latency = EEG.event(i).latency;
                
            elseif isequal(EEG.event(i).type, '295 ')
                t2_end = i;
                t2_end_latency = EEG.event(i).latency;
                
            elseif isequal(EEG.event(i).type, '395 ')
                t3_end = i;
                t3_end_latency = EEG.event(i).latency;
            end
        end
        
        %% Cut EEG file in 3 separate segments
        % everything from 0 until start latency - 5 is removed and
        % everything from end latency + 5 until end of eeg data is removed
        EEG_t1 = eeg_eegrej(EEG, [0 t1_start_latency-5; t1_end_latency+1005 EEG.pnts+1]);
        EEG_t2 = eeg_eegrej(EEG, [0 t2_start_latency-5; t2_end_latency+1005 EEG.pnts+1]);
        EEG_t3 = eeg_eegrej(EEG, [0 t3_start_latency-5; t3_end_latency+1005 EEG.pnts+1]);
        
        %% remove boundary event
        del = 0;
        for i = 1:length(EEG_t1.event)
            if isequal(EEG_t1.event(i-del).type, 'boundary')
                EEG_t1.event(i-del) = [];
                del = del+1;
            end
        end
        
        del = 0;
        for i = 1:length(EEG_t2.event)
            if isequal(EEG_t2.event(i-del).type, 'boundary')
                EEG_t2.event(i-del) = [];
                del = del+1;
            end
        end
        
        del = 0;
        for i = 1:length(EEG_t3.event)
            if isequal(EEG_t3.event(i-del).type, 'boundary')
                EEG_t3.event(i-del) = [];
                del = del+1;
            end
        end
        
        %% save each segment as a new file
        cd([savePath vpNames{zz}]);
                
        EEG = EEG_t1;
        filename = strjoin([vpNames(zz) '_LDT_EEG_t1.mat'], '');
        save(filename, 'EEG');
        
        EEG = EEG_t2;
        filename = strjoin([vpNames(zz) '_LDT_EEG_t2.mat'], '');
        save(filename, 'EEG');
        
        EEG = EEG_t3;
        filename = strjoin([vpNames(zz) '_LDT_EEG_t3.mat'], '');
        save(filename, 'EEG');
        
        disp(['****** Segmentation of ' vpNames{zz} ' is done. ******' ])
        
    end
end

cd(fileparts(tmp.Filename));

