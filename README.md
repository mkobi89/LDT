# LDT word frequancy and language switching with highly proficient bilinguals

This repository is used for the analysis of the LDT data.
For additional information, write to: matthias.kobi@gmail.com / matthias.kobi@uzh.ch

Content description:

data/rawdata:
- answers: contains all LDT behavioral data from participants in .csv (extracted from "Fullresults*.mat" in "matlab_rawdata" by "Matlab/LDT_BHD_preparing.m")
- matlab_rawdata: contains all MATLAB answer files from CLINT task
- psychometrics.csv: psychometrics rawdata for participants in LDT
- stimuli.csv: word frequency ratings for stimuli used in LDT

data/dm_datasets:
- .dat files: behavioral data for all participants after trimming (RT < 0.25 s and RT > 2 s excluded)
- fast-dm.exe: fast-dm modeling used in "R/05_DDM_parameter_estimation.R"
- .ctl files: parameter settings for fast-dm for all 3 tasks
- .txt files: results from fast-dm

data/erp:
- .txt files: results from ERP processing using Liesefeld (2018)

data/outliers:
- RT_exclude_trials.csv: all trials to be excluded as used in "R/04_DDM_generate_datasets_without_outliers.R"

data:
- all_tasks.RData: all data combined from psychometrics, DDM and ERP
- LDT.Rdata: all answers from all trials in all tasks and all participants combined
- LDT_clean.Rdata: same as LDT.Rdata but without outliers (RT < 0.25 s and RT > 2 s excluded)
- psychometrics.Rdata: psychometric data for participants

figures: contains all figures used for publication

Matlab:
- headmodeling: channel location files used for plots
- latency: matlab function (latency.m) to extract ERP features (Liesefeld, 2018)
- EEG_01_taskwise_segmentation_for_automagic.m: prepare raw EEG data of LDT ("All_data/C*/C*_LDT_EEG.mat") and split them into separate files per task for automagic preprocessing
- EEG_02_merge_data_to_fieldtrip_quality_assessment.m: merge preprocessed EEG data and convert them to fieldtrip
- EEG_03_trials_vis_rejc_art.m: use "ft_trialfun_LDT_automagic.m" for trialing ERP data and apply fieldtrip visual artefact rejection
- EEG_04_erp.m: calculate ERPs for all participants and all conditions
- EEG_05_liesefeld.m: apply Liesefeld approach to extract ERP features
- EEG_06_plotting.m: creates topoplot and plots of grand averages
- LDT_BHD_preparing.m: copy all "Fullanswers*.mat" in "All_Data" to "LDT/data/rawdata/matlab_rawdata" and write .csv answer files to "LDT/data/rawdata"
- ft_trialfun_LDT_automagic.m: Fieldtrip trialfunction used for LDT

R: contains all RSkripts
- master.R: controls all R files, gets and organises data and excludes unwanted participants
- 01_answer_files_to_LDT.R: processes .csv files in "LDT/data/rawdata" to "LDT.RData" in "LDT/data"
- 02_psychometrics.R: processes psychometrics.csv in "LDT/data/rawdata" to "psychometrics.RData" in "LDT/data"
- 03_RT_outliers.R: trims RT data from LDT.RData and creates "LDT_clean.RData" in "LDT/data"
- 04_DDM_generate_datasets_without_outliers.R: creates .dat files in "data/dm_datasets" for every participant after trimming (RT < 0.25 s and RT > 2 s excluded)
- 05_DDM_parameter_estimation.R: uses .dat files in "data/dm_datasets" to estimate DDM with "fast-dm.exe" using settings in .ctr files
- 06_final_table.R: collects all relevant data and creates "all_Tasks.RData"
- 07_summarise_tables.R: creates psychometrics table, stimuli table, and behavioral data tables for manuscript in figure folder
- 08_plots.R: creates plots for manuscript in figure folder
- 09_statistics.R: calculates all statistis, creates LMM tables and all combined plot in figure folter

stimulus material: contains stimuli used in LDT, "stimuli_word_frequency_task.csv" for both word frequency tasks, "stimuli_switching_task.csv" for switch task 

participant_selection.xlsx: contains information to exclude participants based on psychometrics and difficulties with erp extraction 

EEG: EGI geodesics
- sampling rate: 500 Hz
- electrodes: 128
- reference electrode: Cz
