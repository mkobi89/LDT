##########################################################
##             DDM DATASETS WITHOUT OUTLIERS            ##
##########################################################
## Description :: removes outlier trials for each subject
##                and saves data as .dat files for each 
##                participant and task for DDM analysis
## Input :::::::: LDT dataframe (01_answer_files_to_LDT.R)
##                RT_exclude_trials (03_RT_outliers.R)
## Libraries :::: dplyr, tidyverse
## Output ::::::: .dat file per task per subject
##########################################################

## libraries
library(tidyverse)

## data path
saveFolder_ddm     <- file.path("data/dm_datasets")

## get subjects and task names
Subjects <- unique(LDT$subjectID)
Tasks <- unique(LDT$which_task)

## load tabel with outliers to exclude
all_trials_exclude <- read.csv("data/outliers/RT_exclude_trials.csv")

## create Dataset for each subject for german and english task
vp_data <- LDT %>% 
  group_by(subjectID) %>%
  select(subjectID, which_task, word_nonword, word_pseudoword_numeric, frequency, correct, RT) #%>% 
#write.table()

## process subjects task 1
for(i in Subjects){
  for(t in Tasks[1]){
    # Get data from a single subject with variables of interest
    dm_data <- vp_data %>%
      filter(subjectID == i, which_task == t)
    
    # create table with outliers
    outliers <- all_trials_exclude %>%
      filter(subjectID == i, which_task == t) %>%
      select(subjectID, which_task, TrialNumber)
    
    # remove trials from dm_data set
    if (nrow(outliers) >= 1) {
      dm_data <- dm_data[-c(outliers$TrialNumber),]
    }
    
    # Create Filename for each dataset
    saveas <- paste(i,t,"_final.dat", sep="")
    
    # Write File for each Subject's data
    savePath = paste(saveFolder_ddm, saveas, sep = "/")
    write.table(dm_data, savePath, row.names = F, col.names = F, sep = " ")
  }
}

## task 2
for(i in Subjects){
  for(t in Tasks[2]){
      # Get data from a single subject with variables of interest
      dm_data <- vp_data %>%
        filter(subjectID == i, which_task == t)
      
      # create table with outliers
      outliers <- all_trials_exclude %>%
        filter(subjectID == i, which_task == t) %>%
        select(subjectID, which_task, TrialNumber)
      
      # remove trials from dm_data set
      if (nrow(outliers) >= 1) {
        dm_data <- dm_data[-c(outliers$TrialNumber),]
      }
      
      # Create Filename for each dataset
      saveas <- paste(i,t,"_final.dat", sep="")
      
      # Write File for each Subject's data
      savePath = paste(saveFolder_ddm, saveas, sep = "/")
      write.table(dm_data, savePath, row.names = F, col.names = F, sep = " ")
  }
}


## task 3
vp_data <- LDT %>% 
  group_by(subjectID) %>%
  select(subjectID, which_task, word_nonword, word_pseudoword_numeric, condition, correct, RT)


for(i in Subjects){
  for(t in Tasks[3]){
    if(i != "C50"){
      # Get data from a single subject with variables of interest
      dm_data <- vp_data %>%
        filter(subjectID == i, which_task == t)
      
      # create table with outliers
      outliers <- all_trials_exclude %>%
        filter(subjectID == i, which_task == t) %>%
        select(subjectID, which_task, TrialNumber)
      
      # remove trials from dm_data set
      if (nrow(outliers) >= 1) {
        dm_data <- dm_data[-c(outliers$TrialNumber),]
      }
      
      # Create Filename for each dataset
      saveas <- paste(i,t,"_final.dat", sep="")
      
      # Write File for each Subject's data
      savePath = paste(saveFolder_ddm, saveas, sep = "/")
      write.table(dm_data, savePath, row.names = F, col.names = F, sep = " ")
    }
  }
}

## clean workspace
rm(all_trials_exclude, dm_data, vp_data, outliers,
   i, saveas, saveFolder_ddm, savePath, t, Tasks, Subjects)
