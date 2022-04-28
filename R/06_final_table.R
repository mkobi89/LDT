##########################################################
##            COLLECT PREPROCESSED DATA                 ##
##########################################################
## Description :: gathers all preprocessed data (EEG, DDM
##                and psychometrics) in long format
##                for further analysis
## Input :::::::: psychometrics (02_psychometrics.R)
##                .txt files for DDM parameters (05_DDM_parameter_estimation.R)
##                .txt files for Liesefeld parameters in data/erp
## Libraries :::: tidyverse
## Output ::::::: all_tasks.Rdata
##########################################################

## libraries
library(tidyverse)

## data path
DDMFolder    <- file.path("data/dm_datasets")
liesefeldFolder    <- file.path("data/erp")

## read DDM parameters
DDM_parms_T1_G_hf_lf <- read.table(paste(DDMFolder, "parms_G_hf_lf_final.txt", sep = "/"), header = TRUE)
DDM_parms_T2_E_hf_lf <- read.table(paste(DDMFolder, "parms_E_hf_lf_final.txt", sep = "/"), header = TRUE)
DDM_parms_T3_SW_Condition <- read.table(paste(DDMFolder, "parms_SW_condition_final.txt", sep = "/"), header = TRUE)

## read ERP parameters
latencies_T1_G_HF <- read.table(paste(liesefeldFolder, "T_latency_T1_HF.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T1_G_LF <- read.table(paste(liesefeldFolder, "T_latency_T1_LF.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T1_G_NW <- read.table(paste(liesefeldFolder, "T_latency_T1_NW.txt", sep = "/"), sep = ",", header = TRUE)

latencies_T2_E_HF <- read.table(paste(liesefeldFolder, "T_latency_T2_HF.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T2_E_LF <- read.table(paste(liesefeldFolder, "T_latency_T2_LF.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T2_E_NW <- read.table(paste(liesefeldFolder, "T_latency_T2_NW.txt", sep = "/"), sep = ",", header = TRUE)

latencies_T3_SW_G_G <- read.table(paste(liesefeldFolder, "T_latency_T3_G_G.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T3_SW_G_E <- read.table(paste(liesefeldFolder, "T_latency_T3_G_E.txt", sep = "/"), sep = ",", header = TRUE)
#latencies_T3_SW_G_NW <- read.table(paste(liesefeldFolder, "T_latency_T3_G_NW.txt", sep = "/"), sep = ",", header = TRUE)

latencies_T3_SW_E_G <- read.table(paste(liesefeldFolder, "T_latency_T3_E_G.txt", sep = "/"), sep = ",", header = TRUE)
latencies_T3_SW_E_E <- read.table(paste(liesefeldFolder, "T_latency_T3_E_E.txt", sep = "/"), sep = ",", header = TRUE)
#latencies_T3_SW_E_NW <- read.table(paste(liesefeldFolder, "T_latency_T3_E_NW.txt", sep = "/"), sep = ",", header = TRUE)

#latencies_T3_SW_NW_G <- read.table(paste(liesefeldFolder, "T_latency_T3_NW_G.txt", sep = "/"), sep = ",", header = TRUE)
#latencies_T3_SW_NW_E <- read.table(paste(liesefeldFolder, "T_latency_T3_NW_E.txt", sep = "/"), sep = ",", header = TRUE)
#latencies_T3_SW_NW_NW <- read.table(paste(liesefeldFolder, "T_latency_T3_NW_NW.txt", sep = "/"), sep = ",", header = TRUE)

## get subjects per task and task names, prepare column names
Subjects_T1 <- unique(latencies_T1_G_HF$vpNames)
Subjects_T2 <- unique(latencies_T2_E_HF$vpNames)
Subjects_T3 <- unique(latencies_T3_SW_E_E$vpNames)

Tasks <- unique(LDT$which_task)

columns <- c("subjectID", "gender", "age", "english_score",
             "task", "word_nonword","frequency","condition",
             "zr", "a", "t0", "v",
             "peakLat","onset", "offset", "areaLat","mean","peakAmp", "area", "width","baseline", "meanAmp")

## process task 1
Task_1 <- matrix(nrow = 0, ncol = 22, byrow = FALSE)
colnames(Task_1) <- columns

for(i in Subjects_T1){
  
  psychometrics_t1 = psychometrics %>% 
    filter(id == i) %>% 
    select(id, gender, age, english_score)
  
  # HF
  DDM_parms_t1 = DDM_parms_T1_G_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_HF)
  
  latencies_t1 = latencies_T1_G_HF %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "German"
  wordtype <- "word"
  freq <- "HF"
  cond <- NA
  
  table <- cbind(psychometrics_t1, task, wordtype, freq, cond, DDM_parms_t1, latencies_t1)
  colnames(table) <- columns
  
  Task_1 <- rbind(Task_1, table)
  
  # LF
  DDM_parms_t1 = DDM_parms_T1_G_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_LF)
  
  latencies_t1 = latencies_T1_G_LF %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "German"
  wordtype <- "word"
  freq <- "LF"
  cond <- NA
  
  table <- cbind(psychometrics_t1, task, wordtype, freq, cond, DDM_parms_t1, latencies_t1)
  colnames(table) <- columns
  
  Task_1 <- rbind(Task_1, table, deparse.level = 0)
  
  # NW
  DDM_parms_t1 = DDM_parms_T1_G_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_NW)
  
  latencies_t1 = latencies_T1_G_NW %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "German"
  wordtype <- "nonword"
  freq <- "PW"
  cond <- NA
  
  table <- cbind(psychometrics_t1, task, wordtype, freq, cond, DDM_parms_t1, latencies_t1)
  colnames(table) <- columns
  
  Task_1 <- rbind(Task_1, table, deparse.level = 0)
  
}

## process task 2
Task_2 <- matrix(nrow = 0, ncol = 22, byrow = FALSE)
colnames(Task_2) <- columns

for(i in Subjects_T2){
  
  psychometrics_t2 = psychometrics %>% 
    filter(id == i) %>% 
    select(id, gender, age, english_score)
  
  # HF
  DDM_parms_t2 = DDM_parms_T2_E_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_HF)
  
  latencies_t2 = latencies_T2_E_HF %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "English"
  wordtype <- "word"
  freq <- "HF"
  cond <- NA
  
  table <- cbind(psychometrics_t2, task, wordtype, freq, cond, DDM_parms_t2, latencies_t2)
  colnames(table) <- columns
  
  Task_2 <- rbind(Task_2, table)
  
  # LF
  DDM_parms_t2 = DDM_parms_T2_E_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_LF)
  
  latencies_t2 = latencies_T2_E_LF %>% 
    filter(vpNames == i) %>%
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "English"
  wordtype <- "word"
  freq <- "LF"
  cond <- NA
  
  table <- cbind(psychometrics_t2, task, wordtype, freq, cond, DDM_parms_t2, latencies_t2)
  colnames(table) <- columns
  
  Task_2 <- rbind(Task_2, table, deparse.level = 0)
  
  # NW
  DDM_parms_t2 = DDM_parms_T2_E_hf_lf %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_NW)
  
  latencies_t2 = latencies_T2_E_NW %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "English"
  wordtype <- "nonword"
  freq <- "PW"
  cond <- NA
  
  table <- cbind(psychometrics_t2, task, wordtype, freq, cond, DDM_parms_t2, latencies_t2)
  colnames(table) <- columns
  
  Task_2 <- rbind(Task_2, table, deparse.level = 0)
}


## process task 3
Task_3 <- matrix(nrow = 0, ncol = 22, byrow = FALSE)
colnames(Task_3) <- columns

for(i in Subjects_T3){
  
  psychometrics_t3 = psychometrics %>% 
    filter(id == i) %>% 
    select(id, gender, age, english_score)
  
  # G-G
  DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_GG)
  
  latencies_t3 = latencies_T3_SW_G_G %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "Switch"
  wordtype <- "word"
  freq <- NA
  cond <- "GG"
  
  table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  colnames(table) <- columns
  
  Task_3 <- rbind(Task_3, table)
  
  # E-E
  DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_EE)
  
  latencies_t3 = latencies_T3_SW_E_E %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "Switch"
  wordtype <- "word"
  freq <- NA
  cond <- "EE"
  
  table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  colnames(table) <- columns
  
  Task_3 <- rbind(Task_3, table)
  
  # G-E
  DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_GE)
  
  latencies_t3 = latencies_T3_SW_G_E %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "Switch"
  wordtype <- "word"
  freq <- NA
  cond <- "GE"
  
  table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  colnames(table) <- columns
  
  Task_3 <- rbind(Task_3, table)
  
  # E-G
  DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
    filter(dataset == i) %>% 
    select(zr, a, t0, v_EG)
  
  latencies_t3 = latencies_T3_SW_E_G %>% 
    filter(vpNames == i) %>% 
    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  task <- "Switch"
  wordtype <- "word"
  freq <- NA
  cond <- "EG"
  
  table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  colnames(table) <- columns
  
  Task_3 <- rbind(Task_3, table)
  
  # NW-G
#  DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
  #    filter(dataset == i) %>% 
  #    select(zr, a, t0, v_nsw_p_d)
  
  #latencies_t3 = latencies_T3_SW_NW_G %>% 
  #  filter(vpNames == i) %>% 
  #  select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  #task <- "Switch"
  #wordtype <- "word"
  #freq <- NA
  #cond <- "NWG"
  
  #table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  #colnames(table) <- columns
  
  #Task_3 <- rbind(Task_3, table)
  
  # NW-E
  #DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
  #  filter(dataset == i) %>% 
  #  select(zr, a, t0, v_nsw_p_e)
  
  #latencies_t3 = latencies_T3_SW_NW_E %>% 
  #  filter(vpNames == i) %>% 
  #  select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  #task <- "Switch"
  #wordtype <- "word"
  #freq <- NA
  #cond <- "NWE"
  
  #table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  #colnames(table) <- columns
  
  #Task_3 <- rbind(Task_3, table)
  
  # E-NW
  #DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
  #  filter(dataset == i) %>% 
  #  select(zr, a, t0, v_nsw_e_p)
  
  #  latencies_t3 = latencies_T3_SW_E_NW %>% 
  #  filter(vpNames == i) %>% 
  #    select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  #task <- "Switch"
  #wordtype <- "nonword"
  #freq <- NA
  #cond <- "ENW"
  
  #  table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  #  colnames(table) <- columns
  
  #  Task_3 <- rbind(Task_3, table)
  
  # G-NW
  #DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
  #  filter(dataset == i) %>% 
  #  select(zr, a, t0, v_nsw_d_p)
  #
  #latencies_t3 = latencies_T3_SW_G_NW %>% 
  #  filter(vpNames == i) %>% 
  #  select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  #task <- "Switch"
  #wordtype <- "nonword"
  #freq <- NA
  #cond <- "GNW"
  
  #table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  #colnames(table) <- columns
  
  #Task_3 <- rbind(Task_3, table)
  
  # NW-NW
  #DDM_parms_t3 = DDM_parms_T3_SW_Condition %>% 
  #  filter(dataset == i) %>% 
  #  select(zr, a, t0, v_nsw_p_p)
  
  #latencies_t3 = latencies_T3_SW_NW_NW %>% 
  #  filter(vpNames == i) %>% 
  #  select(peakLat,onset, offset, areaLat,mean,peakAmp, area, width,baseline,meanAmp)
  
  #task <- "Switch"
  #wordtype <- "nonword"
  #freq <- NA
  #cond <- "NWNW"
  
  #table <- cbind(psychometrics_t3, task, wordtype, freq, cond, DDM_parms_t3, latencies_t3)
  #colnames(table) <- columns
  
  #Task_3 <- rbind(Task_3, table)
  
}

## combine all tasks
all_tasks <- rbind(Task_1, Task_2, Task_3)

## recalculate width based on on- and offsets
for (i in 1:nrow(all_tasks)){
  all_tasks$width[i] <- all_tasks$offset[i] - all_tasks$onset[i]
}

## add accuracy and meanRT values to all_tasks
for (i in 1:nrow(all_tasks)){
  
  if(all_tasks$task[i] != "Switch"){
    subject = all_tasks$subjectID[i]
    task = all_tasks$task[i]
    freq = all_tasks$frequency[i]
    
    vp_data <- LDT_clean %>%
      filter(subjectID == subject, which_task == task, frequency == freq)
    n_right = sum(vp_data$correct)
    mean_rt = mean(vp_data$RT)
    
    all_tasks[i,23] = (n_right/nrow(vp_data))
    all_tasks[i,24] = mean_rt
    
  }
  else{
    subject = all_tasks$subjectID[i]
    task = all_tasks$task[i]
    cond = all_tasks$condition[i]
    
    vp_data <- LDT_clean %>%
      filter(subjectID == subject, which_task == task, condition == cond)
    n_right = sum(vp_data$correct)
    mean_rt = mean(vp_data$RT)
    
    all_tasks[i,23] = (n_right/nrow(vp_data))
    all_tasks[i,24] = mean_rt
  }
}

## add column names
colnames(all_tasks)[23] <- c("accuracy")
colnames(all_tasks)[24] <- c("meanRT")

## define condition as factor
all_tasks$condition <- as.factor(all_tasks$condition)

all_tasks <- all_tasks
## save data
save(all_tasks, file = file.path(saveFolder,"all_Tasks.RData"))
#write.csv(all_tasks,file.path("data/all_tasks.csv"), row.names = FALSE)

## clean workspace
remove(DDM_parms_T1_G_hf_lf, DDM_parms_T2_E_hf_lf, DDM_parms_T3_SW_Condition, DDM_parms_t1, DDM_parms_t2, DDM_parms_t3, latencies_t1, latencies_t2, latencies_T1_G_HF, latencies_T1_G_LF, latencies_T1_G_NW, latencies_T2_E_HF, latencies_T2_E_LF, latencies_T2_E_NW, latencies_t3, latencies_T3_SW_E_E, latencies_T3_SW_E_G, latencies_T3_SW_G_E, latencies_T3_SW_G_G, psychometrics_t1, psychometrics_t2, psychometrics_t3, table, columns, liesefeldFolder, i, DDMFolder, Subjects_T1, Subjects_T2, Subjects_T3, Tasks, task, wordtype, freq, cond, Task_1, Task_2, Task_3, vp_data, mean_rt, n_right, subject)
