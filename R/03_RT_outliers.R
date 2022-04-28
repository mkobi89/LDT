##########################################################
##                PROCESS RT DATA                       ##
##########################################################
## Description :: Removes outlier trials from LDT 
##                dataframe based on trim factor, 
##                summarizes rt data per subject 
##                relevant data
## Input :::::::: LDT dataframe (01_answer_files_to_LDT.R)
## Libraries :::: dyplr, ggplot2
## Output ::::::: LDT_clean.Rdata, descriptives.csv 
##                (per subject), exclude.csv 
##########################################################

## libraries
if (!"tidyverse" %in% installed.packages()[, "Package"]) {
  install.packages("tidyverse")
}

library(tidyverse)

## data path
dataFolder     <- file.path("data")
saveFolder_out <- file.path("data/outliers")

## set trimm factor for removing trials: mean(rt) + trim * sd(rt)
trim <- 3

## get subject and task names
Subjects <- unique(LDT$subjectID)
Tasks <- unique(LDT$which_task)

## preallocate dataframes
exclude <- matrix(nrow = 0, ncol = 18, byrow = FALSE, dimnames = NULL)
descriptives <- matrix(nrow = 0, ncol = 7, byrow = FALSE, dimnames = NULL)

## process subjects --> recognise trim boundaries per subject and task and get trial numbers to exclude from analysis, save it in exlude dataframe. Get RT descriptives with and without trim per subject
for(i in 1:nlevels(Subjects)){

  ## task 1
  
  ldt_t1 = LDT %>% 
    filter(subjectID == Subjects[i], which_task == Tasks[1])
  
  ldt_t1 = ldt_t1 %>%   
    mutate(TrialNumber = 1:nrow(ldt_t1))
  
  #hist(ldt_t1$RT, main = "Response Times German Task", sub= Subjects[i], ylab = "Frequency", xlab = "RT", breaks=50, col = rgb(1, 0, 0, 0.6),
  #     xlim = c(0,10), ylim = c(0, 100))
  
  cutoffs = ldt_t1 %>%
    summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
    mutate(upper = mean+trim*sd)
  
  lower_trim = ldt_t1 %>%
    filter(RT < 0.300)
  
  upper_trim = ldt_t1 %>%
    filter(RT > 2)
  
  trimmedData = ldt_t1 %>%
    filter(RT < 2, RT > 0.300)
  
  trimmedData_cutoffs = trimmedData %>% 
    summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
    mutate(upper = mean+trim*sd)
  
  exclude = rbind(exclude, lower_trim,upper_trim)
  
  cutoffs_summary = rbind(cutoffs, trimmedData_cutoffs)
  
  descript <- data.frame(subjectID = rep(Subjects[i], 2),
                         which_task = c("German", "German"), 
                         Trim = c("No", "Yes"), 
                         cutoffs_summary)
  
  descriptives= rbind(descriptives, descript)
  

  
  ## task 2
  
  ldt_t2 = LDT %>% 
    filter(subjectID == Subjects[i], which_task == Tasks[2])
  
  ldt_t2 = ldt_t2 %>%   
    mutate(TrialNumber = 1:nrow(ldt_t2))
  
  #hist(ldt_t2$RT, main = "Response Times English Task", ylab = "Frequency", xlab = "RT", breaks=50, col = rgb(0, 0, 1, 0.2), add = T)
  #legend("topright", legend=c("German", "English"), col=c(rgb(1, 0, 0, 1), rgb(0, 0, 1, 1)), lty = 1)
  
  cutoffs = ldt_t2 %>%
    summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
    mutate(upper = mean+trim*sd)
  
  lower_trim = ldt_t2 %>%
    filter(RT < 0.300)
  
  upper_trim = ldt_t2 %>%
    filter(RT > 2)
  
  trimmedData = ldt_t2 %>%
    filter(RT < 2, RT > 0.300)
  
  trimmedData_cutoffs = trimmedData %>% 
    summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
    mutate(upper = mean+trim*sd)
  
  exclude = rbind(exclude, lower_trim,upper_trim)
  
  cutoffs_summary = rbind(cutoffs, trimmedData_cutoffs)
  
  descript <- data.frame(subjectID = rep(Subjects[i], 2), which_task = c("English", "English"), Trim = c("No", "Yes"), cutoffs_summary)
  
  descriptives= rbind(descriptives, descript)

  ## task 3
  
  ldt_t3 = LDT %>% 
    filter(subjectID == Subjects[i], which_task == Tasks[3])
  
  if(i != 12){
  
    ldt_t3 = ldt_t3 %>%   
      mutate(TrialNumber = 1:nrow(ldt_t3))
    
    #hist(ldt_t3$RT, main = "Response Times Sswitch Task", ylab = "Frequency", xlab = "RT", breaks=50, col = rgb(0, 1, 0, 0.2), add = T)
    #legend("topright", legend=c("German", "English", "Switch"), col=c(rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 1, 0, 1)), lty = 1)
    
    cutoffs = ldt_t3 %>%
      summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
      mutate(upper = mean+trim*sd)
    
    lower_trim = ldt_t3 %>%
      filter(RT < 0.300)
    
    upper_trim = ldt_t3 %>%
      filter(RT > 2)
    
    trimmedData = ldt_t3 %>%
      filter(RT < 2, RT > 0.300)
    
    trimmedData_cutoffs = trimmedData %>% 
      summarise(mean = mean(RT, na.rm = TRUE), sd = sd(RT, na.rm = TRUE), median = median(RT, na.rm = TRUE)) %>%
      mutate(upper = mean+trim*sd)
    
    exclude = rbind(exclude, lower_trim,upper_trim)
    
    cutoffs_summary = rbind(cutoffs, trimmedData_cutoffs)
    
    descript <- data.frame(subjectID = rep(Subjects[i], 2), which_task = c("Switch", "Switch"), Trim = c("No", "Yes"), cutoffs_summary)
    
    descriptives= rbind(descriptives, descript)}
}

## create LDT_clean dataframe without outliers based on exclude variable
Subjects <- unique(LDT$subjectID)
LDT_clean <- matrix(nrow = 0, ncol = 18, byrow = FALSE)

## task 1
for(i in Subjects){
  # Get data from a single subject with variables of interest
  vp_data <- LDT %>%
    filter(subjectID == i, which_task == "German")
  
  # create table with outliers
  outliers <- exclude %>%
    filter(subjectID == i, which_task == "German") %>%
    select(subjectID, which_task, TrialNumber)
  
  # remove trials from dm_data set
  if (nrow(outliers) >= 1) {
    vp_data <- vp_data[-c(outliers$TrialNumber),]
  }
  
  LDT_clean <- rbind(LDT_clean, vp_data)
}

## task 2
for(i in Subjects){
  # Get data from a single subject with variables of interest
  vp_data <- LDT %>%
    filter(subjectID == i, which_task == "English")
  
  # create table with outliers
  outliers <- exclude %>%
    filter(subjectID == i, which_task == "English") %>%
    select(subjectID, which_task, TrialNumber)
  
  # remove trials from dm_data set
  if (nrow(outliers) >= 1) {
    vp_data <- vp_data[-c(outliers$TrialNumber),]
  }
  
  LDT_clean <- rbind(LDT_clean, vp_data)
}

## task 3
for(i in Subjects){
  # Get data from a single subject with variables of interest
  vp_data <- LDT %>%
    filter(subjectID == i, which_task == "Switch")
  
  # create table with outliers
  outliers <- exclude %>%
    filter(subjectID == i, which_task == "Switch") %>%
    select(subjectID, which_task, TrialNumber)
  
  # remove trials from dm_data set
  if (nrow(outliers) >= 1) {
    vp_data <- vp_data[-c(outliers$TrialNumber),]
  }
  
  LDT_clean <- rbind(LDT_clean, vp_data)
}

## add age and english score data from psychometrics dataframe to LDT_clean
LDT_clean_psycho <- data.frame(matrix(ncol=20,nrow=0))
Subjects <- unique(psychometrics$id)


for (i in Subjects){
  vp_data <- LDT_clean %>%
    filter(subjectID == i)
  
  psycho <- psychometrics %>%
    filter(id == i) %>%
    select(age, english_score)
  
  age <- rep(psycho$age, times=nrow(vp_data))
  english_score <- rep(psycho$english_score, times=nrow(vp_data))

  vp_data <- cbind(vp_data, age, english_score)
  
  LDT_clean_psycho <- rbind(LDT_clean_psycho, vp_data)
}

LDT_clean <- LDT_clean_psycho


## save data
save(LDT_clean, file = paste(dataFolder, "LDT_clean.RData", sep = "/"))
#write.csv(LDT_clean, file.path(dataFolder,'LDT_clean.csv'))

write.csv(descriptives, file.path(saveFolder_out,'RT_descriptives.csv'))
write.csv(exclude, file.path(saveFolder_out,'RT_exclude_trials.csv'))

## clean workspace
remove(cutoffs, cutoffs_summary, descript, ldt_t1, ldt_t2, ldt_t3,
       lower_trim, upper_trim, trimmedData, trimmedData_cutoffs, i, Subjects, 
       Tasks, trim, dataFolder, saveFolder_out, outliers, vp_data, LDT_clean_psycho,
       age, english_score, psycho)


