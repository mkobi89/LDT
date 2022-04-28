##########################################################
##               SUMMARISING TABLES                     ##
##########################################################
## Description :: creates gtsummary tables for psychometrics
##                stimulus material, RT, accuracy,
##                DDM and EEG parameters for each task
## Input :::::::: psychometrics (02_psychometrics.R)
##                all_tasks (06_final_table.R)
##                LDT_clean (03_RT_outliers.R)
## Libraries :::: tidyverse, gtsummary, webshot
## Output ::::::: plots as png in figureFolder
##########################################################

## libraries
if(!"gtsummary" %in% installed.packages()[ ,"Package"]) {
  install.packages("gtsummary")
}

if(!"webshot" %in% installed.packages()[ ,"Package"]) {
  install.packages("webshot")
}

library(tidyverse)
library(gtsummary)
library(webshot)

#### psychometrics ####
# select variables
psychometrics_selected <- psychometrics %>% 
  select(age, gender, english_score, HAWIE_t_value)

# get handedness
for(i in 1:nrow(psychometrics_selected)){
  if(psychometrics$Anz_re_Annett[i] >= 7) 
    {psychometrics_selected$handedness[i] = "right-handed"}
  else if(psychometrics$Anz_re_Annett[i] == 6) 
    {psychometrics_selected$handedness[i] = "ambidextrous"}
  else if(psychometrics$Anz_re_Annett[i] <= 5) 
    {psychometrics_selected$handedness[i] = "left-handed"}
}

psychometrics_selected$handedness <- factor(psychometrics_selected$handedness, levels = c("right-handed","left-handed","ambidextrous"))

# create table
psycho <- tbl_summary(
  psychometrics_selected,
  label = list(age ~ "Age", gender ~ "Gender", english_score ~ "English score", HAWIE_t_value ~ "HAWIE t score", handedness ~ "Handedness"),
  type = list(age ~ 'continuous2', english_score ~ 'continuous2', HAWIE_t_value ~ 'continuous2'),
  statistic = all_continuous() ~ c("{N_nonmiss}",
                                   "{mean}",
                                   "{median} ({p25}, {p75})", 
                                   "{min}, {max}"),) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  modify_caption("**Sample characteristics**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(psycho) %>%
  gt::gtsave(filename = file.path(figureFolder, "psychometrics.png"))

#### stimulus table ####
## data path
dataFolder   <- file.path("data/rawdata")

stimulus <- read_delim(paste(dataFolder,"stimuli.csv", sep="/"), 
                                     ";", escape_double = FALSE, trim_ws = TRUE)

stimulus$task <- as.factor(stimulus$task)
stimulus$condition <- as.factor(stimulus$condition)

stimulus <- stimulus %>% 
  select(task, condition, opm)


# german table

stimulus_ger <- stimulus %>% 
  filter(task == "German") %>% 
  select(condition,opm)
stimulus_ger$condition <- droplevels(stimulus_ger$condition)

s_ger <- tbl_summary(
  stimulus_ger,
  label = list(opm ~ "Occurances per million"),
  by = condition,
  digits = all_continuous() ~ 2,
  statistic = list(all_continuous() ~ "{mean} ({sd})"),
  missing = "no") %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  italicize_levels()

# english table

stimulus_eng <- stimulus %>% 
  filter(task == "English") %>% 
  select(condition,opm)
stimulus_eng$condition <- droplevels(stimulus_eng$condition)

s_eng <- tbl_summary(
  stimulus_eng,
  label = list(opm ~ "Occurances per million"),
  by = condition,
  digits = all_continuous() ~ 2,
  statistic = list(all_continuous() ~ "{mean} ({sd})"),
  missing = "no") %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  italicize_levels()

# switch table

stimulus_switch <- stimulus %>% 
  filter(task == "Switch") %>% 
  select(condition,opm)
stimulus_switch$condition <- droplevels(stimulus_switch$condition)

s_switch <- tbl_summary(
  stimulus_switch,
  label = list(opm ~ "Occurances per million"),
  by = condition,
  digits = all_continuous() ~ 2,
  statistic = list(all_continuous() ~ "{mean} ({sd})"),
  missing = "no") %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  italicize_levels()


## merge tables
t_stimulus <- tbl_merge(
  tbls = list(s_ger,s_eng,s_switch),
  tab_spanner = c("**German Task**","**English Task**", "**Switch Task**")) %>%
  modify_caption("**Stimulus characteristics**")


# save table to figure folder
as_gt(t_stimulus) %>%
  gt::gtsave(filename = file.path(figureFolder, "stimulus.png"))



#### RT tables for task 1 and 2 ####

# select data
RT_1_2 <-  LDT_clean %>%
  filter(which_task == "German" | which_task == "English")  %>%
  select(frequency, which_task, RT)

# create table
rt_1_2 <- tbl_strata(
  RT_1_2,
  strata = c(which_task),
  .tbl_fun =
    ~ .x %>%
    tbl_summary(by = frequency, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no"))    %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords") %>%
  modify_caption("**RT for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(rt_1_2) %>%
  gt::gtsave(filename = file.path(figureFolder, "RT_task_1_2.png"))

#### RT table task 3 ####
# select variables
RT_3 <-  LDT_clean %>%
  filter(which_task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, RT)

# drop unused levels
RT_3$condition <- droplevels(RT_3$condition)
RT_3$condition <- factor(RT_3$condition, levels = c("GG","EG","EE","GE"))

# create table
rt_3 <-    tbl_summary(RT_3,
      by = condition, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no")%>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); GG = German nonswitch, EG = German switch, EE = English nonswitch, and GE = English switch trials") %>%
  modify_caption("**RT for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()   

# save table to figure folder
as_gt(rt_3) %>%
  gt::gtsave(filename = file.path(figureFolder, "RT_task_3.png"))


#### accuracy values tables for task 1 and 2 ####

# select variables
ACC_1_2 <-  all_tasks %>%
  filter(task == "German" | task == "English")  %>%
  select(frequency, task, accuracy)

# create table
acc_1_2 <- tbl_strata(
  ACC_1_2,
  strata = task,
  .tbl_fun = 
    ~ .x %>%
    tbl_summary(
  by = frequency, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no")) %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords") %>%
  modify_caption("**Accuracy for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(acc_1_2) %>%
  gt::gtsave(filename = file.path(figureFolder, "ACC_task_1_2.png"))


#### accuracy values tables for task 3 ####

ACC_3 <-  all_tasks %>%
  filter(task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, accuracy)

ACC_3$condition <- droplevels(ACC_3$condition)
ACC_3$condition <- factor(ACC_3$condition, levels = c("GG","EG","EE","GE"))

    
acc_3 <- tbl_summary(ACC_3,
  by = condition, 
  statistic = list(all_continuous() ~ "{mean} ({sd})"),
  missing = "no") %>% 
  modify_header(label = "**Variable**", update = NULL) %>% # update the column header
  modify_footnote(update = all_stat_cols() ~ "Mean (SD); GG = German nonswitch, EG = German switch, EE = English nonswitch, and GE = English switch trials") %>%
  modify_caption("**Accuracy for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(acc_3) %>%
  gt::gtsave(filename = file.path(figureFolder, "ACC_task_3.png"))

#### Drift diffusion parameters tables for task 1 and 2 ####

# select variables
DDM_1_2_all <-  all_tasks %>%
  filter(task == "German" | task == "English")  %>%
  select(frequency, task, v, zr,a,t0)

# create table
ddm_1_2 <- tbl_strata(
  DDM_1_2_all,
  strata = c(task),
  .tbl_fun =
    ~ .x %>%
    tbl_summary(by = frequency, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no")) %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords") %>%
  modify_caption("**DDM parameters for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_1_2) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_1_2_all.png"))

# v only
DDM_1_2_v <-  all_tasks %>%
  filter(task == "German" | task == "English")  %>%
  select(frequency, task, v)

ddm_1_2_v <- tbl_strata(
  DDM_1_2_v,
  strata = c(task),
  .tbl_fun =
    ~ .x %>%
    tbl_summary(by = frequency, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no")) %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  #  modify_spanning_header() %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords") %>%
  modify_caption("**DDM parameters for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_1_2_v) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_1_2_v.png"))

# all others

DDM_1_2_others <-  all_tasks %>%
  filter(task == "German" | task == "English")  %>%
  filter(frequency == "HF") %>% 
  select(task, zr,a,t0)

DDM_1_2_others$task <- droplevels(DDM_1_2_others$task)

ddm_1_2_others <-    tbl_summary(DDM_1_2_others,
                                 by = task, 
                                 statistic = list(all_continuous() ~ "{mean} ({sd})"),
                                 missing = "no") %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords") %>%
  modify_caption("**DDM parameters for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_1_2_others) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_1_2_others.png"))

#### Drift diffusion parameters tables for task 3 ####

# select variables and drop unused levels
DDM_3 <-  all_tasks %>%
  filter(task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, v, zr,a,t0)

DDM_3$condition <- droplevels(DDM_3$condition)
DDM_3$condition <- factor(DDM_3$condition, levels = c("GG","EG","EE","GE"))

ddm_3 <- tbl_summary(DDM_3,
                     by = condition, 
                     statistic = list(all_continuous() ~ "{mean} ({sd})"),
                     missing = "no") %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); GG = German nonswitch, EG = German switch, EE = English nonswitch, and GE = English switch trials" ) %>%
  modify_caption("**DDM parameters for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_3) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_3.png"))

# v only
DDM_3_v <-  all_tasks %>%
  filter(task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, v)

DDM_3_v$condition <- droplevels(DDM_3_v$condition)
DDM_3_v$condition <- factor(DDM_3_v$condition, levels = c("GG","EG","EE","GE"))

ddm_3_v <-    tbl_summary(DDM_3_v,
                           by = condition, 
                           statistic = list(all_continuous() ~ "{mean} ({sd})"),
                           missing = "no") %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); GG = German nonswitch, EG = German switch, EE = English nonswitch, and GE = English switch trials") %>%
  modify_caption("**DDM parameters for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_3_v) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_3_v.png"))

# all others

DDM_3_others <-  all_tasks %>%
  filter(task == "Switch")  %>%
  filter(condition == "GG") %>% 
  select(zr, a, t0)

#DDM_3_others$condition <- droplevels(DDM_3_others$condition)
#DDM_3_others$condition <- factor(DDM_3_others$condition, levels = c("GG","EG","EE","GE"))

ddm_3_others <- tbl_summary(DDM_3_others,
            statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(update = all_stat_cols() ~ "Mean (SD)") %>%
  modify_caption("**DDM parameters for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(ddm_3_others) %>%
  gt::gtsave(filename = file.path(figureFolder, "DDM_task_3_others.png"))

#### EEG tables for task 1 and 2 ####

# select variables
EEG_1_2 <-  all_tasks %>%
  filter(task == "German" | task == "English")  %>%
  select(frequency, task, meanAmp, areaLat, onset, offset, width)

# create table 
eeg_1_2 <- tbl_strata(
  EEG_1_2,
  strata = c(task),
  .tbl_fun =
    ~ .x %>%
    tbl_summary(by = frequency, 
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                missing = "no")) %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); HF = high-frequency words, LF = low-frequency words, PW = pseudowords, meanAmp = N400 mean amplitude, areaLat = N400 area latency") %>%
  modify_caption("**EEG parameters for German and English wordfrequency task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(eeg_1_2) %>%
  gt::gtsave(filename = file.path(figureFolder, "EEG_task_1_2.png"))

#### EEG tables for task 3 ####
# select variables
EEG_3 <-  all_tasks %>%
  filter(task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, meanAmp, areaLat, onset, offset, width)

EEG_3$condition <- droplevels(EEG_3$condition)
EEG_3$condition <- factor(EEG_3$condition, levels = c("GG","EG","EE","GE"))

eeg_3 <- tbl_summary(EEG_3,
                     by = condition, 
                     statistic = list(all_continuous() ~ "{mean} ({sd})"),
                     missing = "no") %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  modify_footnote(
    update = all_stat_cols() ~ "Mean (SD); GG = German nonswitch, EG = German switch, EE = English nonswitch, and GE = English switch trials, meanAmp = N400 mean amplitude, areaLat = N400 area latency") %>%
  modify_caption("**EEG parameters for switch task**") %>% 
  bold_labels() %>%
  italicize_levels()

# save table to figure folder
as_gt(eeg_3) %>%
  gt::gtsave(filename = file.path(figureFolder, "EEG_task_3.png"))

#### clean workspace ####
remove(acc_1_2, ACC_1_2, acc_3, ACC_3, ddm_1_2, ddm_1_2_others, ddm_1_2_v, DDM_1_2_all, DDM_1_2_others, DDM_1_2_v, ddm_3, ddm_3_others, ddm_3_v, DDM_3, DDM_3_others, DDM_3_v, psycho, psychometrics_selected, rt_1_2, rt_3, RT_1_2, RT_3, eeg_1_2, eeg_3,  EEG_1_2, EEG_3, psycho, psychometrics_selected, s_eng, s_ger, s_switch, stimulus, stimulus_eng, stimulus_ger, stimulus_switch, t_stimulus,i)
