##########################################################
##                       LIBRARIES                      ##
##########################################################
# install tidyverse if not installed already
if(!"tidyverse" %in% installed.packages()[ ,"Package"]) {
  install.packages("tidyverse")
}

if(!"gtsummary" %in% installed.packages()[ ,"Package"]) {
  install.packages("gtsummary")
}

# load tidyverse
library(tidyverse)


##########################################################
##                       SETTINGS                       ##
##########################################################

# set seed for reproducible results
set.seed(1234321)

#dataFolderRaw   <- file.path("data/rawdata")
figureFolder <- file.path("figures")
RFolder      <- file.path("R")
saveFolder    <- file.path("data")


##########################################################
##                      PROCESSING                      ##
##########################################################

# 01: Prepare and preprocess datasets----------------------
# load and merge the data files
source(paste(RFolder, "01_answer_files_to_LDT.R", sep = "/"))

# load psychometrics
source(paste(RFolder, "02_psychometrics.R", sep = "/"))

# calculate outliers --> remove plot script from outliers
source(paste(RFolder, "03_RT_outliers.R", sep = "/"))


# calculate DDM parms
#source(paste(RFolder, "04_DDM_generate_datasets_without_outliers.R", sep = "/"))
#source(paste(RFolder, "05_DDM_parameter_estimation.R", sep = "/"))

# combine LDT, psychometrics, DDM and latencie values in a lfinal data format, calculate accuracy and mean rt for each condition and each person
source(paste(RFolder, "06_final_table.R", sep = "/"))


## remove unwanted  participants
exclude_subjects <- c("C50","C55", "C69", "C80", "CE5", "CE8","CI2","CJ0","CM7","CM8","CN1","CU4")

psychometrics <- psychometrics %>% 
  filter(!(id %in% exclude_subjects))
psychometrics$id <- droplevels(psychometrics$id)

all_tasks <- all_tasks %>% 
  filter(!(subjectID %in% exclude_subjects))
all_tasks$subjectID <- droplevels(all_tasks$subjectID)

LDT_clean <- LDT_clean %>% 
  filter(!(subjectID %in% exclude_subjects))
LDT_clean$subjectID <- droplevels(LDT_clean$subjectID)

LDT <- LDT %>% 
  filter(!(subjectID %in% exclude_subjects))
LDT$subjectID <- droplevels(LDT$subjectID)

exclude <- exclude %>% 
  filter(!(subjectID %in% exclude_subjects))
exclude$subjectID <- droplevels(exclude$subjectID)


# calculate how many trials were excluded
exclude_lower <- exclude %>% 
  filter(RT < 0.3)

critical_trials_switch <- LDT %>%
  filter(condition != "GNW", condition != "ENW", condition != "NWG", condition != "NWE", condition != "NWNW")

trials_wf <- LDT %>% 
  filter(which_task != "Switch")

exclude_wf <- exclude %>% 
  filter(which_task != "Switch")


exclude_switch <- exclude %>% 
  filter(condition != "GNW", condition != "ENW", condition != "NWG", condition != "NWE", condition != "NWNW")

exclude_trial <- (nrow(exclude_switch)+nrow(exclude_wf))/(nrow(critical_trials_switch)+nrow(trials_wf))

# summary tables
source(paste(RFolder, "07_summarise_tables.R", sep = "/"))

# plots
source(paste(RFolder, "08_plots.R", sep = "/"))

#stats
source(paste(RFolder, "09_statistics.R", sep = "/"))

