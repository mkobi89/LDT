###########################################################
##                      PSYCHOMETRICS                    ##
##                    DATA PREPROCESSING                 ##
###########################################################
## Description :: loads and merges psychometrics raw data
## Input :::::::: csv data file 
## Libraries :::: dplyr, readr, eeptools, lubridate
## Output ::::::: psychometrics.Rdata
##########################################################

## libraries
if (!"tidyverse" %in% installed.packages()[, "Package"]) {
  install.packages("tidyverse")
}
if (!"readr" %in% installed.packages()[, "Package"]) {
  install.packages("readr")
}
if (!"eeptools" %in% installed.packages()[, "Package"]) {
  install.packages("eeptools")
}
if (!"lubridate" %in% installed.packages()[, "Package"]) {
  install.packages("lubridate")
}

library(tidyverse)
library(readr)
library(eeptools)
library(lubridate)

## data path
dataFolder   <- file.path("data/rawdata")

## read data
psychometrics <- read_delim(paste(dataFolder,"psychometrics.csv", sep="/"), 
                     ";", escape_double = FALSE, trim_ws = TRUE, skip = 0)

## calculate age of participants
psychometrics$age <- floor(time_length(interval(as.Date(psychometrics$birthdate, "%d.%m.%Y"), as.Date(psychometrics$date, "%d.%m.%Y")), "years"))

## select relevant variables
psychometrics = psychometrics %>% 
  select(id, gender, age, english_score, HAWIE_t_value, t_TMT_A_in_s, F_TMT_A, t_TMT_B_in_s, F_TMT_B, Anz_li_Annett, Anz_re_Annett)

## tell R, which variables in datasets are factors and numeric
psychometrics$HAWIE_t_value <- as.numeric(psychometrics$HAWIE_t_value)
psychometrics$F_TMT_B <- as.numeric(psychometrics$F_TMT_B)
psychometrics$id <- as.factor(psychometrics$id)
psychometrics$gender <- as.factor(psychometrics$gender)

## exchange english_score 0 with NA
psychometrics$english_score[psychometrics$english_score == 0] <- NA


## save data
save(psychometrics, file = file.path(saveFolder,"psychometrcs.RData"))
#write.csv(psychometrics,file.path("data/psychometrics.csv"), row.names = FALSE)

## clean workspace
remove(dataFolder)

