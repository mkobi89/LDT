##########################################################
##                  PREPARE RAWDATA                     ##
##########################################################
## Description :: loads and merges the raw data
## Input :::::::: individual data files for each 
##                Lexical Decision Task and subject
## Libraries :::: dplyr
## Output ::::::: LDT.Rdata
##########################################################


## libraries
if (!"tidyverse" %in% installed.packages()[, "Package"]) {
  install.packages("tidyverse")
}

library(tidyverse)

## data path
dataFolder   <- file.path("data/rawdata/answers")

## file names 
fileNames <- list.files(paste(dataFolder, sep="/"),pattern = "*answers_task*", full.names = TRUE)

## read files
LDT <- lapply(fileNames, read.csv, header = TRUE, sep = ",", stringsAsFactors = FALSE)

## merge list of data.frames into single data.frame
LDT <- bind_rows(LDT)

## adjust subject names
LDT$SubjectID <- gsub("C92", "C50", LDT$SubjectID)
LDT$SubjectID <- gsub("C93", "C64", LDT$SubjectID)
LDT$SubjectID <- gsub("y", "CJ5", LDT$SubjectID)

## adjust condition name
LDT$Condition <- gsub("nsw_d_d", "GG", LDT$Condition)
LDT$Condition <- gsub("sw_d_e", "GE", LDT$Condition)
LDT$Condition <- gsub("nsw_d_p", "GNW", LDT$Condition)
LDT$Condition <- gsub("sw_e_d", "EG", LDT$Condition)
LDT$Condition <- gsub("nsw_e_e", "EE", LDT$Condition)
LDT$Condition <- gsub("nsw_e_p", "ENW", LDT$Condition)
LDT$Condition <- gsub("nsw_p_d", "NWG", LDT$Condition)
LDT$Condition <- gsub("nsw_p_e", "NWE", LDT$Condition)
LDT$Condition <- gsub("nsw_p_p", "NWNW", LDT$Condition)

## remove nonsubjects
LDT <- LDT %>%
  filter(SubjectID != "")

## tell R, which variables are factors and numeric, rearrange levels of two factors
LDT$SubjectID <- as.factor(LDT$SubjectID)
LDT$KeyName <- as.factor(LDT$KeyName)
LDT$Correct <- as.numeric(LDT$Correct)
LDT$AnswerType <- as.factor(LDT$AnswerType)

LDT$WhichTask <- factor(LDT$WhichTask, levels = c("german", "english", "switch"))
levels(LDT$WhichTask) <- c("German", "English", "Switch")

LDT$word_nonword <- factor(LDT$Word_Pseudoword, levels = c(1, 2), labels = c("word", "nonword"))

LDT$Word_Pseudoword <- factor(LDT$Word_Pseudoword, levels = c("1", "0")) #, labels = c("word", "non-word"))
LDT$Word_Pseudoword[is.na(LDT$Word_Pseudoword)] <- 0

LDT$Wordfrequency <- factor(LDT$Wordfrequency, levels = c(2, 1, 3))
levels(LDT$Wordfrequency) <- c("HF", "LF", "PW")

LDT$Language <- as.factor(LDT$Language)
levels(LDT$Language) <- c("german", "english", "nonword")

LDT$Sequenz <- as.factor(LDT$Sequenz)
LDT$Condition <- as.factor(LDT$Condition)

## rename column names
colnames(LDT) <-
  c(
    "subjectID",
    "which_task",
    "order_index",
    "word",
    "key",
    "RT",
    "correct",
    "answer_type",
    "trigger",
    "word_pseudoword_numeric",
    "frequency",
    "language",
    "test",
    "sequence",
    "condition",
    "d_e_sw",
    "e_d_sw",
    "word_nonword")

## set dataframe as tibble
LDT <- LDT %>% 
  as_tibble()

## save data
save(LDT, file = paste(saveFolder, "LDT.RData", sep = "/"))
#write.csv(LDT,file.path("data/LDT.csv"), row.names = FALSE)

## clean workspace
remove(fileNames, dataFolder)
