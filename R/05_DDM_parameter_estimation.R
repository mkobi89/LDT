##########################################################
##     ESTIMATING DDM PARAMETERS WITH FAST_DM           ##
##########################################################
## Description :: define desired parameters for the ddm
##                for estimating ddm parameters per task
## Input :::::::: .dat dataset for each person 
##                (04_DDM_generate_datasets_without_outliers)
## Libraries :::: uses fast-dm.exe in dm_datasets - folder
## Output ::::::: .ctl file for each model with configurations
##                .txt file with the estimated parameters
##                for each task
##########################################################

## data path
dataFolder <- file.path("data/dm_datasets")

## Each dataset contains 5 Variables
# subjectID
# which_task: levels "german", "english", "switch"
# word_nonword: levels "word", "nonword"
# word_pseudoword_numeric: levels "1", "0"
# frequency: levels "HF", "LF", "NW" for German and English task
# condition: 9 switch and nonswitch conditions for switch task
# correct response: levels "1", "0"
# RT: Reaction time in seconds

## set directroy to .dat files
setwd(dataFolder)


## Create Fast-DM Control-File for Parameter-Estimation for task 1
# everything set to 0
# zr not set, because it will be estimated and varying for every participant
# v depends on Wordfrquency/Condition
# RESPONSE is word_pseudoword_numeric
cat(c(
  "precision ", 3, "\n", #Precision of Parameter Estimation 
  "method ", "ml", "\n", #Optimierungskriterium (Maximum Likelihood)
  "set szr 0", "\n", #Fixation of szr
  "set sv 0", "\n", #Fixation of sv
  "set st0 0.2", "\n", #Variation of t0
  "set d 0", "\n", #Fixation of d
  "set p 0", "\n", #Fixation of p
  "depends v frequency", "\n", #This command lets drift rates vary over Wordfrequency
  "format subjectID which_task word_nonword RESPONSE frequency correct TIME", "\n",
  "load *German_final.dat", "\n",
  "log parms_G_hf_lf_final.txt"), 
  file=paste("parms_G_hf_lf_final.ctl",sep=""), sep="", fill=F)

##Executing Fast-DM - depending on operating system
unlink("parms_G_hf_lf_final.txt")
if(Sys.info()[['sysname']] == "Darwin"){
  system("./fast-dm parms_G_hf_lf_final.ctl")
}else if(Sys.info()[['sysname']] == "Windows"){
  system("fast-dm.exe parms_G_hf_lf_final.ctl")
}


## Create Fast-DM Control-File for Parameter-Estimation for task 2
# everything set to 0
# z not set, because it will be estimated and varying for every participant
# v depends on Wordfrquency/Condition
# RESPONSE is word_pseudoword_numeric
cat(c(
  "precision ", 3, "\n", #Precision of Parameter Estimation 
  "method ", "ml", "\n", #Optimierungskriterium (Maximum Likelihood)
  "set szr 0", "\n", #Fixation of szr
  "set sv 0", "\n", #Fixation of sv
  "set st0 0.2", "\n", #Variation of t0
  "set d 0", "\n", #Fixation of d
  "set p 0", "\n", #Fixation of p
  "depends v frequency", "\n", #This command lets drift rates vary over Wordfrequency
  "format subjectID which_task word_nonword RESPONSE frequency correct TIME", "\n",
  "load *English_final.dat", "\n",
  "log parms_E_hf_lf_final.txt"), 
  file=paste("parms_E_hf_lf_final.ctl",sep=""), sep="", fill=F)

##Executing Fast-DM - depending on operating system
unlink("parms_E_hf_lf_final.txt")
if(Sys.info()[['sysname']] == "Darwin"){
  system("./fast-dm parms_E_hf_lf_final.ctl")
}else if(Sys.info()[['sysname']] == "Windows"){
  system("fast-dm.exe parms_E_hf_lf_final.ctl")
}


## Create Fast-DM Control-File for Parameter-Estimation for task 3
# everything set to 0
# z not set, because it will be estimated and varying for every participan
# v depends on Wordfrquency/Condition
# RESPONSE is word_pseudoword_numeric
cat(c(
  "precision ", 3, "\n", #Precision of Parameter Estimation 
  "method ", "ml", "\n", #Optimierungskriterium (Maximum Likelihood)
  "set szr 0", "\n", #Fixation of szr
  "set sv 0", "\n", #Fixation of sv
  "set st0 0.2", "\n", #Variation of t0
  "set d 0", "\n", #Fixation of d
  "set p 0", "\n", #Fixation of p
  "depends v condition", "\n", #This command lets drift rates vary over Wordfrequency
#  "depends t0 condition", "\n", #This command lets drift rates vary over Wordfrequency
  "format subjectID which_task word_nonword RESPONSE condition correct TIME", "\n",
  "load *Switch_final.dat", "\n",
  "log parms_SW_condition_final.txt"), 
  file=paste("parms_SW_condition_final.ctl",sep=""), sep="", fill=F)

##Executing Fast-DM - depending on operating system
unlink("parms_SW_condition_final.txt")
if(Sys.info()[['sysname']] == "Darwin"){
  system("./fast-dm parms_SW_condition_final.ctl")
}else if(Sys.info()[['sysname']] == "Windows"){
  system("fast-dm.exe parms_SW_condition_final.ctl")
}


## set working directory back to project
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")


##  clean workspace 
remove(dataFolder)
