##########################################################
##                    STATISTICS                        ##
##########################################################
## Description :: calculate statistics for DDM drift rate
##                and EEG area latency and mean amplitude 
##                for each task
## Input :::::::: all_tasks (06_final_table.R)
## Libraries :::: tidyverse, lme4, gtsummary, png
## Output ::::::: plots as png in figureFolder
##########################################################

## libraries
if(!"lme4" %in% installed.packages()[ ,"Package"]) {
  install.packages("lme4")
}

library(tidyverse)
library(lme4)
library(gtsummary)
library(png)

#### task 1 and 2 #### 
# select variables and calculate mean age
task_1_2_stats <- all_tasks %>%
  filter(task == "German" | task == "English" ) %>%
  filter(!is.na(english_score)) %>% 
  select(subjectID, task, age, english_score, frequency, a, v, areaLat, meanAmp) %>%
  mutate(centered_age = age-mean(age))
task_1_2_stats$v <- abs(task_1_2_stats$v)

## drift rate

drift_null <- lmer(v ~ (1|subjectID), data= task_1_2_stats)

drift_1 <- lmer(v ~ (1|subjectID) + (1|subjectID:frequency), data= task_1_2_stats)


anova(drift_null,drift_1)

drift_2 <- lmer(v ~ (1|subjectID) + (1|subjectID:frequency) + (1|subjectID:task), data= task_1_2_stats)
anova(drift_1,drift_2)


drift_3 <- lmer(v ~ task + (1|subjectID) + (1|subjectID:frequency)+ (1|subjectID:task), data= task_1_2_stats)
anova(drift_2,drift_3)
summary(drift_3)

drift_4 <- lmer(v ~ task + frequency + (1|subjectID) + (1|subjectID:frequency)+ (1|subjectID:task), data= task_1_2_stats)
anova(drift_3,drift_4)
summary(drift_4)

drift_5 <- lmer(v ~ task * frequency + (1|subjectID) + (1|subjectID:frequency)+ (1|subjectID:task), data= task_1_2_stats)
anova(drift_4,drift_5)



# fullest model with significant result
summary(drift_5)

# create table
t1 <- tbl_regression(drift_5, 
                     exponentiate = FALSE,
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                     label = list(task ~ "Block", frequency ~ "Frequency"),
                     intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## area latency

areaLat_null <- lmer(areaLat ~ (1|subjectID), data= task_1_2_stats)

areaLat_1 <- lmer(areaLat ~ (1|subjectID) + (1|subjectID:frequency) , data= task_1_2_stats)

anova(areaLat_null,areaLat_1)

areaLat_2 <- lmer(areaLat ~ (1|subjectID) +(1|subjectID:task), data= task_1_2_stats)

anova(areaLat_null,areaLat_2)

areaLat_3 <- lmer(areaLat ~ task + (1|subjectID), data= task_1_2_stats)

anova(areaLat_null,areaLat_3)
summary(areaLat_3)


areaLat_4 <- lmer(areaLat ~ task + frequency + (1|subjectID), data= task_1_2_stats)

anova(areaLat_3,areaLat_4)


# fullest model with significant result
summary(areaLat_3)


# create table
t2 <- tbl_regression(areaLat_3, 
                     exponentiate = FALSE,
                     label = list(task ~ "Block"),
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                     intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## mean amplitude
meanAmp_null <- lmer(meanAmp ~  (1|subjectID), data= task_1_2_stats)

meanAmp_1 <- lmer(meanAmp ~  (1|subjectID) + (1|subjectID:frequency), data= task_1_2_stats)

anova(meanAmp_null,meanAmp_1)

meanAmp_2 <- lmer(meanAmp ~  (1|subjectID) + (1|subjectID:task), data= task_1_2_stats)
anova(meanAmp_null,meanAmp_2)

meanAmp_3 <- lmer(meanAmp ~ task + (1|subjectID), data= task_1_2_stats)
anova(meanAmp_null,meanAmp_3)

meanAmp_4 <- lmer(meanAmp ~ frequency + (1|subjectID), data= task_1_2_stats)
anova(meanAmp_null,meanAmp_4)


# fullest model with significant result
summary(meanAmp_4)

# create table
t3 <- tbl_regression(meanAmp_4, 
                     exponentiate = FALSE,
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                     label = list(frequency ~ "Frequency"),
                     intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## merge tables to final table for task 1 and 2
table_task_1_2 <- tbl_merge(
  tbls = list(t1, t2, t3),
  tab_spanner = c("**Drift rate**", "**Area latency**", "**Mean amplitude**")) %>%
  modify_caption("**LMM Fixed Effects of the German and English word frequency task**")


## save table
as_gt(table_task_1_2) %>%
  gt::tab_source_note(gt::md("*Intercept values for mean age; HF = high-frequency words, LF = low-frequency words, PW = pseudowords*")) %>%
  gt::gtsave(filename = file.path(figureFolder, "Statistics_task_1_2.png"))


#### task 3 ####
# select variables
task_3_stats <- all_tasks %>%
  filter(task == "Switch") %>%
  filter(!is.na(english_score)) %>% 
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(subjectID, age, english_score, condition, a, v,t0, areaLat, meanAmp) %>% 
  mutate(centered_age = age-mean(age))

# calculate factors switch and language
for(i in 1:nrow(task_3_stats)){
  if(task_3_stats$condition[i] == "GG" | task_3_stats$condition[i] == "EE"){
    task_3_stats$switch[i] = "nonswitch"
  }
  else{
    task_3_stats$switch[i] = "switch"
  }
}
for(i in 1:nrow(task_3_stats)){
if(task_3_stats$condition[i] == "GG" | task_3_stats$condition[i] == "EG"){
  task_3_stats$language[i] = "German"
}
else{
  task_3_stats$language[i] = "English"
}  
}

# set as factors
task_3_stats$switch <- factor(task_3_stats$switch, levels = c("nonswitch", "switch"))
task_3_stats$language <- factor(task_3_stats$language, levels = c("German", "English"))



## drift rate

drift3_null <- lmer(v ~ (1|subjectID), data = task_3_stats)
drift3_1 <- lmer(v ~ (1|subjectID)+ (1|subjectID:switch), data = task_3_stats)

anova(drift3_null,drift3_1)

drift3_null <- lmer(v ~ (1|subjectID), data = task_3_stats)
drift3_1_2 <- lmer(v ~ (1|subjectID)+ (1|subjectID:language) , data = task_3_stats)


anova(drift3_null, drift3_1_2)

drift3_1_2 <- lmer(v ~ (1|subjectID)+ (1|subjectID:language), data = task_3_stats)
drift3_2 <- lmer(v ~ language + (1|subjectID)+ (1|subjectID:language), data = task_3_stats)

anova(drift3_1_2, drift3_2)
summary(drift3_2)

drift3_2 <- lmer(v ~ language + (1|subjectID)+ (1|subjectID:language), data = task_3_stats)
drift3_3 <- lmer(v ~ language + switch + (1|subjectID)+ (1|subjectID:language), data = task_3_stats)

anova(drift3_2, drift3_3)

summary(drift3_3)

drift3_3 <- lmer(v ~ language + switch + (1|subjectID)+ (1|subjectID:language), data = task_3_stats)
drift3_4 <- lmer(v ~ language * switch + (1|subjectID)+ (1|subjectID:language), data = task_3_stats)

anova(drift3_3, drift3_4)


# create table
task3_t1 <- tbl_regression(drift3_3, 
                     exponentiate = FALSE,
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                     label = list(language ~ "Language", switch ~ "Switch"),
                     intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## area latency

areaLat3_null <- lmer(areaLat ~ (1|subjectID), data = task_3_stats)
areaLat3_1 <- lmer(areaLat ~ (1|subjectID) + (1|subjectID:language), data = task_3_stats)
areaLat3_1_2 <- lmer(areaLat ~ (1|subjectID) + (1|subjectID:switch), data = task_3_stats)

anova(areaLat3_null, areaLat3_1_2)

areaLat3_2 <- lmer(areaLat ~ language + (1|subjectID), data = task_3_stats)

anova(areaLat3_null, areaLat3_2)
summary(areaLat3_2)

areaLat3_3 <- lmer(areaLat ~ switch + (1|subjectID), data = task_3_stats)

anova(areaLat3_null, areaLat3_3)
summary(areaLat3_3)


# create table
task3_t2 <- tbl_regression(areaLat3_null, 
                     exponentiate = FALSE,
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                     intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**") 


## mean amplitude

meanAmp3_null <- lmer(meanAmp ~ (1|subjectID), data= task_3_stats)

meanAmp3_1 <- lmer(meanAmp ~  (1|subjectID) + (1|subjectID:language), data= task_3_stats)

anova(meanAmp3_null, meanAmp3_1)

meanAmp3_2 <- lmer(meanAmp ~  (1|subjectID) + (1|subjectID:language)+(1|subjectID:switch), data= task_3_stats)

anova(meanAmp3_1, meanAmp3_2)

meanAmp3_3 <- lmer(meanAmp ~  language + (1|subjectID) + (1|subjectID:language) + (1|subjectID:switch), data= task_3_stats)

anova(meanAmp3_2, meanAmp3_3)
summary(meanAmp3_3)

meanAmp3_4 <- lmer(meanAmp ~  language + switch + (1|subjectID) + (1|subjectID:language)+ (1|subjectID:switch), data= task_3_stats)

anova(meanAmp3_3, meanAmp3_4)
summary(meanAmp3_4)

meanAmp3_5 <- lmer(meanAmp ~  language * switch + (1|subjectID) + (1|subjectID:language)+ (1|subjectID:switch), data= task_3_stats)

anova(meanAmp3_4, meanAmp3_5)
summary(meanAmp3_4)

# fullest model with significant result
summary(meanAmp3_3)


task3_t3 <- tbl_regression(meanAmp3_3,
                           label = list(language ~ "Language"),
                     exponentiate = FALSE,
                     pvalue_fun = ~style_pvalue(.x, digits = 2),
                                         intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## merge tables to final table for task 1 and 2
table_task_3 <- tbl_merge(
  tbls = list(task3_t1, task3_t2, task3_t3),
  tab_spanner = c("**Drift rate**", "**Area latency**", "**Mean amplitude**")
) %>%
  modify_caption("**LMM Fixed Effects of the switch task**")


## save table
as_gt(table_task_3) %>%
  gt::tab_source_note(gt::md("*Intercept values for mean age*")) %>%
  gt::gtsave(filename = file.path(figureFolder, "Statistics_task_3.png"))



#### combining v and ERP ####
task_all_stats <- full_join(task_1_2_stats, task_3_stats)

task_all_stats <- task_all_stats %>% 
  select(subjectID,age, centered_age,english_score,a,v,areaLat,meanAmp)


##plots combined tasks
task_all_stats_gg <- task_all_stats %>% 
  select(age,english_score,a,v,areaLat,meanAmp)

filename <- file.path(figureFolder,"ggairs_task_all.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")
ggpairs(task_all_stats_gg, title = "Variable overview", cardinality_threshold = NULL) + theme_bw()
dev.off()

# plot age and english_score
  bxp_a <- ggplot(task_all_stats, aes(x = age, y = a, color = english_score)) +
    geom_point() +
    theme_bw() +
    xlab("Age") +
    ggtitle("Boundary separation a") 
    bxp_a$labels$colour <- "English score"


bxp_v <- ggplot(task_all_stats, aes(x = age, y = v, color = english_score)) +
  geom_point() +
  theme_bw() +
  xlab("Age") +
  ggtitle("Drift rate v")
bxp_v$labels$colour <- "English score"

bxp_areaLat <- ggplot(task_all_stats, aes(x = age, y = areaLat, color = english_score)) +
  geom_point() +
  theme_bw() +
  xlab("Age") +
  ylab("Area latency (s)")+
  ggtitle("N400 area latency")
bxp_areaLat$labels$colour <- "English score"

bxp_meanAmp <- ggplot(task_all_stats, aes(x = age, y = meanAmp, color = english_score)) +
  geom_point() +
  theme_bw() +
  xlab("Age") +
  ylab("Mean amplitude (µV)")+
  ggtitle("N400 mean amplitude")
bxp_meanAmp$labels$colour <- "English score"

filename <- file.path(figureFolder,"All_task_age_english_score.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")

ggarrange(bxp_a, bxp_v, bxp_areaLat, bxp_meanAmp, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

dev.off()

## stats drivt rate all tasks
comb_null_v <- lmer(v ~ (1|subjectID), data = task_all_stats)
comb_1_v <- lmer(v ~ (1|subjectID), data = task_all_stats)

anova(comb_null_v , comb_1_v)

comb_2_v <- lmer(v ~ english_score + (1|subjectID), data = task_all_stats)
anova(comb_1_v , comb_2_v)
summary(comb_2_v)

comb_3_v <- lmer(v ~ english_score + centered_age + (1|subjectID), data = task_all_stats)
anova(comb_2_v , comb_3_v)
summary(comb_3_v)

comb_4_v <- lmer(v ~  english_score + centered_age + meanAmp + (1|subjectID), data = task_all_stats)
anova(comb_3_v  , comb_4_v )
summary(comb_4_v)


comb_5_v <- lmer(v ~ english_score + centered_age + meanAmp + areaLat + (1|subjectID), data = task_all_stats)
anova(comb_4_v , comb_5_v)


comb_t1 <- tbl_regression(comb_4_v, 
                           exponentiate = FALSE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(english_score ~ "English score", centered_age ~ "Age", meanAmp ~"Mean amplitude"),
                           intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## stats area latency all tasks
comb_null_areaLat <- lmer(areaLat ~ (1|subjectID), data = task_all_stats)
comb_1_areaLat <- lmer(areaLat ~ (1|subjectID), data = task_all_stats)

anova(comb_null_areaLat , comb_1_areaLat)

comb_2_areaLat <- lmer(areaLat ~ english_score + (1|subjectID), data = task_all_stats)
anova(comb_1_areaLat , comb_2_areaLat)
summary(comb_2_areaLat)

comb_3_areaLat <- lmer(areaLat ~ english_score + centered_age + (1|subjectID), data = task_all_stats)
anova(comb_2_areaLat , comb_3_areaLat)

comb_4_areaLat <- lmer(areaLat ~ english_score + centered_age + meanAmp + (1|subjectID), data = task_all_stats)
anova(comb_3_areaLat , comb_4_areaLat)

comb_t2 <- tbl_regression(comb_2_areaLat, 
                          exponentiate = FALSE,
                          pvalue_fun = ~style_pvalue(.x, digits = 2),
                          label = list(english_score ~ "English score"),
                          intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")

## stats mean amplitude all tasks

comb_null_meanAmp <- lmer(meanAmp ~ (1|subjectID), data = task_all_stats)
comb_1_meanAmp <- lmer(meanAmp ~ english_score + (1|subjectID), data = task_all_stats)

anova(comb_null_meanAmp , comb_1_meanAmp)

comb_2_meanAmp <- lmer(meanAmp ~ centered_age + (1|subjectID), data = task_all_stats)
anova(comb_null_meanAmp , comb_2_meanAmp)

summary(comb_2_meanAmp)

comb_3_meanAmp <- lmer(meanAmp ~ centered_age + areaLat + (1|subjectID), data = task_all_stats)
anova(comb_2_meanAmp , comb_3_meanAmp)

comb_t3 <- tbl_regression(comb_2_meanAmp, 
                          exponentiate = FALSE,
                          pvalue_fun = ~style_pvalue(.x, digits = 2),
                         label = list(centered_age ~ "Age"),
                          intercept= TRUE) %>% 
  add_global_p() %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% 
  modify_header(label = "**Variable**")


## merge tables to final table for all tasks
table_task_all <- tbl_merge(
  tbls = list(comb_t1, comb_t2, comb_t3),
  tab_spanner = c("**Drift rate**", "**Area latency**", "**Mean amplitude**")) %>%
  modify_caption("**LMM Fixed Effects for all combined tasks**")

as_gt(table_task_all) %>%
  gt::tab_source_note(gt::md("*Intercept values for mean age*")) %>%
  gt::gtsave(filename = file.path(figureFolder, "Statistics_task_all.png"))



## clean workspace
remove(areaLat_1, areaLat_2, areaLat_3, areaLat_4, areaLat_null, areaLat3_null, areaLat3_1, areaLat3_1_2, areaLat3_2, areaLat3_3,
       comb_1_areaLat, comb_1_meanAmp, comb_1_v, comb_2_areaLat, comb_2_meanAmp, comb_2_v, comb_3_areaLat, comb_3_meanAmp,
       comb_3_v, comb_4_areaLat, comb_4_v, comb_5_v, comb_null_areaLat, comb_null_meanAmp, comb_null_v, comb_t1, comb_t2,
       comb_t3, drift_1, drift_2, drift_3, drift_4, drift_5, drift_null, drift3_1, drift3_1_2, drift3_2, drift3_3, drift3_4,
       drift3_null, meanAmp_1, meanAmp_2, meanAmp_3, meanAmp_4, meanAmp_null, meanAmp3_1, meanAmp3_2, meanAmp3_3, meanAmp3_4,
       meanAmp3_5, meanAmp3_null, t1, t2, t3, table_task_1_2, table_task_3, table_task_all, task_1_2_stats, task_3_stats,
       task_all_stats, task_all_stats_gg, task3_t1, task3_t2, task3_t3, bxp_a, bxp_areaLat, bxp_meanAmp, bxp_v)
