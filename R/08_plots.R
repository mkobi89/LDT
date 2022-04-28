##########################################################
##                    PLOTS                             ##
##########################################################
## Description :: creates plots for DDM and EEG parameters 
#                 for each task
## Input :::::::: all_tasks (06_final_table.R)
## Libraries :::: tidyverse, graphics, GGally, ggpubr, png
## Output ::::::: plots as png in figureFolder
##########################################################

## libraries
if(!"GGally" %in% installed.packages()[ ,"Package"]) {
  install.packages("GGally")
} 

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

if(!"png" %in% installed.packages()[ ,"Package"]) {
  install.packages("png")
}

library(tidyverse)
library(graphics)
library(GGally)
library(ggpubr)
library(png)


#### task 1 & 2 ####
## select important variables 
task_1_2_plot <- all_tasks %>%
  filter(task == "German" | task == "English" ) %>%
  select(task, age, english_score, frequency, a, v, areaLat, meanAmp)
task_1_2_plot$v <- abs(task_1_2_plot$v)


#### plot ggpairs task 1 & 2 ####
filename <- file.path(figureFolder,"ggairs_task_1_2.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")
ggpairs(task_1_2_plot, title = "Variable overview task 1 and 2", cardinality_threshold = NULL) + theme_bw()
dev.off()


#### plot word frequency task ####

rt_wf <-  LDT_clean %>%
  filter(which_task != "Switch")  %>%
  select(which_task, frequency, RT) %>% 
  rename(Language = which_task, Frequency = frequency)


library(plotrix)
rt_wf_summary <- rt_wf %>% 
  group_by(Language, Frequency) %>% 
  summarise(mRT = mean(RT), sd = sd(RT), se = std.error(RT))


## Plot frequency new
pd <- position_dodge(0)
bxp_rt <- ggplot(rt_wf_summary,aes(x = Frequency, y = mRT, linetype = Language, group = Language)) +
  geom_point(position = pd, size=1.5) +
  geom_line() +
  geom_errorbar(aes(ymin = mRT - se, ymax = mRT + se), width = .1, position = pd) + 
  theme_bw() +
  ggtitle("Mean RT") 


bxp_a <- ggplot(task_1_2_plot, aes(x = frequency, y = a)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Frequency") + 
  ggtitle("boundary separation a") +
  facet_wrap( ~ task)

bxp_v <- ggplot(task_1_2_plot, aes(x = frequency, y = v)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Frequency") + 
  ggtitle("Drift rate") +
  facet_wrap( ~ task)

bxp_areaLat <- ggplot(task_1_2_plot, aes(x = frequency, y = areaLat)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Area latency") +
  xlab("Frequency") + 
  facet_wrap( ~ task)

bxp_meanAmp <- ggplot(task_1_2_plot, aes(x = frequency, y = meanAmp)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Mean amplitude") +
  xlab("Frequency") + 
  facet_wrap( ~ task)

filename <- file.path(figureFolder,"Task_1_2_frequency.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")

ggarrange(bxp_rt, bxp_v, bxp_areaLat, bxp_meanAmp, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

dev.off()


#### task 3 ####
## select important variables, calculate switch and language factor
task_3_plot <- all_tasks %>%
  filter(task == "Switch") %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(age, english_score, condition, a, v, t0, zr, areaLat, meanAmp, onset, offset, width, meanRT)

for(i in 1:nrow(task_3_plot)){
  if(task_3_plot$condition[i] == "GG" | task_3_plot$condition[i] == "EE"){
    task_3_plot$Switch[i] = "Nonswitch"
  }
  else{
    task_3_plot$Switch[i] = "Switch"
  }
}
for(i in 1:nrow(task_3_plot)){
  if(task_3_plot$condition[i] == "GG" | task_3_plot$condition[i] == "EG"){
    task_3_plot$Language[i] = "German"
  }
  else{
    task_3_plot$Language[i] = "English"
  }  
}

task_3_plot$Switch <- factor(task_3_plot$Switch, levels = c("Nonswitch", "Switch"))
task_3_plot$Language <- factor(task_3_plot$Language, levels = c("German", "English"))




#### plot ggpairs task 3 ####

## subsets for ggpairs
task_3_plot_important <- task_3_plot %>%
  select(age, english_score, Language, Switch, a, v, areaLat, meanAmp)

filename <- file.path(figureFolder,"ggairs_task_3.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")
ggpairs(task_3_plot_important, title = "Variable overview task 3", cardinality_threshold = NULL) + theme_bw()
dev.off()

#### plot switch task ####
rt_switch <-  LDT_clean %>%
  filter(which_task == "Switch")  %>%
  filter(condition == "GG" | condition == "GE" | condition == "EE" | condition == "EG" ) %>% 
  select(condition, RT)

# drop unused levels
rt_switch$condition <- droplevels(rt_switch$condition)
rt_switch$condition <- factor(rt_switch$condition, levels = c("GG","EG","EE","GE"))

# get factors for switch and laguage
for(i in 1:nrow(rt_switch)){
  if(rt_switch$condition[i] == "GG" | rt_switch$condition[i] == "EE"){
    rt_switch$Switch[i] = "Nonswitch"
  }
  else{
    rt_switch$Switch[i] = "Switch"
  }
}
for(i in 1:nrow(rt_switch)){
  if(rt_switch$condition[i] == "GG" | rt_switch$condition[i] == "EG"){
    rt_switch$Language[i] = "German"
  }
  else{
    rt_switch$Language[i] = "English"
  }  
}

rt_switch$Switch <- factor(rt_switch$Switch, levels = c("Nonswitch", "Switch"))
rt_switch$Language <- factor(rt_switch$Language, levels = c("German", "English"))

# get sd and se for RT distribution
library(plotrix)
rt_switch_summary <- rt_switch %>% 
  group_by(Language,Switch) %>% 
  summarise(mRT = mean(RT), sd = sd(RT), se = std.error(RT))

v_switch_summary <- task_3_plot_important %>% 
  group_by(Language,Switch) %>% 
  summarise(Drift = mean(v), sd = sd(v), se = std.error(v))
  

## Plot switch
pd <- position_dodge(0)
bxp3_rt <- ggplot(rt_switch_summary,aes(x = Switch, y = mRT, linetype = Language, group = Language)) +
  geom_point(position = pd, size=1.5) +
  geom_line() +
  geom_errorbar(aes(ymin = mRT - se, ymax = mRT + se), width = .1, position = pd) + 
  theme_bw() +
  ggtitle("Mean RT") 

bxp3_v <- ggplot(v_switch_summary,aes(x = Switch, y = Drift, linetype = Language, group = Language)) +
  geom_point(position = pd, size=1.5) +
  geom_line() +
  geom_errorbar(aes(ymin = Drift - se, ymax = Drift + se), width = .1, position = pd) + 
  theme_bw() +
  ggtitle("Drift rate") 

bxp3_areaLat <- ggplot(task_3_plot_important, aes(x = Switch, y = areaLat)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Area latency") +
  facet_wrap( ~ Language)

bxp3_meanAmp <- ggplot(task_3_plot_important, aes(x = Switch, y = meanAmp)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Mean amplitude") +
  facet_wrap( ~ Language)

filename <- file.path(figureFolder,"Task_3_switch_language.png")
png(filename,pointsize = 20,width=1000, height=600,units = "px")

ggarrange(bxp3_rt, bxp3_v, bxp3_areaLat, bxp3_meanAmp, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

dev.off()



## clean workspace
remove(bxp_rt, bxp_a, bxp_areaLat, bxp_meanAmp, bxp_v, bxp3_rt, bxp3_areaLat, bxp3_meanAmp, bxp3_v,
       task_1_2_plot, task_3_plot, task_3_plot_important, i, pd, rt_switch, rt_switch_summary, rt_wf, rt_wf_summary, 
       v_switch_summary)
