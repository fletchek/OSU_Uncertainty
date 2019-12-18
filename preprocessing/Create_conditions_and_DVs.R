rm(list = ls())
library(here)
library(tidyverse)
data_folder = here("experiments/exp_gauges/data/Good/")
block03 <- file.path(data_folder, list.files(data_folder, pattern = 'BLOCK03'))
block04 <- file.path(data_folder, list.files(data_folder, pattern = 'BLOCK04'))
files <- c(block03, block04)

raw_d <-NULL # map_df(files, read_tsv)

for (this_file in files) {
  # include the line below to help with debugging if there are hangups loading the files
  #print (paste("Loading ", this_file,"."))
  temp <- read_tsv(this_file, skip = 1) 
  raw_d = rbind(raw_d,temp)
}
#Identify which participants were in Static, Current, and Cumulative conditions
static_pid = raw_d %>% filter(task == "Static", aid_format == "By_event", trial == 5) %>% pull(ID)
current_pid = raw_d %>% filter(task == "Dynamic", aid_format == "By_event", trial == 5) %>% pull(ID)
cumulative_pid = raw_d %>% filter(task == "Dynamic", aid_format == "Cumulative", trial == 5) %>% pull(ID)

#Create condition and aid factors
raw_d = raw_d %>% mutate(condition = case_when((ID %in% static_pid)  ~ "Static", 
                                               (ID %in% current_pid)  ~ "Current", 
                                               (ID %in% cumulative_pid)  ~ "Cumulative"))

raw_d = raw_d %>% mutate(aided = case_when((aid_format == "None") ~ "Unaided", 
                                           (aid_format == "Cumulative" | aid_format == "By_event") ~ "Aided"))

# create a column to indicate what the choice would have been based on 
# the maximum probability given the current status of the display at the 
# time of response
curr_aid_pick = raw_d %>%  
  select(p_friend:p_hostile) %>% 
  mutate(x = max.col(.),  
         curr_aid_pick = case_when(x == 1 ~ "Friendly", 
                                   x == 2 ~ "Civilian", 
                                   x == 3 ~ "Hostile"))  %>%  select(curr_aid_pick)

# create a column to indicate what the choice would have been based on 
# the maximum probability given the cumulative status of the display at the 
# time of response
cum_aid_pick = raw_d %>% 
  select(cum_p_friend:cum_p_hostile) %>% 
  mutate(x = max.col(.),  
         cum_aid_pick = case_when(x == 1 ~ "Friendly", 
                                  x == 2 ~ "Civilian", 
                                  x == 3 ~ "Hostile")) %>% select(cum_aid_pick)

# create a column to indicate what the curresponding current guage probability was for the chosen ID
resp_curr_p = raw_d %>%  
  mutate(resp_curr_p = case_when(response == "Friendly" ~ p_friend, 
                                 response == "Civilian" ~ p_unknown, 
                                 response == "Hostile" ~ p_hostile)) %>% select(resp_curr_p)

# create a column to indicate what the curresponding cumulative guage probability was for the chosen ID
resp_cum_p = raw_d %>%  
  mutate(resp_cum_p = case_when(response == "Friendly" ~ cum_p_friend, 
                                response == "Civilian" ~ cum_p_unknown, 
                                response == "Hostile" ~ cum_p_hostile))  %>%  select(resp_cum_p)

raw_d = cbind(raw_d,curr_aid_pick,resp_curr_p,cum_aid_pick,resp_cum_p)
raw_d = raw_d %>% mutate(curr_aid_match = as.integer(response == curr_aid_pick))
raw_d = raw_d %>% mutate(cum_aid_match = as.integer(response == cum_aid_pick))

# Making Wordgauge variable
raw_d$Fcumgauge <- ifelse(raw_d$cum_p_friend > raw_d$cum_p_hostile & raw_d$cum_p_friend > raw_d$cum_p_unknown, 1, 0)
raw_d$Hcumgauge <- ifelse(raw_d$cum_p_hostile > raw_d$cum_p_friend & raw_d$cum_p_hostile > raw_d$cum_p_unknown, 1, 0)
raw_d$Ccumgauge <- ifelse(raw_d$cum_p_unknown > raw_d$cum_p_friend & raw_d$cum_p_unknown > raw_d$cum_p_hostile, 1, 0)
raw_d$Allcumgauge <- ifelse(raw_d$Fcumgauge == 1, 1, ifelse(raw_d$Ccumgauge == 1, 2, ifelse(raw_d$Hcumgauge == 1, 3, 0)))
raw_d$wordgauge <- ifelse(raw_d$Allcumgauge == 1, "Friend", ifelse(raw_d$Allcumgauge == 2, "Civilian", ifelse(raw_d$Allcumgauge == 3, "Hostile", 0)))

# create columns for trial, hits, and false alarm counts
raw_d = raw_d %>% mutate(Friend_trial = as.integer(signal == "Friendly"))
raw_d = raw_d %>% mutate(Not_Friend_trial = as.integer(signal != "Friendly"))
raw_d = raw_d %>% mutate(Friend_hit = as.integer((response == "Friendly") & (signal == "Friendly")))
raw_d = raw_d %>% mutate(Friend_FA = as.integer((response == "Friendly") & (signal != "Friendly")))
raw_d = raw_d %>% mutate(Hostile_trial = as.integer(signal == "Hostile"))
raw_d = raw_d %>% mutate(Not_Hostile_trial = as.integer(signal != "Hostile"))
raw_d = raw_d %>% mutate(Hostile_hit = as.integer((response == "Hostile") & (signal == "Hostile")))
raw_d = raw_d %>% mutate(Hostile_FA = as.integer((response == "Hostile") & (signal != "Hostile")))
raw_d = raw_d %>% mutate(Civilian_trial = as.integer(signal == "Civilian"))
raw_d = raw_d %>% mutate(Not_Civilian_trial = as.integer(signal != "Civilian"))
raw_d = raw_d %>% mutate(Civilian_hit = as.integer((response == "Civilian") & (signal == "Civilian")))
raw_d = raw_d %>% mutate(Civilian_FA = as.integer((response == "Civilian") & (signal != "Civilian")))


# # Ensure that relevent variables are factors to support ANOVA analysis
# raw_d$wordgauge = factor(raw_d$wordgauge, levels = c("Friend","Hostile", "Civilian"), ordered = TRUE)
# raw_d$condition = factor(raw_d$condition, levels = c("Static","Current","Cumulative"), ordered = TRUE)
# raw_d$aided = factor(raw_d$aided, levels = c("Unaided","Aided"), ordered = TRUE)
# raw_d$signal = factor(raw_d$signal, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE) 
# raw_d$response = factor(raw_d$response, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE) 
# raw_d$cum_aid_pick = factor(raw_d$cum_aid_pick, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE) 
# raw_d$curr_aid_pick = factor(raw_d$curr_aid_pick, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE) 
# raw_d$ID = factor(raw_d$ID)


save(raw_d, file = here("/Preprocessing/Raw_d.Rdata"))
