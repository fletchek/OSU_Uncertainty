rm(list = ls())
load(here("/preprocessing/Raw_d.Rdata"))

# Remove participants with fewer Friend or Hostile signals in a block
sig_trial_limit = 5
low_sig_trials <- raw_d %>% 
  group_by(ID, condition, aided) %>% 
  summarize(F_trials = sum(Friend_trial),
            H_trials = sum(Hostile_trial)) %>% 
  filter(F_trials < sig_trial_limit | H_trials < sig_trial_limit) %>% 
  pull(ID) %>% 
  unique()
low_sig_num = length(low_sig_trials)
filt_d = raw_d %>% filter(!(ID %in% low_sig_trials))

# filter out participants who responded Civilian every trial in one block or the other
# this screening was NOT preregistered, because we didn't anticipate it, but affects small # of subjects
mindless <- filt_d %>% 
  group_by(ID, condition, aided) %>% 
  summarize(F_resp = sum(response == "Friendly"),
            C_resp = sum(response == "Civilian"),
            H_resp = sum(response == "Hostile")) %>% 
  filter(F_resp == 0 & H_resp == 0) %>% 
  pull(ID) %>% 
  unique()
mindless_num = length(mindless)
filt_d = filt_d %>% filter(!(ID %in% mindless))

#Remove slow trials
threshold_rt = 60
long_rt = filt_d %>% filter(RT > threshold_rt) %>% pull(RT)
long_rt_num = length(long_rt)
filt_d = filt_d %>% filter(RT <= threshold_rt)

# create factors

filt_d$condition = factor(filt_d$condition, levels = c("Static","Current", "Cumulative"), ordered = TRUE)
filt_d$aided = factor(filt_d$aided, levels = c("Unaided","Aided"), ordered = TRUE)
filt_d$signal = factor(filt_d$signal, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE)
filt_d$response = factor(filt_d$response, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE)
filt_d$cum_aid_pick = factor(filt_d$cum_aid_pick, levels = c("Friendly","Hostile","Civilian"), ordered = TRUE)
filt_d$ID = as.factor(filt_d$ID)


save(filt_d, file = here("/preprocessing/Filt_d.Rdata"))
rm(raw_d,filt_d)
save(list = ls(), file = here("/preprocessing/Filt_variables.Rdata"))
