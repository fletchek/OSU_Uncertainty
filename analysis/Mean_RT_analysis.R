rm(list = ls())
load(here("/preprocessing/Filt_d.Rdata"))

RT_cond <- filt_d %>%
  group_by(ID, condition, aided) %>%
  summarize(RT_mean = mean(RT)) %>%
  ungroup()

RT_signal <- filt_d %>%
  group_by(ID, condition, aided, signal) %>%
  filter(signal != "Civilian") %>%
  summarize(RT_mean = mean(RT)) %>%
  ungroup()

save(RT_cond, file = here("/analysis/RT_cond.Rdata"))
save(RT_signal, file = here("/analysis/RT_signal.Rdata"))