library (tidyverse)
library (coda)
library(dplyr)

rm(list = ls())

load(here("analysis/RT_cond_Aid_x_Condition.Rdata"))


#### Means and 95% highest-density intervals of MCMC sample chains.
# Max_Rhat =max(samples$BUGSoutput$summary[,8])`   
# DIC =samples$BUGSoutput$DIC`

summarise <- function(my.samples, theta, diff = FALSE) {
  # my.samples should be BUGSoutput produced by jags
  # theta is a string that names the parameter that we want summarised
  theta.samples <- as.mcmc(my.samples$sims.matrix[,theta]) # get samples for the parameter we are interested in
  M = mean(theta.samples)
  lb = HPDinterval(theta.samples)[1] # lower bound of HDI
  ub = HPDinterval(theta.samples)[2] # upper bound
  
  # if these data are difference scores
  if (diff == TRUE)
  {
    p.neg = length(theta.samples[theta.samples < 0])/length(theta.samples) * 100 # percentage of posterior distribution < 0
    p.pos = length(theta.samples[theta.samples > 0])/length(theta.samples) * 100 # percentage > 0
    sprintf ('M = %.2f, HDI = [%.2f, %.2f], %.0f%% < 0 < %1.0f%%', M, lb, ub, p.neg, p.pos)
  } else
  {
    p.neg = length(theta.samples[theta.samples < 0])/length(theta.samples) * 100 # percentage of posterior distribution < 0
    p.pos = length(theta.samples[theta.samples > 0])/length(theta.samples) * 100 # percentage > 0
    sprintf ('M = %.2f, HDI = [%.2f, %.2f]', M, lb, ub)
  }
}

p_data <- function(my.samples, theta) {
  # my.samples should be BUGSoutput produced by jags
  # theta is a string that names the parameter that we want summarised
  theta.samples <- as.mcmc(my.samples$sims.matrix[,theta]) # get samples for the parameter we are interested in
  M = mean(theta.samples)
  lb = HPDinterval(theta.samples)[1] # lower bound of HDI
  ub = HPDinterval(theta.samples)[2] # upper bound
  
  params = data.frame(M,lb,ub)
  return(params)
  
}


Model_RT_cond_cell = data.frame(condition = c(rep("Static, By Event",2),rep("Dynamic, By Event",2),rep("Dynamic, Cumulative",2)),
                      aided = c("Unaided","Aided"),
                      RT_cond_cell_M = rep(NA,6),
                      RT_cond_cell_lb = rep(NA,6),
                      RT_cond_cell_ub = rep(NA,6),
                      RT_cond_cell_M = rep(NA,6),
                      RT_cond_cell_lb = rep(NA,6),
                      RT_cond_cell_ub = rep(NA,6))
Model_RT_cond_cell$condition = factor(Model_RT_cond_cell$condition, levels = c("Static, By Event","Dynamic, By Event", "Dynamic, Cumulative"), ordered = TRUE)
Model_RT_cond_cell$aided = factor(Model_RT_cond_cell$aided, levels = c("Unaided","Aided"), ordered = TRUE)

Model_RT_cond_cell[1,3:5] = p_data(samples$BUGSoutput, 'd_m[1,1]')
Model_RT_cond_cell[2,3:5] = p_data(samples$BUGSoutput, 'd_m[1,2]')
Model_RT_cond_cell[3,3:5] = p_data(samples$BUGSoutput, 'd_m[2,1]')
Model_RT_cond_cell[4,3:5] = p_data(samples$BUGSoutput, 'd_m[2,2]')
Model_RT_cond_cell[5,3:5] = p_data(samples$BUGSoutput, 'd_m[3,1]')
Model_RT_cond_cell[6,3:5] = p_data(samples$BUGSoutput, 'd_m[3,2]')
Model_RT_cond_cell[1,6:8] = p_data(samples$BUGSoutput, 'c_m[1,1]')
Model_RT_cond_cell[2,6:8] = p_data(samples$BUGSoutput, 'c_m[1,2]')
Model_RT_cond_cell[3,6:8] = p_data(samples$BUGSoutput, 'c_m[2,1]')
Model_RT_cond_cell[4,6:8] = p_data(samples$BUGSoutput, 'c_m[2,2]')
Model_RT_cond_cell[5,6:8] = p_data(samples$BUGSoutput, 'c_m[3,1]')
Model_RT_cond_cell[6,6:8] = p_data(samples$BUGSoutput, 'c_m[3,2]')

save(Model_RT_cond_cell,file = here("analysis/Model_RT_cond_cell_CI.Rdata"))


# comparison = c("Aid_Static","Static.Current_Unaided","Static.Current_Aided","Aid_Current","Current.Cumulative_Unaided","CurrentCumulative_Aided")
Model_RT_cond_inf =data.frame(condition = c("Aid_Static","Static.Current_Unaided","Static.Current_Aided",
                                    "Aid_Current","Current.Cumulative_Unaided","CurrentCumulative_Aided"),
                      RT_cond_inf_M = rep(NA,6),
                      RT_cond_inf_lb = rep(NA,6),
                      RT_cond_inf_ub = rep(NA,6),
                      RT_cond_inf_M = rep(NA,6),
                      RT_cond_inf_lb = rep(NA,6),
                      RT_cond_inf_ub = rep(NA,6))
# Model_RT_cond_inf$condition = factor(Model_RT_cond_inf$condition, levels = c("Static, By Event","Dynamic, By Event", "Dynamic, Cumulative"), ordered = TRUE)
# Model_RT_cond_inf$aided = factor(Model_RT_cond_inf$aided, levels = c("Unaided","Aided"), ordered = TRUE)

Model_RT_cond_inf[1,2:4] = p_data(samples$BUGSoutput, 'd_aid.static')
Model_RT_cond_inf[2,2:4] = p_data(samples$BUGSoutput, 'd_display.static.current')
Model_RT_cond_inf[3,2:4] = p_data(samples$BUGSoutput, 'd_static.current.aided')
Model_RT_cond_inf[4,2:4] = p_data(samples$BUGSoutput, 'd_aid.current')
Model_RT_cond_inf[5,2:4] = p_data(samples$BUGSoutput, 'd_display.current.cumulative')
Model_RT_cond_inf[6,2:4] = p_data(samples$BUGSoutput, 'd_current.cumulative.aided')
Model_RT_cond_inf[1,5:7] = p_data(samples$BUGSoutput, 'c_aid.static')
Model_RT_cond_inf[2,5:7] = p_data(samples$BUGSoutput, 'c_display.static.current')
Model_RT_cond_inf[3,5:7] = p_data(samples$BUGSoutput, 'c_static.current.aided')
Model_RT_cond_inf[4,5:7] = p_data(samples$BUGSoutput, 'c_aid.current')
Model_RT_cond_inf[5,5:7] = p_data(samples$BUGSoutput, 'c_display.current.cumulative')
Model_RT_cond_inf[6,5:7] = p_data(samples$BUGSoutput, 'c_current.cumulative.aided')

save(Model_RT_cond_inf,file = here("analysis/Model_RT_cond_inf_CI.Rdata"))

