rm(list = ls())
load(here("/preprocessing/Filt_d.Rdata"))

pid_list = unique(filt_d$ID)

static_list = filt_d %>% filter(condition == "Static", trial == 3) %>% select(ID,condition,aided)
current_list = filt_d %>% filter(condition == "Current", trial == 3) %>% select(ID,condition,aided)
cum_list = filt_d %>% filter(condition == "Cumulative", trial == 3) %>% select(ID,condition,aided)


#Create tibbles to hold accuracy data
acc_static = data.frame(matrix(data = NA, nrow = dim(static_list)[1],ncol = 14))
colnames(acc_static) = c("resp_acc","norm_acc",
                         "resp_FH","resp_HH","resp_CH","resp_FFA","resp_HFA","resp_CFA",
                         "norm_FH","norm_HH","norm_CH","norm_FFA","norm_HFA","norm_CFA")
cmat_acc_static = cbind(static_list,acc_static)

acc_current = data.frame(matrix(data = NA, nrow = dim(current_list)[1],ncol = 14))
colnames(acc_current) = c("resp_acc","norm_acc",
                          "resp_FH","resp_HH","resp_CH","resp_FFA","resp_HFA","resp_CFA",
                          "norm_FH","norm_HH","norm_CH","norm_FFA","norm_HFA","norm_CFA")
cmat_acc_current = cbind(current_list,acc_current)

acc_cum = data.frame(matrix(data = NA, nrow = dim(cum_list)[1],ncol = 14))
colnames(acc_cum) = c("resp_acc","norm_acc",
                      "resp_FH","resp_HH","resp_CH","resp_FFA","resp_HFA","resp_CFA",
                      "norm_FH","norm_HH","norm_CH","norm_FFA","norm_HFA","norm_CFA")
cmat_acc_cum = cbind(cum_list,acc_cum)


for(t in 1:dim(static_list)[1]){
  temp_list = filt_d %>% filter(ID == static_list[t,1], condition == static_list[t,2], aided == static_list[t,3]) %>% select(ID,condition,aided,signal,response,cum_aid_pick)
  cmat_static_resp = confusionMatrix(temp_list$response,temp_list$signal)
  cmat_static_cum = confusionMatrix(temp_list$cum_aid_pick,temp_list$signal)
  cmat_acc_static[t,4] = cmat_static_resp$overall[1]
  cmat_acc_static[t,5] = cmat_static_cum$overall[1]
  cmat_acc_static[t,6:8] = cmat_static_resp$byClass[1:3,1]
  cmat_acc_static[t,9:11] =  rep(1,3) - cmat_static_resp$byClass[1:3,2]
  cmat_acc_static[t,12:14] = cmat_static_cum$byClass[1:3,1]
  cmat_acc_static[t,15:17] =  rep(1,3) - cmat_static_cum$byClass[1:3,2]
}

for(t in 1:dim(current_list)[1]){
  temp_list = filt_d %>% filter(ID == current_list[t,1], condition == current_list[t,2], aided == current_list[t,3]) %>% select(ID,condition,aided,signal,response,cum_aid_pick)
  cmat_current_resp = confusionMatrix(temp_list$response,temp_list$signal)
  cmat_current_cum = confusionMatrix(temp_list$cum_aid_pick,temp_list$signal)
  cmat_acc_current[t,4] = cmat_current_resp$overall[1]
  cmat_acc_current[t,5] = cmat_current_cum$overall[1]
  cmat_acc_current[t,6:8] = cmat_current_resp$byClass[1:3,1]
  cmat_acc_current[t,9:11] =  rep(1,3) - cmat_current_resp$byClass[1:3,2]
  cmat_acc_current[t,12:14] = cmat_current_cum$byClass[1:3,1]
  cmat_acc_current[t,15:17] =  rep(1,3) - cmat_current_cum$byClass[1:3,2]
}

for(t in 1:dim(cum_list)[1]){
  temp_list = filt_d %>% filter(ID == cum_list[t,1], condition == cum_list[t,2], aided == cum_list[t,3]) %>% select(ID,condition,aided,signal,response,cum_aid_pick)
  cmat_cum_resp = confusionMatrix(temp_list$response,temp_list$signal)
  cmat_cum_cum = confusionMatrix(temp_list$cum_aid_pick,temp_list$signal)
  cmat_acc_cum[t,4] = cmat_cum_resp$overall[1]
  cmat_acc_cum[t,5] = cmat_cum_cum$overall[1]
  cmat_acc_cum[t,6:8] = cmat_cum_resp$byClass[1:3,1]
  cmat_acc_cum[t,9:11] =  rep(1,3) - cmat_cum_resp$byClass[1:3,2]
  cmat_acc_cum[t,12:14] = cmat_cum_cum$byClass[1:3,1]
  cmat_acc_cum[t,15:17] =  rep(1,3) - cmat_cum_cum$byClass[1:3,2]
}

cmat_accuracy = rbind(cmat_acc_static,cmat_acc_current,cmat_acc_cum)
cmat_accuracy = cmat_accuracy %>%mutate(norm_dif = (resp_acc - norm_acc),
                                        norm_dif_FH = (resp_FH - norm_FH), norm_dif_HH = (resp_HH - norm_HH),norm_dif_CH = (resp_CH - norm_CH),
                                        norm_dif_FFA = (resp_FFA - norm_FFA), norm_dif_HFA = (resp_HFA - norm_HFA),norm_dif_CFA = (resp_CFA - norm_CFA))

save(cmat_accuracy, file = here("/analysis/Accuracy.Rdata"))

# rearrange accuracy file to create signal colum
acc_friend = cmat_accuracy %>% select(ID,condition,aided,resp_FH,resp_FFA,norm_dif_FH,norm_dif_FFA) %>%
  rename(resp_HR = resp_FH, resp_FAR = resp_FFA, norm_dif_HR = norm_dif_FH, norm_dif_FAR = norm_dif_FFA)
acc_friend = cbind(acc_friend, "signal" = rep("Friendly",dim(acc_friend)[1]))
acc_hostile = cmat_accuracy %>% select(ID,condition,aided,resp_HH,resp_HFA,norm_dif_HH,norm_dif_HFA) %>%
  rename(resp_HR = resp_HH, resp_FAR = resp_HFA, norm_dif_HR = norm_dif_HH, norm_dif_FAR = norm_dif_HFA)
acc_hostile = cbind(acc_hostile, "signal" = rep("Hostile",dim(acc_friend)[1]))
acc_signal = rbind(acc_friend,acc_hostile)

save(acc_signal, file = here("/analysis/Acc_signal.Rdata"))