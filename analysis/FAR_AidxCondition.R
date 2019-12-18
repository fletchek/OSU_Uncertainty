library(tidyverse)
library(readr)
library(here)

# Script to generate estimates of d' and model predictions
# for Bartlett & McCarley study of effects of aid format
# on use of signal detection aid. 
# Fits a 2 (BLOCK: Aided vs Unaided) x 3 (FORMAT: Raw, Llikelihood ratio, Confidence) model.
# UPDATED: 2018 AUG 19

## First, load & preprocess data

### clear the workspace
rm(list = ls())

expt = "FAR"

# p = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(p)

load(here("/analysis/Accuracy.Rdata"))
cmat_accuracy$ID.num <- match(cmat_accuracy$ID, unique(cmat_accuracy$ID))

### create a column with numeric codes for the display format
# dat$format <- as.factor(dat$format) # order will be confidence, likelihood, raw
# dat$format <- factor(dat$format, levels(dat$format)[c(3, 2, 1)]) # put in order raw, likelihood, confidence
cmat_accuracy$condition.code <- match(cmat_accuracy$condition, levels(cmat_accuracy$condition))

### same thing for block
# dat$block <- as.factor(dat$block)
# dat$block <- factor(dat$block, levels(dat$block)[c(2, 1)]) # put in order unaided, aide
cmat_accuracy$aided.code <- match(cmat_accuracy$aided, levels(cmat_accuracy$aided))


### calculate d' and c, then filter out subjects with d' below cutoff of 0.5
# cmat_accuracy <- cmat_accuracy %>% mutate(resp_HR = (resp_FH + resp_HH)/2+.001)
# cmat_accuracy <- cmat_accuracy %>% mutate(resp_FAR = (resp_FFA + resp_HFA)/2+.001)
# cmat_accuracy <- cmat_accuracy %>%
#   mutate(d = qnorm(resp_HR) - qnorm(resp_FAR), c = -.5 * (qnorm(resp_HR) + qnorm(resp_FAR))) #%>%
#   # filter(d >= 0.5)

# # we will need efficiency scores to calculate correlations laters
# # get those and save them now
# eff <- cmat_accuracy %>% filter(aided == 'Aided') %>% mutate(score = (d/sqrt(d^2 + 3^2))^2) %>% select(ID, score) %>% arrange(ID)
# write_csv(eff, "Efficiency_scores.csv")

## Now set up modeling
library(runjags)
library(rjags)
library(R2jags)

nBurnin <- 10000
thinSteps <- 1
nChains <- 4
samplesToKeep <- 100000
nIter <- ceiling((samplesToKeep * thinSteps) / nChains) + nBurnin

modelString <-
  "
model{
for (i in 1:k){
d[i] ~ dnorm(d_mu + d_form[form[i]] + d_block[block[i]] + d_formxblock[form[i], block[i]] + d_s[ID[i]], 1/d_ySigma^2)
c[i] ~ dnorm(c_mu + c_form[form[i]] + c_block[block[i]] + c_formxblock[form[i], block[i]] + c_s[ID[i]], 1/c_ySigma^2)
}

###############################################
# Priors on hyperparameters
d_mu ~ dnorm(d_M, 1/(100 * d_SD)^2)
c_mu ~ dnorm(c_M, 1/(100 * c_SD)^2)


for (i in 1:3){
d_form[i] ~ dnorm(0,1/d_form.sigma^2)
c_form[i] ~ dnorm(0,1/c_form.sigma^2)
}

for (i in 1:2){
d_block[i] ~ dnorm(0, 1/d_block.sigma^2)
c_block[i] ~ dnorm(0, 1/c_block.sigma^2)
}

for (i in 1:3){
for (j in 1:2){
d_formxblock[i, j]  ~ dnorm(0, 1/d_formxblock.sigma^2) 
c_formxblock[i, j]  ~ dnorm(0, 1/c_formxblock.sigma^2) 
}}

for (i in 1:n){ 
d_s[i] ~ dnorm(0, 1/d_s.sigma^2)
c_s[i] ~ dnorm(0, 1/c_s.sigma^2)
}

d_ySigma ~ dunif(d_SD/1000, d_SD*1000)
c_ySigma ~ dunif(c_SD/1000, c_SD*1000)


# set priors on deflections following Kruschke (2015)
# prior is gamma distribution with rate 2 * SD
# and mode SD/2
d_gMode <- d_SD/2; d_gSD <- 2 * d_SD
d_gRate <- (d_gMode + sqrt(d_gMode^2 + 4 * d_gSD^2))/(2 * d_gSD ^ 2)
d_gShape <- 1 + d_gMode * d_gRate
d_form.sigma ~ dgamma(d_gShape, d_gRate)
d_block.sigma ~ dgamma(d_gShape, d_gRate)
d_formxblock.sigma ~ dgamma(d_gShape, d_gRate)
d_s.sigma ~ dgamma(d_gShape, d_gRate)

c_gMode <- c_SD/2; c_gSD <- 2 * c_SD
c_gRate <- (c_gMode + sqrt(c_gMode^2 + 4 * c_gSD^2))/(2 * c_gSD ^ 2)
c_gShape <- 1 + c_gMode * c_gRate
c_form.sigma ~ dgamma(c_gShape, c_gRate)
c_block.sigma ~ dgamma(c_gShape, c_gRate)
c_formxblock.sigma ~ dgamma(c_gShape, c_gRate)
c_s.sigma ~ dgamma(c_gShape, c_gRate)

for (i in 1:3){ 
for (j in 1:2) {
d_m[i, j] <- d_mu + d_form[i] + d_block[j] + d_formxblock[i, j]
c_m[i, j] <- c_mu + c_form[i] + c_block[j] + c_formxblock[i, j]
}}


d_b0 <- mean(d_m[1:3,1:2])
for (i in 1:3){
d_b_form[i] <- mean(d_m[i, 1:2]) - d_b0
}
for (j in 1:2){
d_b_block[j] <- mean(d_m[1:3, j]) - d_b0
}
for (i in 1:3){
for (j in 1:2){
d_b_formxblock[i, j] <- d_m[i, j] - (d_b0 + d_b_form[i] + d_b_block[j])
}
}

c_b0 <- mean(c_m[1:3,1:2])
for (i in 1:3){
c_b_form[i] <- mean(c_m[i, 1:2]) - c_b0
}
for (j in 1:2){
c_b_block[j] <- mean(c_m[1:3, j]) - c_b0
}
for (i in 1:3){
for (j in 1:2){
c_b_formxblock[i, j] <- c_m[i, j] - (c_b0 + c_b_form[i] + c_b_block[j])
}
}

# Calculate the aid benefit for each of the three formats
for(i in 1:3){
  d_aid.effect[i] <- d_b_block[2] + d_b_formxblock[i, 2] -  (d_b_block[1] + d_b_formxblock[i, 1])
  c_aid.effect[i] <- c_b_block[2] + c_b_formxblock[i, 2] -  (c_b_block[1] + c_b_formxblock[i, 1])
}

# Calculate the benefit of aiding in each display condition 
d_aid.static = d_m[1,2] - d_m[1,1]
d_aid.current = d_m[2,2] - d_m[2,1]
d_aid.cumulative = d_m[3,2] - d_m[3,1]
c_aid.static = c_m[1,2] - c_m[1,1]
c_aid.current = c_m[2,2] - c_m[2,1]
c_aid.cumulative = c_m[3,2] - c_m[3,1]

# Calculate the benefit of unaided current and cumulative displays 
d_display.static.current = d_m[2,1] - d_m[1,1]
d_display.current.cumulative = d_m[3,1] - d_m[2,1]
d_display.static.cumulative = d_m[3,1] - d_m[1,1]
c_display.static.current = c_m[2,1] - c_m[1,1]
c_display.current.cumulative = c_m[3,1] - c_m[2,1]
c_display.static.cumulative = c_m[3,1] - c_m[1,1]

# Calculate the aid x display interaction 
d_static.current.aided = d_aid.effect[2] - d_aid.effect[1]
d_current.cumulative.aided = d_aid.effect[3] - d_aid.effect[2]
d_static.cumulative.aided = d_aid.effect[3] - d_aid.effect[1]
c_static.current.aided = c_aid.effect[2] - c_aid.effect[1]
c_current.cumulative.aided = c_aid.effect[3] - c_aid.effect[2]
c_static.cumulative.aided = c_aid.effect[3] - c_aid.effect[1]

}
"
writeLines(modelString, con = "TEMP.txt")
parameters <- c(
  "d_mu", "d_form", "d_block", "d_m", 
  "d_aid.effect","d_aid.static","d_aid.current","d_aid.cumulative",
  "d_display.static.current", "d_display.current.cumulative", "d_display.static.cumulative",
  "d_static.current.aided", "d_current.cumulative.aided", "d_static.cumulative.aided",
  "c_mu", "c_form", "c_block", "c_m", 
  "c_aid.effect","c_aid.static","c_aid.current","c_aid.cumulative",
  "c_display.static.current", "c_display.current.cumulative", "c_display.static.cumulative",
  "c_static.current.aided", "c_current.cumulative.aided", "c_static.cumulative.aided"
)

dataList <- list(
  k = length(cmat_accuracy$ID),
  n = length(unique(cmat_accuracy$ID)),
  ID = cmat_accuracy$ID.num,
  form = cmat_accuracy$condition.code,
  block = cmat_accuracy$aided.code,
  d = cmat_accuracy$resp_FFA, c = cmat_accuracy$resp_HFA,
  d_M = mean(cmat_accuracy$resp_FFA), d_SD = sd(cmat_accuracy$resp_FFA),
  c_M = mean(cmat_accuracy$resp_HFA), c_SD = sd(cmat_accuracy$resp_HFA)#,
)

samples <- jags(dataList,
  parameters.to.save = parameters,
  model.file = "TEMP.txt",
  n.chains = nChains, n.iter = nIter, n.burnin = nBurnin, n.thin = thinSteps
)
samples
save(samples,file = paste0(here("analysis"),"/",expt,"_Aid_x_Condition.Rdata"))






