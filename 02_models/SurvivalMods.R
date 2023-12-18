### Survival Models ###

### Initialize Workspace -------------------------------------------------------

rm(list = ls())

library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(lme4) #runs models
library(lmerTest) #model testing, gives p values
library(jtools) #checking out cloglog results
library(performance) #diagnostics
library(DHARMa) #additional diagnostics
library(sjPlot)
library(effects)

dat <- read.csv("Data/mod_data_allpots.csv")

### Data Manipulation ----------------------------------------------------------

#Add live vs. sterile soil column
dat <- dat %>% mutate(microbes = if_else(soil == "control", "sterile", "live"))

### Basic Survival Model -------------------------------------------------------
dat$soil <- factor(dat$soil, levels =c("control", "conspecific", "heterospecific"))

surv <- glmer(status ~ soil*Treatment + scale(Init_Height) + as.factor(Cohort) +
                 (1|Seedling)+ (1|Inoc_Sp) + (1|Table) + (1|begin),
               family = binomial(link = "logit"),
               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)),
               data = dat)
testDispersion(surv)
survOutput <- simulateResiduals(surv, plot=F)
plot(survOutput) #looks way better with these changes. Cohort is a much better
#fit than duration here. 

summary(surv)

plotResiduals(survOutput, dat$soil)

saveRDS(surv, "02_models/modoutput/survival.RDS")
