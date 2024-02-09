### Accounting for seedling deaths in glmmTMB

### Initializing Workspace -----------------------------------------------------

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
library(car)
library(ggeffects)
library(glmmTMB)

dat <- read.csv("Data/mod_data_allpots.csv")



#error on line 114- all Guapst Controls were harvested on day 212
#perf_dat[114,24] <- 212

#if a seedling is dead, give it biomass 0 instead of NA
dat$Tot_Biomass <- if_else(dat$status == 0, 0, dat$Tot_Biomass)

dat2 <- dat %>%
    mutate(Init_Height = replace_na(Init_Height, 0),
           check = log1p((Init_Height)),
           trn_init_height = scale(check),
           bcheck = log1p(Tot_Biomass),
           trn_tot_bm = scale(bcheck))
          



## All treatments
dat2$soil <- factor(dat2$soil,
                        levels = c(  "heterospecific", "conspecific","control"))


# mod1 <- glmmTMB(log1p(Tot_Biomass) ~ scale(log1p(Init_Height)) +
#                   Treatment*soil + as.factor(Cohort) +
#   (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
# ziformula = ~ scale(log1p(Init_Height)) + Treatment*soil + as.factor(Cohort),
# data = dat2)


# 
# mod1_a <- glmmTMB(Tot_Biomass ~ scale(log1p(Init_Height)) + ## ok, problem is the log - switch to Gamma?
#                     Treatment*soil + as.factor(Cohort) +
#                     (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
#                   ziformula = ~ scale(log1p(Init_Height)) + 
#                     Treatment*soil + as.factor(Cohort),
#                   data = dat2, family = lognormal(link = "log")) ## doesn't work
summary(dat2$Tot_Biomass)
mod1_b <- glmmTMB(Tot_Biomass ~ scale(log1p(Init_Height)) + ## ok, problem is the log - switch to Gamma?
                    Treatment*soil + as.factor(Cohort) +
                    (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
                  ziformula = ~ scale(log1p(Init_Height)) + 
                    Treatment*soil + as.factor(Cohort),
                  data = dat2, family = ziGamma(link = "log"))

summary(mod1_b)


