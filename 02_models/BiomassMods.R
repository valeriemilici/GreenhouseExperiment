##### Biomass Models ########

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

perf_dat <- dat %>% filter(status == 1) #living plants only

#error on line 114- all Guapst Controls were harvested on day 212
perf_dat[114,24] <- 212

## All treatments
perf_dat$soil <- factor(perf_dat$soil,
                        levels = c( "conspecific", "heterospecific","control"))
### Above vs. Belowground Allocation -------------------------------------------

bioratio.mod <- lmer(log(Root_Mass/Stem_Mass) ~ Treatment*soil + duration +
                                 log(Init_Height) + log(Tot_Biomass) +
                                 (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin), 
                               data = perf_dat)

check_model(bioratio.mod) #looks ok
summary(bioratio.mod)
Anova(bioratio.mod)
#no real patterns to report. To be expected based on the data. Only differences
#are between the sterile and live treatments. Microbes affect allocation, but 
#not between microbe treatments. 

saveRDS(bioratio.mod, file = "02_models/modoutput/bmratio.RDS")

#quick prediction plotting
rat_pred<- as.data.frame(ggemmeans(bioratio.mod, c("Treatment", "soil")))

plot(ggpredict(bioratio.mod, c("Treatment", "soil")))

ggplot(rat_pred) + 
  geom_pointrange(aes(x = x, y = predicted, ymin = conf.low,
                      ymax = conf.high, color = group),
                  position = position_dodge(width = 0.5))+
  theme_bw(16)

### Total Biomass --------------------------------------------------------------
bm.mod0 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment*soil + duration +
                                 scale(log(Init_Height)) +
                                 (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin), 
                               data = perf_dat)

bm.mod <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment*soil + duration +
                           scale(Init_Height) +
                           (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin), 
                         data = perf_dat)

check_model(bm.mod) 
anova(bm.mod0, bm.mod) 
summary(bm.mod)
#In the wet treatment, seedlings in heterospecific soil grew larger than seedlings
#in conspecific soil. 
saveRDS(bm.mod, file = "02_models/modoutput/biomass.RDS")