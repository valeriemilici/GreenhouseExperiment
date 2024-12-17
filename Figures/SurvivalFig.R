## Survival Fig###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())
library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(sjPlot) #forest plots
library(ggthemes)
library(lme4) #extracting fixed effects
library(lmerTest) #p-values

dat <- read.csv("Data/mod_data_allpots.csv")
surv.mod <- readRDS("02_models/modoutput/survival.RDS")
load("03_Bootstrapping/BootOutput/surv.boot")

# Create a function for getting the mean model prediction
surv.fun <- function(.){
  preddat <- expand.grid(Treatment = c("D", "W"),
                         soil = c("conspecific", "heterospecific", "control"),
                         Init_Height = 7,
                         Cohort = 1)
  predict(., newdata = preddat, re.form = ~0)
}

# Create a data frame of preds + CIs to plot
surv.dat <- data.frame(preds =surv.fun(surv.mod),
                    confint(surv.boot),
                    trt = c("50% WHC", "100% WHC"),
                    soil = c("Conspecific", "Conspecific",
                             "Heterospecific", "Heterospecific",
                             "Control", "Control"))
names(surv.dat)[2:3] <- c("lwr", "upr")

# Relevel treatments 
surv.dat$trt <- factor(surv.dat$trt, levels = c("50% WHC", "100% WHC"))
surv.dat$soil <- factor(surv.dat$soil, levels = c("Control", "Conspecific", 
                                                  "Heterospecific"))

# Plot the figure
fig4 <- ggplot(surv.dat, aes(soil, plogis(preds))) +
  geom_pointrange(mapping = aes(ymin = plogis(lwr), ymax = plogis(upr))) + 
  facet_wrap(~trt) +
  labs(x = "Soil Inoculum Source",
       y = "Probability of Survival") +
  theme_bw(16) 

fig4
# This looks kind of dumb. The mean model estimate may vary, but survival
# was so high in all treatments that the CIs all touch 1.00. 

ggsave(plot = fig4, filename = "04_figures/figures/Survival.png",
       height = 3, width = 6, units = "in")
