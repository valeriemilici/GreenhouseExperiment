## biomass Figs###

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
bm.mod <- readRDS("02_models/modoutput/biomass.RDS")
bm.ratio <- readRDS("02_models/modoutput/bmratio.RDS")
load("03_Bootstrapping/BootOutput/bm.boot")
load("03_Bootstrapping/BootOutput/bm.PSF")
load("03_Bootstrapping/BootOutput/bmratio.boot")

### Biomass Model Output -------------------------------------------------------

#Remove intercept variation
### Get the SEs without intercept ---------------
preddat <- expand.grid(Treatment = c("D", "W"),
                       soil = c( "control", "conspecific", "heterospecific"))


#variance covariance matrix of the model without nuisance variables
mod.vcov <- vcov(bm.mod)[-c(1,5:6),][,-c(1,5:6)]
#create the design matrix
Xmat <- model.matrix(~Treatment*soil, data = preddat)
XB <- model.matrix(~Treatment*soil, data = preddat)
Xmat<- Xmat[-c(1,2),][,-1] #remove intercept
XBeta <- XB %*% fixef(bm.mod)[-c(5:6)] #remove duration and init height

#matrix algebra between the two
bm.se <- data.frame(se = sqrt(diag(Xmat %*% mod.vcov %*% t(Xmat))))


dry.est <- data.frame(var = c("Heterospecific", "Conspecific"),
                      est = XBeta[c(5,3),],
                      trt = c("D", "D"),
                      se = c(bm.se[3,],bm.se[1,]))

wet.est <- data.frame(var = c("Heterospecific", "Conspecific"),
                      est = XBeta[c(6,4),],
                      trt = c("W", "W"),
                      se = c(bm.se[4,],bm.se[2,]))

mod.est <- rbind(wet.est, dry.est)

mod.est <- mutate(mod.est, upr = est + se, lwr = est - se)

#create descriptive labels for the facts
trt.labs <- as_labeller(c("D" = "50% WHC",
                          "W" = "100% WHC"))
#plot the output --------------------------------
fig1 <- ggplot(mod.est, aes(var, exp(est))) + 
  geom_pointrange(aes(ymin = exp(lwr), ymax = exp(upr)),
                  position = position_dodge(width = 0.4)) +
  facet_wrap(~trt,
             labeller = labeller(trt = trt.labs)) +
  labs(x = "Soil Inoculum Source",
       y = "Total Dry Biomass (g)",
       col = "") +
  theme_bw(16) 
  

fig1 #looks good

ggsave(plot = fig1, filename = "04_figures/figures/bmresponse.pdf",
       height = 3, width = 6, units = "in")

### NPSF Estimates -------------------------------------------------------------

bm.fun <- function(.){
  con.dat <- expand.grid(Treatment = c("D", "W"),
                         soil = "conspecific",
                         duration = 212,
                         Init_Height = 7)
  het.dat <- expand.grid(Treatment = c("D", "W"),
                         soil = "heterospecific",
                         duration = 212,
                         Init_Height = 7)
  con.est <- predict(., newdata = con.dat, re.form = ~0)
  het.est <- predict(., newdata = het.dat, re.form = ~0)
  
  log((exp(con.est)/exp(het.est)))
}

PSF <- data.frame(preds =bm.fun(bm.mod),
                  confint(bm.PSF),
                  trt = c("50% WHC", "100% WHC"))
names(PSF)[2:3] <- c("lwr", "upr")
PSF$trt <- factor(PSF$trt, levels = c("50% WHC", "100% WHC"))

fig2 <- ggplot(PSF, aes(trt, (preds))) +
  geom_pointrange(mapping = aes(ymin = (lwr), ymax = (upr))) + 
  labs(x = "Soil Moisture Treatment",
       y = "Log Response Ratio") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw(16) 

fig2

ggsave(plot = fig2, filename = "04_figures/figures/bmPSF.pdf",
       height = 3, width = 3, units = "in")

### Biomass Ratio Figure -------------------------------------------------------

bmratio.fun <- function(.){
  preddat <- expand.grid(Treatment = c("D", "W"),
                         soil = c("conspecific", "heterospecific", "control"),
                         duration = 212,
                         Init_Height = 7,
                         Tot_Biomass = 0.30)
  predict(., newdata = preddat, re.form = ~0)
}

ratio <- data.frame(preds =bmratio.fun(bm.ratio),
                  confint(bmratio.boot),
                  trt = c("50% WHC", "100% WHC"),
                  soil = c("Conspecific", "Conspecific",
                           "Heterospecific", "Heterospecific",
                           "control", "control"))
names(ratio)[2:3] <- c("lwr", "upr")

ratio <- ratio %>% filter(soil != "control")
ratio$trt <- factor(ratio$trt, levels = c("50% WHC", "100% WHC"))

fig3 <- ggplot(ratio, aes(soil, exp(preds))) +
  geom_pointrange(mapping = aes(ymin = exp(lwr), ymax = exp(upr))) + 
  facet_wrap(~trt) +
  labs(x = "Soil Inoculum Source",
       y = "Biomass Ratio (R:S)") +
  theme_bw(16) 

fig3

ggsave(plot = fig3, filename = "04_figures/figures/bmRatio.pdf",
       height = 3, width = 6, units = "in")



### Species-Specific effects on conspecific density
con_re <- data.frame(ranef(bm.mod)) #extract random effects from model

con.sp_re <-filter(con_re, grpvar == "Seedling") 
con.sp_re <- con.sp_re[-(1:2)] #remove unneccesary columns

con.sp_re <- con.sp_re %>% mutate( lwr = condval - condsd) %>%
  mutate(upr = condval + condsd)%>%
  dplyr::rename(Sp = grp) #add columns for upper and lower CI values


