## RGR Figs###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())
library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(sjPlot) #forest plots
library(ggthemes)

RGR.mod <- readRDS("modoutput/RGRmod.RDS") #specialized treeplot
RGR.mod2 <- readRDS("modoutput/RGRwanted.RDS") #generalized treeplot
concntrl <- readRDS("modoutput/RGRcc.RDS") #spec estimate fig
het<- readRDS("modoutput/RGRhet.RDS") #spec estimate fig
live <- readRDS("modoutput/RGRlive.RDS") #live ests
dead <- readRDS("modoutput/RGRdead.RDS") #dead ests

### Plotting Specialized Effects (NSPF) ---------------------------------------
### Tree Plot -----------------------------------------------------------------

RGR.forest <- plot_model(RGR.mod,
                         rm.terms = c("log_Init_Height_mm",
                                      "as.factor(Cohort)2"),
                         vline.color = "grey",
                         title = "",
                         axis.labels = c("  Heterospecific Soil x\n100% Soil Moisture",
                                         "   Conspecific Soil x\n100% Soil Moisture",
                                         "100% Soil Moisture", 
                                         "Heterospecific Soil",
                                         "Conspecific Soil"),
                         axis.title = "Treatment effects on RGR log(mm/mm/mo)",
                         show.values = T) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.05, 0.06)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))


#summary(RGR.mod)
RGR.forest

ggsave("04_figures/RGR.forest.jpeg", width = 6, height = 4, units = "in")

### RGR model Sum to Zero Contrasts### -----------------------------------------
### Get the SEs without intercept ---------------
preddat <- expand.grid(soil = c( "heterospecific", "control", "conspecific"),
                       Treatment = c("D", "W"))

hetdat <- expand.grid(soil = c(  "control", "conspecific", "heterospecific"),
                      Treatment = c("D", "W"))
#preddat$soil <- factor(preddat$soil,
 #                      levels = c("control", "conspecific", "heterospecific"))

#variance covariance matrix of the model without nuisance variables
mod.vcov <- vcov(concntrl)[-c(1:2,6),][,-c(1:2,6)]
het.vcov <- vcov(het)[-c(1:2,6),][,-c(1:2,6)]
#create the design matrix
Xmat <- model.matrix(~soil*Treatment, data = preddat)
hetmat <- model.matrix(~soil*Treatment, data = hetdat)
XB <- model.matrix(~soil*Treatment, data = preddat)
Xmat<- Xmat[-c(1,4),][,-1] #remove intercept
hetmat <- hetmat[-c(1,4),][,-1]
XBeta <- XB %*% fixef(concntrl)[-c(2,6)] #1 = het, 2 = consp, 3 = cntrl, 4 = whet, 5 = wconsp, 6 = wcntrl

#matrix algebra between the two
concntrl.se <- data.frame(se = sqrt(diag(Xmat %*% mod.vcov %*% t(Xmat))))
het.se <- data.frame(se = sqrt(diag(hetmat %*% het.vcov %*% t(hetmat))))


dry.est <- data.frame(var = c("Heterospecific", "Conspecific", "Control"),
                      est = XBeta[c(1:3),],
                      trt = c("D", "D", "D"),
                      se = c(het.se[2,],concntrl.se[1,], concntrl.se[2,]))

wet.est <- data.frame(var = c("Heterospecific", "Conspecific", "Control"),
                      est = XBeta[c(4:6),],
                      trt = c("W", "W", "W"),
                      se = c(het.se[4,],concntrl.se[3,], concntrl.se[4,]))

mod.est <- rbind(wet.est, dry.est)

mod.est <- mutate(mod.est, upr = est + se, lwr = est - se)

#relevel so they plot in the best order
mod.est$var <- factor(mod.est$var,
                      levels = c("Control", "Conspecific", "Heterospecific"))

#create descriptive labels for the facts
trt.labs <- as_labeller(c("D" = "50% Soil Moisture",
                          "W" = "100% Soil Moisture"))

#plot the output --------------------------------
RGR.final <- ggplot(mod.est, aes(var, est)) + 
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.4)) +
  facet_wrap(~trt,
             labeller = labeller(trt = trt.labs)) +
  labs(x = "Soil Type",
       y = "RGR \nlog(mm/mm/mo)",
       col = "") +
  theme_bw()+
  ylim(0.419, 0.45) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(angle = 15, hjust = 0.7))

RGR.final #looks good

#save output
ggsave("04_figures/RGRfinal.jpeg", width = 6, height = 4, units = "in") 


## Combined plot of estimates and predictions ---------------------------------

 spec.fig <- RGR.forest + RGR.final +
  plot_annotation(tag_levels = "A")

spec.fig 

ggsave("04_figures/RGRspec.jpeg", width = 12, height = 5, units = "in")
### RGR Figure: Generalized effects -------------------------------------------
### Get the SEs without intercept ---------------
preddat <- expand.grid(wanted = c("dead", "live"),
                       Treatment = c("D", "W"))

#variance covariance matrix of the model without nuisance variables
live.vcov <- vcov(live)[-c(1:2,5),][,-c(1:2,5)] 
dead.vcov <- vcov(dead)[-c(1:2,5),][,-c(1:2,5)] 

#create the design matrix
Xmat <- model.matrix(~wanted*Treatment, data = preddat)
deadmat <- model.matrix(~wanted*Treatment, data = preddat)
XB <- model.matrix(~wanted*Treatment, data = preddat)
Xmat<- Xmat[-1,][,-1] #remove intercept
deadmat <- deadmat[-1,][,-1]
XBeta <- XB %*% fixef(live)[-c(2,5)] #dead, live, wdead, wlive

#matrix algebra between the two
live.se <- data.frame(se = sqrt(diag(Xmat %*% live.vcov %*% t(Xmat))))
dead.se <- data.frame(se = sqrt(diag(deadmat %*% dead.vcov %*% t(deadmat))))
est <- data.frame(var = c("Sterile", "Live", "Sterile", "Live"),
                  est = XBeta[c(1:4),],
                  trt = c("D", "D", "W", "W"),
                  se = c(dead.se[1,], live.se[1,], dead.se[3,], live.se[3,]))

est <- est %>% mutate(lwr = est - se,
                      upr = est + se)


#Plot ----------------------------
#create descriptive labels for the facts
trt.labs <- as_labeller(c("D" = "50% Soil Moisture",
                          "W" = "100% Soil Moisture"))

est$var <- factor(est$var,
                  levels = c("Sterile", "Live"))

#plot the output --------------------------------
RGR.general <- ggplot(est, aes(var, est)) + 
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(~trt,
             labeller = labeller(trt = trt.labs)) +
  labs(x = "Soil Type",
       y = "RGR \nlog(mm/mm/mo)",
       col = "") +
 # scale_x_discrete(expand = c(1,0.05)) +
  theme_bw()+
  ylim(0.42, 0.45) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 12))

ggsave("04_figures/RGRgeneral.jpeg", height = 3, width = 4.5, units = "in")

### Forest Plot ---------------------------------

gen.forest <- plot_model(RGR.mod2,
                         rm.terms = c("log_Init_Height_mm",
                                      "as.factor(Cohort)2"),
                         vline.color = "grey",
                         title = "",
                         axis.labels = c("          Live Soil x \n100% Soil Moisture",
                                         "100% Soil Moisture", 
                                         "Live Soil"),
                         axis.title = "Treatment effects on RGR log(mm/mm/mo)",
                         show.values = T) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.05, 0.06)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

ggsave("04_figures/gen.forest.jpeg", height = 4, width = 6, units = "in")

### Patchwork Plot ---------------------------------

gen.forest + RGR.general + plot_annotation(tag_levels = "A")

ggsave("04_figures/genplot.jpeg", height = 5, width = 12, units = "in")
