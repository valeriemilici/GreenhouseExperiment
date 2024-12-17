### Height Change Figures ###

### Initialize Workspace -------------------------------------------------------
rm(list = ls())

library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(sjPlot) #forest plots

### Model Output --------------------------------------------------------------
htmod <- readRDS("modoutput/htmod.RDS")

### Forest Plot ---------------------------------------------------------------

Ht.forest <- plot_model(htmod,
                         rm.terms = c("log_Init_Height_mm",
                                      "as.factor(Cohort)2",
                                      "as.factor(Cohort)1",
                                      "t_S0",
                                      "duration.mo"),
                         vline.color = "grey",
                         #order.terms = c(2,1,3,5,4),
                         title = "",
                         # axis.labels = c("Sterile Soil x Wet", "Heterospecific Soil x Wet",
                         #               "Wet Treatment",  "Sterile Soil",
                         #              "Heterospecific Soil" ),
                         axis.title = "Effect on Growth (mm)",
                         show.values = T,
                        #transform = exp) + work on this bit
                         theme_bw())

Ht.forest

### Height Change -------------------------------------------------------------

preddat <- data.frame(t_S0 =c(21,21,21,21,21,21),
                          duration.mo = c(5,5,5,5,5,5), 
                          het = c(5,0,5,0,0,0),
                          con = c(0,5,0,5,0,0),
                          W = c(5, 5, 0, 0,5,0),
                          Whet = c(5,0,0,0,0,0), 
                          Wcon = c(0,5,0,0,0,0))


#variance covariance matrix of the model without nuisance variables
mod.vcov <- vcov(height.mod)[-c(3:4),][,-c(3:4)]

#create the design matrix
#for the SE matrix multiplication
Xmat <- model.matrix(~0 + t_S0 + duration.mo + het + con + W + Whet + Wcon , data = preddat)

#for the main effects
XB <- model.matrix(~0 + t_S0 + duration.mo + het + con + W + Whet + Wcon , data = preddat)
XBeta <- XB %*% fixef(height.mod)[-c(3:4)] 

#matrix algebra between the two
est.se <- data.frame(se = sqrt(diag(Xmat %*% mod.vcov %*% t(Xmat))))

dry.est <- data.frame(var = c("heterospecific" ,"conspecific"),
                      est = XBeta[c(3:4),],
                      trt = c("D", "D"),
                      se = est.se[c(3:4),])
wet.est <- data.frame(var = c("heterospecific" ,"conspecific"),
                      est = XBeta[c(6:7),],
                      trt = c("W", "W"),
                      se = est.se[c(6:7),])


dry.est <- dry.est %>% mutate(lwr = est - se,
                              upr = est + se)

wet.est <- wet.est %>% mutate(lwr = est - se,
                              upr = est + se)
#Dry Treatment -----------------------------------------------------------------
dryplotH <- ggplot(dry.est, aes(trt, exp(est), col = var, group = var)) + 
  geom_pointrange(aes(ymin = exp(lwr), ymax = exp(upr)),
                  position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = exp(XBeta[1,1]), color = " dark grey") +
  labs(x = "",
       y = "Height Change (mm)",
       col = "") +
  scale_x_discrete(labels = c( "D" = "50% Saturation"),
                   limits = c("D")) +
  scale_fill_viridis_d(option = "viridis", begin = 0.1, end = 0.6) +
  scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.6) +
  theme_classic()+
  #ylim(0.415, 0.460) +
  theme(legend.position = 'bottom') +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 16))