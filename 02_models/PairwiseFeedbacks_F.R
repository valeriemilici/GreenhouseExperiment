###Pairwise Biomass --------------
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

dat <- read.csv("Data/mod_data_allpots.csv")

perf_dat <- dat %>% filter(status == 1) %>%
  mutate(Wanted = ifelse(soil == "control", "dead", "live"))

#error on line 114- all Guapst Controls were harvested on day 212
perf_dat[114,24] <- 212

### Guapst and HEISCO ----- This is the only one we can model this for
dat2 <- perf_dat %>% filter(Seedling == c("GUAPST" , "HEISCO")) %>%
  filter(Inoc_Sp == c("GUAPST", "HEISCO")) %>%
           filter(soil != "control")

m1 <- lmer(log(Tot_Biomass) ~ Seedling + Inoc_Sp * Treatment + log(Init_Height) + (1|Table),
           data = dat2)
summary(m1)

library(lme4)
library(parallel)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(lme4))

PSF.fun <- function(.){
  Gg.dat <- expand.grid(Seedling = "GUAPST",
                         Inoc_Sp = "GUAPST",
                        Treatment = c("D", "W"),
                         Init_Height = 7)
  Gh.dat <- expand.grid(Seedling = "GUAPST",
                        Inoc_Sp = "HEISCO",
                        Treatment = c("D", "W"),
                        Init_Height = 7)
  Hg.dat <- expand.grid(Seedling = "HEISCO",
                        Inoc_Sp = "GUAPST",
                        Treatment = c("D", "W"),
                        Init_Height = 7)
  Hh.dat <- expand.grid(Seedling = "HEISCO",
                        Inoc_Sp = "HEISCO",
                        Treatment = c("D", "W"),
                        Init_Height = 7)
  
  Gg.est <- predict(., newdata = Gg.dat, re.form = ~0)
  Gh.est <- predict(., newdata = Gh.dat, re.form = ~0)
  Hg.est <- predict(., newdata = Hg.dat, re.form = ~0)
  Hh.est <- predict(., newdata = Hh.dat, re.form = ~0)
  
  log(exp(Gg.est)) - log(exp(Gh.est)) - log(exp(Hg.est)) + log(exp(Hh.est))
}

# bootstrap the models ------------

PSF <- bootMer(m1, nsim = 1000, FUN = PSF.fun,
                  parallel = "snow", ncpus = detectCores(),
                  cl =cl)
### Stop the cluster ----------------------------
stopCluster(cl = cl)

confint(PSF)



