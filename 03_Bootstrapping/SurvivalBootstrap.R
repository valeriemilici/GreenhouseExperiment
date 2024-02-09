rm(list = ls())

library(lme4)
library(parallel)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(lme4))

surv.mod <- readRDS("02_models/modoutput/survival.RDS")
###Bootstrap survival model -----------------------------------------------

surv.fun <- function(.){
  preddat <- expand.grid(Treatment = c("D", "W"),
                         soil = c("conspecific", "heterospecific", "control"),
                         Init_Height = 7,
                         Cohort = 1)
  predict(., newdata = preddat, re.form = ~0)
}

surv.boot <- bootMer(surv.mod, nsim = 10000, FUN = surv.fun,
                   parallel = "snow", ncpus = detectCores(),
                   cl = cl)
save(surv.boot, file = "03_Bootstrapping/BootOutput/surv.boot")
