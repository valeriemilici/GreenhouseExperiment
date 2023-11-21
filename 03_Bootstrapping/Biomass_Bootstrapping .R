rm(list = ls())

library(lme4)
library(parallel)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(lme4))

bm.mod <- readRDS("02_models/modoutput/biomass.RDS")
bmratio.mod <- readRDS("02_models/modoutput/bmratio.RDS")
###Bootstrap Total Biomass model -----------------------------------------------

bm.fun <- function(.){
  preddat <- expand.grid(Treatment = c("D", "W"),
                         soil = c("conspecific", "heterospecific", "control"),
                         duration = 212,
                         Init_Height = 7)
  predict(., newdata = preddat, re.form = ~0)
}

bm.boot <- bootMer(bm.mod, nsim = 10000, FUN = bm.fun,
                    parallel = "snow", ncpus = detectCores(),
                    cl = cl)
save(bm.boot, file = "03_Bootstrapping/BootOutput/bm.boot")

###Bootstrap Biomass Ratio Model -----------------------------------------------
ratio.fun <- function(.){
  preddat <- expand.grid(Treatment = c("D", "W"),
                         soil = c("conspecific", "heterospecific", "control"),
                         duration = 212,
                         Init_Height = 7,
                         Tot_Biomass = 0.30)
  predict(., newdata = preddat, re.form = ~0)
}

bmratio.boot <- bootMer(bmratio.mod, nsim = 10000, FUN = ratio.fun,
                   parallel = "snow", ncpus = detectCores(),
                   cl = cl)

save(bmratio.boot, file = "03_Bootstrapping/BootOutput/bmratio.boot")

###Estimate the CIs for NPSF ---------------------------------------------------
bmPSF.fun <- function(.){
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

# bootstrap the models ------------

bm.PSF <- bootMer(bm.mod, nsim = 10000, FUN = bmPSF.fun,
                   parallel = "snow", ncpus = detectCores(),
                   cl =cl)

save(bm.PSF, file = "03_Bootstrapping/BootOutput/bm.PSF")

### Stop the cluster ----------------------------
stopCluster(cl = cl)
