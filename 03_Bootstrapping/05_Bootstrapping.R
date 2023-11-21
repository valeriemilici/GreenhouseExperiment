### Model Bootstraps ###

### Initialize Workspace -------------------------------------------------------
library(lme4)
library(parallel)

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(lme4))

### Bootstraps -----------------------------------------------------------------

### Wet and Dry RGR Models ----------------------

RGR.fun <- function(.){
  preddat <- expand.grid(soil = c("conspecific", "heterospecific", "control"),
                         Treatment = c("W", "D"),
                         Cohort = "2",
                         log_Init_Height_mm_s = 0)
  predict(., newdata=preddat, re.form = ~0)
}

RGRw.boot <- bootMer(wetRGR, nsim=1000, FUN = RGR.fun,
                     parallel = "snow", ncpus = detectCores(),
                     cl = cl)

RGRd.boot <- bootMer(dryRGR, nsim=1000, FUN = RGR.fun,
                     parallel = "snow", ncpus = detectCores(),
                     cl = cl)


### Height Model --------------------------------
height.fun <- function(.){
  preddat <- data.frame(cbind(Init_Height_mm <- c(60,60,60,60,60,60), duration.mo <- c(5,5,5,5,5,5), 
                              W <- c(5, 5, 0, 0,5,0), het <- c(5,0,5,0,0,0), cntrl <- c(0,5,0,5,0,0),
                              Whet <- c(5,0,0,0,0,0), Wcntrl <- c(0,5,0,0,0,0)))
  names(preddat)[1:7] <- c("Init_Height_mm", "duration.mo", "W", "het", "cntrl", "Whet", "Wcntrl")
  predict(., newdata = preddat, re.form = ~0)
}

height.mod.boot <- bootMer(height.mod, nsim = 1000, FUN = height.fun,
                           parallel = "snow", ncpus = detectCores(),
                           cl = cl)

save(height.mod.boot, file = "BootOutput/height.mod.boot")

### Survival Model ------------------------------
surv.fun <- function(.){
  preddat <- expand.grid(Treatment = c("W", "D"),
                         soil = c("conspecific", "heterospecific", "control"),
                         htprevcensus = 7)
  predict(., newdata = preddat, re.form = ~0)
}

surv.mod.boot <- bootMer(surv.ts, nsim = 1000, FUN = surv.fun,
                         parallel = "snow", ncpus = detectCores(),
                         cl =cl)

save(surv.mod.boot, file = "BootOutput/surv.mod.boot")

### Stop the cluster ----------------------------
stopCluster(cl = cl)
