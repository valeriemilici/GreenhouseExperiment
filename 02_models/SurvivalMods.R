### Survival Models ###

### Initialize Workspace -------------------------------------------------------

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
library(effects)

dat.ts <- read.csv("Data/TimeSeriesData.csv")
dat <- read.csv("Data/mod_data_allpots.csv")

### Data Manipulation ----------------------------------------------------------

##step 1: find the troublesome rows and change them so they're appropriate
dat.check <- subset(dat.ts, htprevcensus != 999, CensusNo !=1)
dat.check<- subset(dat.check, Height == 0)
dat.check$Height = 999
dat.check$status = 1

##step 2: remove those rows from dat.ts
dat.ts <- dat.ts[-c(608, 1745, 2007, 2224, 2274, 2418),]

##step 3: add in the correct rows
dat.ts <- rbind(dat.ts, dat.check)

## standardize height from prev census
dat.ts.s <- dat.ts %>% mutate(
  across(htprevcensus, log, .names = "log_{.col}"))

dat.s <- filter(dat.ts.s, log_htprevcensus >= 0) #removes NAs and inf

dat.s$log_htprevcensus_s <- ((dat.s$log_htprevcensus -
                      mean(dat.s$log_htprevcensus))/sd(dat.s$log_htprevcensus))

#Add Cohort Column

cohort <- select(dat, c("Location", "Cohort"))

dat.s <- left_join(dat.s, cohort, by = "Location")

#Add live vs. sterile soil column
dat <- dat %>% mutate(microbes = if_else(soil == "control", "sterile", "live"))

### Times Series Survival Model ------------------------------------------------

#We don't report the time series data,
#so consider removing this code from the script
surv.ts <- glmer(status ~ soil*Treatment + log_htprevcensus_s + log(Duration) +
                    (1|Seedling)+ (1|Inoc_Sp) + (1|Location) + (1|Begin),
                  data = subset(dat.s, htprevcensus != 999 &
                                  htprevcensus != 0 &
                                  CensusNo != 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)),
                  offset = log(censusint))
#note- Duration is a much better fit than Cohort (dAIC = 12)
summary(surv.ts)
#this tells me that in the dry treatment, sterile control plants were 22% less likely
#to die each day than plants in live soil. This pattern changes in the wet treatment
# control plants then are ~just as likely to die as plants in live soil treatments. 
#so microbes improve survival outcomes in lower moisture but don't affect survival in 
# higher moisture?
Anova(surv.ts)
## Diagnostics
testDispersion(surv.ts) #good
modelOutput <- simulateResiduals(surv.ts, plot=F)
plot(modelOutput) #looks good

testdat <- subset(dat.s, htprevcensus != 999 & htprevcensus != 0 &
                    !is.na(htprevcensus) & CensusNo != 1)

plotResiduals(modelOutput, testdat$soil)

plot(allEffects(surv.ts))

### Basic Survival Model -------------------------------------------------------
dat$soil <- factor(dat$soil, levels =c("control", "conspecific", "heterospecific"))

surv <- glmer(status ~ soil*Treatment + scale(Init_Height) + as.factor(Cohort) +
                 (1|Seedling)+ (1|Inoc_Sp) + (1|Table) + (1|begin),
               family = binomial(link = "logit"),
               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)),
               data = dat)
testDispersion(surv)
survOutput <- simulateResiduals(surv, plot=F)
plot(survOutput) #looks way better with these changes. Cohort is a much better
#fit than duration here. 

summary(surv)

plotResiduals(survOutput, dat$soil)
