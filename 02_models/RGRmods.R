### RGR Models ###

### Initialize Workspace -------------------------------------------------------

rm(list = ls())

library(tidyverse) #always
library(lme4) #runs models
library(lmerTest) #model testing, gives p values
library(performance) #diagnostics
library(DHARMa) #additional diagnostics

dat <- read.csv("Data/mod_data_allpots.csv")


### Data Manipulation ----------------------------------------------------------
perf_dat <- filter(dat,status == 1)

#Useful Function --------------------------------

RGRf <- function(h1, h2, t){
  h1_l = log(h1)
  h2_l = log(h2)
  RGR = (h2_l - h1_l)/t 
  return(RGR)
}

#Add RGR columns --------------------------------
perf_dat3 <- perf_dat %>% 
  mutate(duration.mo = duration/30,
         Final_Height_mm = Final_Height*10,
         Init_Height_mm = Init_Height*10,
         RGR = RGRf(Init_Height_mm, Final_Height_mm, duration.mo),
         log_Init_Height_mm = log(Init_Height_mm),
        wanted = ifelse(soil == "control", "dead", "live"),
        RGRn = (RGR - mean(RGR)) / sd(RGR)
) 

### RGR ------------------------------------------------------------------

#relevel treatment contrasts
perf_dat3$soil <- factor(perf_dat3$soil,
                          levels = c("control", "conspecific", "heterospecific"))

perf_dat3$soil <- factor(perf_dat3$soil,
                         levels = c("conspecific", "heterospecific", "control"))


# use for effects and forest plot
RGR.mod<- lmer(RGR ~ soil*Treatment + 
                 log_Init_Height_mm + 
                 duration + 
                  (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
               data = perf_dat3)

#It is better to use duration rather than cohort to capture differences based on
#when treatments began/were harvested. 

check_model(RGR.mod)
anova(RGR.mod)
summary(RGR.mod)
saveRDS(RGR.mod, file = "02_models/modoutput/RGRmod.RDS")

## compare observed to the fitted values --------------------------------------

perf_dat3$fit <- fitted(RGR.mod2)

dHt <- function(RGR, H1, t){
  Ht.change = exp(RGR*t)*H1 - H1
  
  return(Ht.change) 
  }

perf_dat3 <- mutate(perf_dat3, Ht.ch = dHt(fit, Init_Height_mm, duration.mo))
plot(perf_dat3$dH, perf_dat3$Ht.ch)
#the model isn't doing a good job with the extreme observations
#thinning the data doesn't really change much with performance. 


# Patterns when ecto-mycorrhizal species are removed ---------------------------
RGR.mod2<- lmer(RGR ~ soil*Treatment + 
                 log_Init_Height_mm + 
                 duration + 
                 (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
               data = perf_dat4)

check_model(RGR.mod2)
anova(RGR.mod2)
summary(RGR.mod2)

### Sum to Zero Contrast RGR ---------------------------------------------------

#set contrasts
perf_dat3$soil <- as.factor(perf_dat3$soil)
#perf_dat3$Treatment <- as.factor(perf_dat3$Treatment)
contrasts(perf_dat3$soil) <- contr.sum(3)
#contrasts(perf_dat3$Treatment) <- contr.sum(2)

RGR.mod3<- lmer(RGR ~ log_Init_Height_mm + 
                  soil*Treatment + as.factor(Cohort) +
                  (1|Seedling) + (1|Inoc_Sp) + (1|Table),
                data = perf_dat3)


summary(RGR.mod3) 
#soil1 is conspecific, soil 2 is control


saveRDS(RGR.mod3, file = "modoutput/RGRcc.RDS")

#get the het estimate----------------------------
perf_dat3$soil <- factor(perf_dat3$soil,
                         levels = c("conspecific", "heterospecific", "control"))
contrasts(perf_dat3$soil) <- contr.sum(3)

RGR.mod4<- lmer(RGR ~ log_Init_Height_mm + 
                  soil*Treatment + as.factor(Cohort) +
                  (1|Seedling) + (1|Inoc_Sp) + (1|Table),
                data = perf_dat3)


summary(RGR.mod4)

saveRDS(RGR.mod4, file = "modoutput/RGRhet.RDS")