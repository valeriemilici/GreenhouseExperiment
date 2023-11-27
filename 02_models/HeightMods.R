### Change in Height Models ###

## Notes- although there are some intriguing patterns between the height and 
## root models, the patterns aren't consistent or strong enough to add clarity
## to the story. They show some differences that indicate that conspecific 
## microbes may have a negative effect on plant growth relative to control and
## heterospecific plants, and although that is promising, I think it would not
## add enough to the story. 

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

dat <- read.csv("Data/mod_data_allpots.csv")

### Data Manipulation ----------------------------------------------------------
perf_dat <- filter(dat,status == 1) 
perf_dat <- mutate(perf_dat, Wanted = ifelse(soil == "control", "dead", "live"))

### Manual Matrix -------------------------------

perf_dat2 <- perf_dat
perf_dat2 <- mutate(perf_dat2, duration.mo = duration/30)


perf_dat2$W <- ifelse(perf_dat2$Treatment == "W", perf_dat2$duration.mo , 0)
perf_dat2$het <- ifelse(perf_dat2$soil == "heterospecific",perf_dat2$duration.mo , 0)
perf_dat2$con <- ifelse(perf_dat2$soil == "conspecific", perf_dat2$duration.mo, 0)

perf_dat2$Whet <- ifelse(perf_dat2$soil == "heterospecific", perf_dat2$W, 0)

perf_dat2$Wcon <- ifelse(perf_dat2$soil == "conspecific", perf_dat2$W, 0)

perf_dat2 <- mutate(perf_dat2, Final_Height_mm = Final_Height*10)
perf_dat2 <- mutate(perf_dat2, Init_Height_mm = Init_Height*10)

perf_dat2 <- mutate(perf_dat2, t_S0 = duration.mo * log(Init_Height_mm))

### Change in Height -----------------------------------------------------------
height.mod <- lmer(log(Final_Height_mm) ~ 0 + t_S0 + duration.mo +
                      W + het + con + Whet  + Wcon +
                      (0 + duration.mo | Seedling) + 
                      (0 + duration.mo | Inoc_Sp) + 
                      (0 + duration.mo |Table) +
                     offset(log(Init_Height_mm)),
                    data = perf_dat2)
summary(height.mod)

plot_model(height.mod, rm.terms = c("log(Init_Height_mm)")) + coord_flip(ylim = c(-0.05, 0.05)) + scale_y_continuous(breaks=seq(-0.05, 0.05, 0.025))

check_model(height.mod) #looks good
saveRDS(height.mod, "modoutput/htmod.RDS") 

mod1 <- lmer(log(Final_Height) ~ log(Init_Height) + soil*Treatment +
               as.factor(Cohort) + 
               (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
             data = perf_dat)

anova(mod1)
check_model(mod1)
#looks pretty good
#note- cohort is better than duration (dAIC = 11)
summary(mod1)

### Root Length ----------------------------------------------------------------

#Root Length
Root.mod1<- lmer(scale(Root_LG) ~ 0 + scale(Init_Height) + 
                   soil*Treatment + as.factor(Cohort) +
                   (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
                 data = perf_dat)
check_model(Root.mod1)
summary(Root.mod1) #absolutely no effect
anova(Root.mod1)
Root.mod2<- lmer(scale(Root_LG) ~ scale(Init_Height) + 
                   soil*Treatment + duration +
                   (1|Seedling) + (1|Inoc_Sp) + (1|Table) + (1|begin),
                 data = perf_dat)
check_model(Root.mod1)
summary(Root.mod2)

AIC(Root.mod1, Root.mod2)
