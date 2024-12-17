### Species Specific Effects
### Initializing Workspace -----------------------------------------------------
rm(list = ls())
library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(ggeffects)
library(stringr)
library(brms)
library(broom.mixed)
library(tidybayes)
library(bayesplot)

### Read in and Prep Data ------
dat <- read.csv("data/mod_data_allpots.csv")

#summary(dat)

dat2 <- mutate(dat, Tot_Biomass_0 = replace_na(Tot_Biomass, 0),
               soil = relevel(factor(soil), "heterospecific"),
               Cohort = factor(Cohort))
## There are 2 PSIDFR that have NA initial heights. Replacing with mean
## PSIDFR intial height

psidfr_ht <- mean(dat2$Init_Height[dat2$Seedling == "PSIDFR"], na.rm=T)
dat2 <- mutate(dat2, Init_Height = replace_na(Init_Height, psidfr_ht))
#summary(dat2)
write.csv(dat2, "Output/Models/dat2.csv") #save for plotting

### The Hurdle Model (Survival and Biomass effects combined) ----------

prior_bm_mod <- c(prior(student_t(7, 0, 2), class = b),
                  prior(student_t(7, 0, 2), class = b,
                        dpar = "hu"))# set t priors with 7 df,
# zero mean and sd of 2 on all population-level parameters

## The Model ---
brms_bm.mod <- brm(bf(Tot_Biomass_0 ~ scale(duration) +
                        scale(Init_Height) + soil.dim1 + 
                        Treatment * soil + 
                        (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin),
                      hu ~  Cohort + scale(Init_Height) + soil.dim1 +
                        Treatment * soil +
                        (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin)),
                   data = dat2, cores=4, 
                   family = hurdle_lognormal(link_hu="logit"),
                   control=list(adapt_delta = 0.99, max_treedepth = 12), 
                   prior=prior_bm_mod)
## diagnostics
pp_check(brms_bm.mod, ndraws=100) +
  pp_check(brms_bm.mod, resp = "mu", ndraws=100) +
  pp_check(brms_bm.mod, resp="hu", ndraws=100) ## look reasonable

summary(brms_bm.mod)
tidy(brms_bm.mod)|> print(n = 18)

## Save the Model
save(brms_bm.mod, file = "Output/Models/brms_bm")

### Biomass Ratio Model ------

prior_bmrat_mod <- prior(student_t(7, 0, 2), class = b) # set t priors with 7 df,
# zero mean and sd of 2 on all population-level parameters

#root mass fraction ---
brms_RMF.mod <- brm(log(Root_Mass/Tot_Biomass) ~ scale(duration) +
                        scale(Init_Height) + soil.dim1 + 
                        Treatment * soil + 
                        (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin),
                      data = filter(dat2, status == 1), cores=4, 
                      family = gaussian(),
                      control=list(adapt_delta = 0.99), 
                      prior=prior_bmrat_mod)

summary(brms_RMF.mod) #this model is easier to interpet.
#although there are some significant differences among the treatments, overall
#RMF among treatments are very similar (within 4% of each other).

pp_check(brms_RMF.mod, ndraws=100) +
  pp_check(brms_RMF.mod, resp = hu, ndraws=100) ## good

#Save model
save(brms_RMF.mod, file = "Output/Models/brms_RMF.mod")



### Species Specific Model -----

prior_bm_mod_sp <- c(prior(student_t(7, 0, 2), class = b),
                     
                     prior(student_t(7, 0, 2), class = b, dpar = "hu"))
brms_bm.mod_sp <- brm(bf(Tot_Biomass_0 ~ Treatment * soil + scale(duration) +
                           scale(Init_Height) +  soil.dim1 +
                           (1 + soil*Treatment||Seedling) +
                        (1|Inoc_Tag) + (1|Table) + (1|begin),
                      hu ~ Treatment * soil + Cohort + 
                        scale(Init_Height) + soil.dim1 +
                        (1 + soil*Treatment||Seedling) +
                        (1|Inoc_Tag) + (1|Table) + (1|begin)),
                      data = dat2, cores=4, family = hurdle_lognormal(),
                      control=list(adapt_delta = 0.99, max_treedepth = 12),
                      prior=prior_bm_mod_sp)

pp_check(brms_bm.mod_sp, ndraws=100) +
  pp_check(brms_bm.mod_sp, resp = hu, ndraws=100) ## not bad

summary(brms_bm.mod_sp)

save(brms_bm.mod_sp, file = "Output/Models/brms_bm.mod_sp")