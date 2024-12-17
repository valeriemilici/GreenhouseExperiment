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
library(ggeffects)
library(glmmTMB)
library(pls)
library(stringr)

dat <- read.csv("Data/mod_data_allpots.csv")
### Summary --------
# We assessed whether soil nutrient variation caused by the inoculum (20% of pot
# volume was field soil) could affect our primary response, total biomass. We
# used two methods to test this: first, by created a PCA of soil nutrient vars
# and using PC1 coordinates in the model of biomass. PC1 did not explain variation
# in total biomass (below). We also used a partial least squares regression 
# approach to identify specific soil nutrient variables that may be associated 
# with biomass, and then added those into the model. PLSR indicated that none of
# the nutrient variables were useful, but we included the top 2 vars in the model
# to be certain. The model indicated that the specific nutrients (mehlich Mg and
# mehlich Al) did not predict biomass. Models that include soil nutrients 
# tended to perform worse (assessed via AIC) than models that exclude nutrients,
# the exception is one model that included PC1. This model performs equally to
# the model that did not include nutrients. Nevertheless, soil nutrients to not
# improve model performance. 

### Data Prep --------
perf <- dat %>% filter(status == 1)#living plants only

#error on line 114- all Guapst Controls were harvested on day 212
perf[114,24] <- 212

## All treatments
perf$soil <- factor(perf$soil,
                        levels = c(  "heterospecific", "conspecific","control"))

# the raw soil moisture data
nutrients <- read.csv("Data/Soil_Geochem.csv")

# add nutrients to perf
nutrients1 <- nutrients %>% rename(Inoc_Tag = ID_2,
                                  Inoc_Sp = Sample.ID) %>% 
  group_by(Inoc_Sp) %>%
  mutate(NO2.NO3 = replace_na(NO2.NO3, mean(NO2.NO3, na.rm = T)),
         NH4 = replace_na(NH4, mean(NH4, na.rm = T)),
         pH_BaCl2 = replace_na(pH_BaCl2, mean(pH_BaCl2, na.rm = T)),
         Al = replace_na(Al, mean(Al, na.rm = T)),
         Ca = replace_na(Ca, mean(Ca, na.rm = T)),
         Fe = replace_na(Fe, mean(Fe, na.rm = T)),
         K = replace_na(K, mean(K, na.rm = T)),
         Mg = replace_na(Mg, mean(Mg, na.rm = T)),
         Mn = replace_na(Mn, mean(Mn, na.rm = T)),
         Na = replace_na(Na, mean(Na, na.rm = T)),
         TEB = replace_na(TEB, mean(TEB, na.rm = T)),
         ECEC = replace_na(ECEC, mean(ECEC, na.rm = T)),
         EBS = replace_na(EBS, mean(EBS, na.rm = T))) %>%
  ungroup() %>%
  mutate(Inoc_Tag = case_when(Inoc_Tag == "46965" ~ "59,46;1,1",
                              Inoc_Tag == "48187" ~ "57,43;4,3",
                              Inoc_Tag == "47905" ~ "57,55;3,3",
                              Inoc_Tag == "control" ~ "Control",
                              TRUE ~ Inoc_Tag)) %>%
  select(!c(Date, Lab.ID))
  
  
nutrients1[3,2] <- "20507_2"
nutrients1[15,2] <- "47183_2"
nutrients1[17,1] <- "PSIDFR"
nutrients1[21:23,1] <- "LACIAG"

nutrients1$Inoc_Sp <- toupper(nutrients1$Inoc_Sp)
nutrients1$m_B <- replace_na(as.numeric(nutrients1$m_B),0)
nutrients1$m_P <- replace_na(as.numeric(nutrients1$m_P), 0)

test <- left_join(perf, nutrients1, by = c("Inoc_Sp", "Inoc_Tag"))

test1 <- test %>% group_by(Inoc_Sp) %>%
  mutate(NO2.NO3 = replace_na(NO2.NO3, mean(NO2.NO3, na.rm = T)),
         NH4 = replace_na(NH4, mean(NH4, na.rm = T)),
         pH_BaCl2 = replace_na(pH_BaCl2, mean(pH_BaCl2, na.rm = T)),
         Al = replace_na(Al, mean(Al, na.rm = T)),
         Ca = replace_na(Ca, mean(Ca, na.rm = T)),
         Fe = replace_na(Fe, mean(Fe, na.rm = T)),
         K = replace_na(K, mean(K, na.rm = T)),
         Mg = replace_na(Mg, mean(Mg, na.rm = T)),
         Mn = replace_na(Mn, mean(Mn, na.rm = T)),
         Na = replace_na(Na, mean(Na, na.rm = T)),
         TEB = replace_na(TEB, mean(TEB, na.rm = T)),
         ECEC = replace_na(ECEC, mean(ECEC, na.rm = T)),
         EBS = replace_na(EBS, mean(EBS, na.rm = T)),
         LOI. = replace_na(LOI., mean(LOI., na.rm = T)),
         pH_H2O = replace_na(pH_H2O, mean(pH_H2O, na.rm = T)),
         pHCaCl2 = replace_na(pHCaCl2, mean(pHCaCl2, na.rm = T)),
         m_Al = replace_na(m_Al, mean(m_Al, na.rm = T)),
         m_B = replace_na(m_B, mean(m_B, na.rm = T)),
         m_Ca = replace_na(m_Ca, mean(m_Ca, na.rm = T)),
         m_Cu = replace_na(m_Cu, mean(m_Cu, na.rm = T)),
         m_Fe = replace_na(m_Fe, mean(m_Fe, na.rm = T)),
         m_K = replace_na(m_K, mean(m_K, na.rm = T)),
         m_Mg = replace_na(m_Mg, mean(m_Mg, na.rm = T)),
         m_Mn = replace_na(m_Mn, mean(m_Mn, na.rm = T)),
         m_P = replace_na(m_P, mean(m_P, na.rm = T)),
         m_Zn = replace_na(m_Zn, mean(m_Zn, na.rm = T))) %>%
  ungroup() 

  colSums(is.na(test1)) #Fab
  
  #remove the mis-identified inoculum source
  
  perf <- perf %>% filter(Inoc_Tag != "22401")
  
### Method 1: PCA --------------
  
# The original biomass model
bm.mod1 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + duration +
                           scale(Init_Height) +
                           (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin), 
                         data = perf)

# With both PCA options
bm.mod2 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + scale(duration) +
                            scale(Init_Height) + pca2 +
                            (1|Seedling) + (1|Inoc_Sp/Inoc_Tag) + (1|Table) + (1|begin), 
                          data = perf)

bm.mod3 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + scale(duration) +
                            scale(Init_Height) + pca2 +
                            (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin), 
                          data = perf)

AIC(bm.mod2, bm.mod3, bm.mod1)
# all variants are very close to one another. Great sign. The pca that has the 
# slightly lower aic (dAIC = ~0.8) is the pca that excludes nitrogen. Makes some
# sense because this PC1 explains 42.5% variation vs the other one explaines
# 39.8% variation. So slightly better. 

Anova(bm.mod3) #pca does not explain anything in the model
Anova(bm.mod1)
# otherwise the two models are ~ identical 
# again, they're identical

summary(bm.mod2)
summary(bm.mod3)

### Method 2: Partial Least Squares Regression --------------------------------

test2 <- test1 %>% select(22,31:56) #only test soil variables

pls.model = plsr(Tot_Biomass ~ ., data = test2, validation = "CV")

ncomp.onesigma <- selectNcomp(pls.model, method = "randomization",
                              plot = TRUE, ylim = c(.25, .35))
#WOW. The PLSR is suggesting that none of the soil variables should
# be included in the model. The "selection" is at zero.

# Find the number of dimensions with lowest cross validation error
cv = RMSEP(pls.model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1

# Rerun the model
pls.model = plsr(Tot_Biomass ~ ., data = test2, ncomp = best.dims)

#extract the identities of the coefficients that are the best
coefficients = coef(pls.model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(tail(coefficients, 10)) #they appear in reverse order
#m_Mg and m_Al score the highest. But the ncomp score also suggested that
# none of these vars are worthwhile. Let's put m_Mg in a model to test. 

bm.mod <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + duration +
                           scale(Init_Height) + scale(m_Mg) + scale(m_Al) + 
                           (1 + soil:Treatment||Seedling) + + (1|Inoc_Tag) +
                           (1|Table) + (1|begin), 
                         data = test1)

bm.mod2 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + duration +
                           scale(Init_Height) + scale(m_Mg) + scale(m_Al) + 
                           (soil|Seedling) + (1|Inoc_Tag) +
                           (1|Table) + (1|begin), 
                         data = test1)

bm.mod3 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + duration +
                            scale(Init_Height)  + 
                            (soil|Seedling) + (1|Inoc_Tag) +
                            (1|Table) + (1|begin), 
                          data = test1)

bm.mod4 <- lmerTest::lmer(log(Tot_Biomass) ~ Treatment * soil + duration +
                            scale(Init_Height)  + 
                            (1|Seedling) + (1|Inoc_Tag) +
                            (1|Table) + (1|begin), 
                          data = test1)
summary(bm.mod4)

AIC(bm.mod, bm.mod2, bm.mod3, bm.mod4)
# the simplest model is the best model. 



