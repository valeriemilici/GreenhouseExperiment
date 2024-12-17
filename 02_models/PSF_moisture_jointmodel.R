### Species Specific Effects
### Initializing Workspace -----------------------------------------------------

library(tidyverse) #always
library(ggplot2) #for plotting
library(patchwork) #joins plots nicely
library(stringr)
library(brms)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


### Read in and Prep Data ------
dat <- read.csv("mod_data_allpots.csv")


summary(dat)
### Set Contrasts -------
dat2 <- mutate(dat, Tot_Biomass_0 = replace_na(Tot_Biomass, 0),
               soil = relevel(factor(soil), "heterospecific"),
               Cohort = factor(Cohort))
## There are 2 PSIDFR that have NA initial heights. Replacing with mean
## PSIDFR intial height
## 
psidfr_ht <- mean(dat2$Init_Height[dat2$Seedling == "PSIDFR"], na.rm=T)
dat2 <- mutate(dat2, Init_Height = replace_na(Init_Height, psidfr_ht))
summary(dat2)



## Fit the population-level model
prior_bm_mod <- prior(student_t(7, 0, 2), class = b) # set t priors with 7 df,
# zero mean and sd of 2 on all population-level parameters

with(dat2, table(duration,  Cohort))

brms_bm.mod <- brm(bf(Tot_Biomass_0 ~ scale(duration) +
                        scale(Init_Height) + soil.dim1 + 
                        Treatment * soil + 
                        (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin),
                      hu ~  Cohort + scale(Init_Height) + soil.dim1 +
                        Treatment * soil +
                        (1|Seedling) + (1|Inoc_Tag) + (1|Table) + (1|begin)),
                   data = dat2, cores=4, 
                   family = hurdle_lognormal(link_hu="logit"),
                   control=list(adapt_delta = 0.99), 
                   prior=prior_bm_mod)
                   
## diagnostics
pp_check(brms_bm.mod, ndraws=100) +
  pp_check(brms_bm.mod, resp = "mu", ndraws=100) +
  pp_check(brms_bm.mod, resp="hu", ndraws=100) ## look reasonable

summary(brms_bm.mod)
tidy(brms_bm.mod)|> print(n = 18)

pred_dat_pop <- with(dat2,
                 expand_grid(Treatment = unique(Treatment),
                        soil = levels(soil), 
                        soil.dim1 = 0, duration = 212, Init_Height = 7,
                        Cohort = levels(Cohort)[1]))
## get posterior draws of biomass (mu) and mortality (hu)
pred_dat_pop <- add_epred_draws(brms_bm.mod, newdata = pred_dat_pop, 
                             dpar = c("mu", "hu"),
                             re_formula = ~ 0)

## quick plot
pred_dat_pop |> ggplot(aes(x = soil, y = .epred)) +
  facet_wrap(~Treatment) + stat_pointinterval(point_interval = "median_hdci")

## calculate credible intervals
plot_dat_pop <- pred_dat_pop |> group_by(soil, Treatment) |> 
  median_hdci(mu, hu, .epred, .width=c(0.95)) |> 
  rename("mu.est" = "mu", "hu.est" = "hu", "epred.est" = ".epred", 
         "epred.lower" = ".epred.lower", "epred.upper" = ".epred.upper") |> 
  pivot_longer(cols=mu.est:epred.upper) |> 
  separate_wider_delim(name, delim =".", names=c("vital_rate", "stat")) |> 
  mutate(value = case_when(
    vital_rate == "hu" ~ 1-value,
    vital_rate == "mu" ~ exp(value),
    .default = value)) |> 
  pivot_wider(names_from=stat, values_from = value)

ggplot(plot_dat_pop,
       aes(x = soil, y=est, colour = Treatment,
           ymin = lower, ymax = upper)) +
  facet_wrap(~vital_rate) +
  geom_pointrange(position = position_dodge2(width=0.3)) +
  theme_bw()

pl_pop <- group_by(plot_dat_pop, vital_rate) |> 
  group_map(\(.x, ...) 
            ggplot(.x,
                   aes(x = soil, y=est, colour = Treatment,
                       ymin = lower, ymax = upper)) +
              geom_pointrange(position = position_dodge2(width=0.3)) +
              theme_bw())
((pl_pop[[1]] + labs(y = "Biomass combined")) | 
    ((pl_pop[[2]] + labs(y = "survival"))/(pl_pop[[3]] + labs(y = "biomass")))) +
  plot_layout(guides = "collect")

## error bars are a bit large, but this is probably partly due to the large
## uncertainty in the intercept. 

## calculating the log response ratio
plot_psf_pop <- ungroup(pred_dat_pop) |> select(-(soil.dim1:`.iteration`)) |> 
  mutate(.epred = log(.epred), hu = log(1-hu)) |> 
  rename("combined" = ".epred", "survival" = "hu", "biomass" = "mu") |> 
  pivot_longer(cols = combined:survival, names_to = "vital_rate") |> 
  pivot_wider(names_from=soil, values_from=value) |> 
  pivot_longer(cols=conspecific:control, names_to="soil") |> 
  mutate(psf = value - heterospecific,
         Treatment = factor(Treatment, labels = c("50% WHC", "100% WHC")),
         vital_rate = factor(vital_rate, 
                             levels = c("combined", "biomass", "survival")),
         soil = factor(ifelse(soil == "conspecific", "live", "sterile")))|> 
  group_by(Treatment, soil, vital_rate) |>
  median_hdci(psf, .width = c(0.66, 0.95))


## log response ratio results
pl_psf_pop <- group_by(plot_psf_pop, vital_rate) |> 
  group_map(\(.x, ...) 
            ggplot(.x,
                   aes(x = soil, y=psf, colour = Treatment,
                       ymin = .lower, ymax = .upper)) +
              geom_hline(yintercept=0, linetype = "dashed") +
              geom_pointrange(data = filter(.x, .width == 0.95),
                              position = position_dodge2(width = 0.2)) +
              geom_pointrange(data = filter(.x, .width == 0.66),linewidth = 1,
                              position = position_dodge2(width = 0.2)) +
              #labs(x = NULL, y = NULL) +
              labs(x = "Soil Inoculum Source", y = "Log Response Ratio") +
              theme_bw())

(((pl_psf_pop[[1]] + labs(tag = "A", x= NULL)) | 
    ((pl_psf_pop[[2]] + labs(tag = "B", y = NULL, x = NULL))/
       (pl_psf_pop[[3]] + labs(tag = "C", y = NULL, x = NULL)))) +
  plot_layout(guides = "collect"))  +
  plot_annotation(caption = "Soil Inoculum Source") &
  theme(plot.caption = element_text(hjust = 0.5, size = 12))

## plot the parameters 
## Might be a more striking way of presenting the model coefficients
draws <- tidy_draws(brms_bm.mod) |> 
  dplyr::select(starts_with("."), starts_with("b_"))
draws <- pivot_longer(draws, starts_with("b_"), names_to = "par") |> 
   mutate(par = str_sub(par, start = 3))

## identify mortality and growth parameters.
draws <- draws |> separate_wider_delim(par, delim="u_", names=c("par", "variable"),
                                       too_few = "align_end"  ) |>
  mutate(par = factor(ifelse(is.na(par), "biomass", "survival"),
                      levels = c("survival", "biomass")),
         value = value * ifelse(par == "survival", -1, 1)) ## convert motality to survival

tree_plot_pop <- filter(draws, !(variable  %in% c("Intercept"))) |> 
  group_by(par, variable) |> median_hdci(.width = c(.95)) |>
  ggplot(aes(y = variable, x = value, xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept=0, linetype = "dashed") +
  geom_pointinterval() +
  facet_wrap(~par, scales="free_x") +
  theme_bw(base_size = 14) +
  labs(x = "Estimate", y = NULL)

## needs a little tidying up to get better names and perhaps line up the 
## two duration parameters
## demonstrates a couple of things. 1. which effects have evidence of being
## non-zero (w, W:conspecific for growth, height); 
## (2) the correspondence between the effects on survival and growth. Note 
## how W:control is negative in both cases.
pivot_wider(filter(draws, !(variable  %in% c("Intercept"))), 
            names_from=par, values_from = value)|> 
  group_by(variable) |> median_hdci(biomass, survival, .width = c(.95)) |>
  ggplot(aes(y = biomass, ymin = biomass.lower, ymax=biomass.upper,
             x = survival, xmin = survival.lower, xmax = survival.upper)) +
  geom_pointinterval() +
  theme_bw(base_size = 14)



## now the variances
draws_var <- tidy_draws(brms_bm.mod) |> 
  dplyr::select(starts_with("."), starts_with("sd_"))
draws_var <- pivot_longer(draws_var, starts_with("sd_"), names_to = "par") |> 
  mutate(par = str_sub(par, start = 4)) |> 
  mutate(par = str_remove(par, "_Intercept")) |>
  mutate(par = str_remove(par, "_$")) |> 
  separate_wider_delim(par, delim = "_", names=c("par", "vital_rate"), 
                       too_few="align_start", too_many = "merge") |> 
  mutate(vital_rate = replace_na(vital_rate, "biomass"),
         vital_rate = case_when(
           vital_rate %in% c("Tag", "biomass") ~ "biomass",
           vital_rate %in% c("Tag__hu", "_hu") ~ "survival"))

table(draws_var$vital_rate, draws_var$par) ## looks ok

tree_plot_var <- draws_var |> group_by(par, vital_rate) |> 
  median_hdci(value, .width = c(0.95)) |>
  ggplot(aes(y = par, x = value, xmin = .lower, xmax = .upper)) +
  facet_wrap(~vital_rate, scales="free_x") +
  geom_pointinterval() +
  theme_bw(base_size = 14) +
  labs(x = "Standard deviation", y = NULL)

tree_plot_pop / tree_plot_var +
  plot_layout(heights = c(0.7, 0.3)) # looks okay, but might want to 
## switch out the names to something more understandable to a reader.


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
                      prior=prior_bm_mod)

pp_check(brms_bm.mod_sp, ndraws=100) +
  pp_check(brms_bm.mod_sp, resp = hu, ndraws=100) ## not bad

summary(brms_bm.mod_sp)

pred_dat_sp <- with(dat2, 
                 expand_grid(Treatment = unique(Treatment),
                             soil = levels(soil), soil.dim1 = 0,
                             Cohort = levels(Cohort)[1],
                             Seedling = unique(dat2$Seedling),
                             duration = 212, Init_Height = 7))


pred_dat_sp <- add_epred_draws(brms_bm.mod_sp, newdata = pred_dat_sp, 
                             dpar = TRUE,
                             re_formula = ~ (1 + soil*Treatment||Seedling) )

plot_dat_sp <- pred_dat_sp |> group_by(soil, Treatment, Seedling) |> 
  median_hdci(mu, hu, .epred)

## combined biomass and survival
ggplot(plot_dat_sp, aes(x = soil, colour = Treatment, 
             y=.epred, ymin = .epred.lower, ymax = .epred.upper)) + 
  facet_wrap(~Seedling, scales = "free_y") + 
  geom_pointrange(position = position_dodge2(width = 0.3)) +
  theme_bw()


plot_dat_sp <- 
  ungroup(pred_dat_sp) |> select(Treatment, soil, Seedling, `.draw`:hu)|> 
  mutate(.epred = log(.epred), hu = log(1-hu)) |> 
  rename("combined" = ".epred", "survival" = "hu", "biomass" = "mu") |> 
  pivot_longer(cols = combined:survival, names_to = "vital_rate") |> 
  pivot_wider(names_from=soil, values_from=value) |> 
  pivot_longer(cols=conspecific:control, names_to="soil") |> 
  mutate(psf = value - heterospecific,
         Treatment = factor(Treatment, labels = c("50% WHC", "100% WHC")),
         vital_rate = factor(vital_rate, 
                             levels = c("combined", "biomass", "survival")),
         soil = factor(ifelse(soil == "conspecific", "live", "sterile")))|> 
  group_by(Treatment, Seedling, soil, vital_rate) |>
  median_hdci(psf, .width = c(0.66, 0.95))

ggplot(filter(plot_dat_sp, vital_rate == "combined"),
       aes(x = soil, y=psf, colour = Treatment,
           ymin = .lower, ymax = .upper)) +
  facet_wrap(~Seedling, scales = "free_y") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "combined",
                                .width == 0.66), linewidth = 1,
                  position = position_dodge2(width=0.2) ) +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "combined", 
                                .width == 0.95), 
                  position = position_dodge2(width=0.2)) +
  theme_bw()



ggplot(filter(plot_dat_sp, vital_rate == "survival"),
       aes(x = soil, y=psf, colour = Treatment,
           ymin = .lower, ymax = .upper)) +
  facet_wrap(~Seedling, scales = "free_y") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "survival",
                                .width == 0.66), linewidth = 1,
                  position = position_dodge2(width=0.2) ) +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "survival", 
                                .width == 0.95), 
                  position = position_dodge2(width=0.2)) +
   theme_bw()

  
ggplot(filter(plot_dat_sp, vital_rate == "biomass"),
       aes(x = soil, y=psf, colour = Treatment,
           ymin = .lower, ymax = .upper)) +
  facet_wrap(~Seedling, scales = "free_y") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "biomass",
                                .width == 0.66), linewidth = 1,
                  position = position_dodge2(width=0.2) ) +
  geom_pointrange(data = filter(plot_dat_sp, vital_rate == "biomass", 
                                .width == 0.95), 
                  position = position_dodge2(width=0.2)) +
  theme_bw()

