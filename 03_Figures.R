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

dat2 <- read.csv("Output/Models/dat2.csv")
load("Output/Models/brms_bm.mod") # hurdle model
load("Output/Models/brms_RMF.mod") # biomass allocation model
load("Output/Models/brms_bm.mod_sp") # species-specific model

### Figure 2 ------------

# Forest Plot
## present all model coefficients
draws <- tidy_draws(brms_bm.mod) |> 
  dplyr::select(starts_with("."), starts_with("b_"))
draws <- pivot_longer(draws, starts_with("b_"), names_to = "par") |> 
  mutate(par = str_sub(par, start = 3))

#pop_ranef_est <- read.csv(file = "Data/bm_pop_ranef.csv")
## identify mortality and growth parameters.
draws <- draws |> separate_wider_delim(par, delim="u_",
                                       names=c("par", "variable"),
                                       too_few = "align_end"  ) |>
  mutate(par = factor(ifelse(is.na(par), "Biomass", "Survival"),
                      levels = c("Survival", "Biomass")),
         value = value * ifelse(par == "Survival", -1, 1),
         variable = case_when(variable == "TreatmentW:soilcontrol" ~
                                "Sterile Soil:\n 100% WHC",
                              variable == "TreatmentW:soilconspecific" ~
                                "Conspecific Soil:\n 100% WHC",
                              variable == "TreatmentW" ~ "100% WHC",
                              variable == "soilcontrol" ~ "Sterile Soil",
                              variable == "soilconspecific" ~
                                "Conspecific Soil",
                              variable == "soil.dim1" ~ "Geochemistry",
                              variable == "scaleInit_Height" ~ "Initial Height",
                              variable == "scaleduration" ~ "Duration",
                              variable == "Cohort2" ~ "Cohort",
                              variable == "Intercept" ~ "Intercept")) ## convert mortality to survival

draws$variable <- factor(draws$variable,
                         levels = c("Cohort",
                                    "Duration",
                                    "Initial Height",
                                    "Geochemistry",
                                    "Conspecific Soil:\n 100% WHC",
                                    "Sterile Soil:\n 100% WHC",
                                    "100% WHC",
                                    "Conspecific Soil",
                                    "Sterile Soil",
                                    "Intercept"))

Figure2A <- filter(draws, !(variable  %in% c("Intercept"))) |> 
  group_by(par, variable) |> median_hdci(.width = c(.95)) |>
  ggplot(aes(y = variable, x = value, xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept=0, linetype = "dashed") +
  geom_pointinterval() +
  facet_wrap(~par, scales="free_x") +
  theme_bw(base_size = 14) +
  labs(x = "Estimate", y = NULL)

Figure2A

## now the variances
draws_var <- tidy_draws(brms_bm.mod) |> 
  dplyr::select(starts_with("."), starts_with("sd_"))
draws_var <- pivot_longer(draws_var, starts_with("sd_"), names_to = "par") |> 
  mutate(par = str_sub(par, start = 4)) |> 
  mutate(par = str_remove(par, "_Intercept")) |>
  mutate(par = str_remove(par, "_$")) |> 
  separate_wider_delim(par, delim = "_", names=c("par", "vital_rate"), 
                       too_few="align_start", too_many = "merge") |> 
  mutate(vital_rate = replace_na(vital_rate, "Biomass"),
         vital_rate = case_when(
           vital_rate %in% c("Tag", "Biomass") ~ "Biomass",
           vital_rate %in% c("Tag__hu", "_hu") ~ "Survival"),
         par = case_when(par == "Table" ~ "Greenhouse \nLocation",
                         par == "Seedling" ~ "Species",
                         par == "Inoc" ~ "Inoc. Source",
                         par == "begin" ~ "Experiment \nStart Date"))

#table(draws_var$vital_rate, draws_var$par) ## looks ok
draws_var$vital_rate <- factor(draws_var$vital_rate,
                               levels = c("Survival", "Biomass"))

## Combine the figures
Figure2B <- draws_var |> group_by(par, vital_rate) |> 
  median_hdci(value, .width = c(0.95)) |>
  ggplot(aes(y = par, x = value, xmin = .lower, xmax = .upper)) +
  facet_wrap(~vital_rate, scales="free_x") +
  geom_pointinterval() +
  theme_bw(base_size = 14) +
  labs(x = "Standard Deviation", y = NULL)

Figure2B

Figure2 <- Figure2A / Figure2B +
  plot_layout(heights = c(0.7, 0.3)) + plot_annotation(tag_levels = "A") 

Figure2

ggsave(plot = Figure2, filename = "Output/MSFigs/Fig2.png")

### Figure 3 -----------
## Plot the hurdle model output
#create prediction matrix
pred_dat_pop <- with(dat2,
                     expand_grid(Treatment = unique(Treatment),
                                 soil = levels(soil), 
                                 soil.dim1 = 0, duration = 212, Init_Height = 7,
                                 Cohort = levels(Cohort)[1]))
## get posterior draws of biomass (mu) and mortality (hu)
pred_dat_pop <- add_epred_draws(brms_bm.mod, newdata = pred_dat_pop, 
                                dpar = c("mu", "hu"),
                                re_formula = ~ 0)

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

plot_psf_pop2 <- plot_psf_pop |> filter(soil == "live")

ggplot(plot_psf_pop2,
       aes(x = Treatment, y=psf, 
           ymin = .lower, ymax = .upper)) +
  facet_wrap(~vital_rate, scales = "free_y") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_pointrange(data = filter(plot_psf_pop2, .width == 0.95),
                  position = position_dodge2(width = 0.2)) +
  geom_pointrange(data = filter(plot_psf_pop2, .width == 0.66),
                  position = position_dodge2(width = 0.2), linewidth = 1.) +
  labs(x = "Soil Moisture Treatment", y = "Log Response Ratio") +
  theme_bw()

## arrange more usefully
## log response ratio results
pl_psf_pop <- group_by(plot_psf_pop2, vital_rate) |> 
  group_map(\(.x, ...) 
            ggplot(.x,
                   aes(x = Treatment, y=psf,
                       ymin = .lower, ymax = .upper)) +
              geom_hline(yintercept=0, linetype = "dashed") +
              geom_pointrange(data = filter(.x, .width == 0.95),
                              position = position_dodge2(width = 0.2)) +
              geom_pointrange(data = filter(.x, .width == 0.66),linewidth = 1,
                              position = position_dodge2(width = 0.2)) +
              #labs(x = NULL, y = NULL) +
              labs(x = "Soil Moisture Treatment", y = "Log Response Ratio") +
              theme_bw())

Figure3 <- (((pl_psf_pop[[1]] + labs(tag = "A", x= NULL)) | 
               ((pl_psf_pop[[2]] + labs(tag = "B", y = NULL, x = NULL))/
                  (pl_psf_pop[[3]] + labs(tag = "C", y = NULL, x = NULL)))) +
              plot_layout(guides = "collect"))  +
  plot_annotation(caption = "Soil Moisture Treatment") &
  theme(plot.caption = element_text(hjust = 0.5, size = 12))

Figure3

ggsave(plot = Figure3, filename = "Output/MSFigs/Fig3.png")
### Figure 4 -------
pred_dat_RMF <- with(filter(dat2, status == 1),
                     expand_grid(Treatment = unique(Treatment),
                                 soil = levels(soil), 
                                 soil.dim1 = 0, duration = 212, Init_Height = 7,
                                 Cohort = levels(Cohort)[1]))
## get posterior draws of RMF (mu)
pred_dat_RMF <- add_epred_draws(brms_RMF.mod, newdata = pred_dat_RMF, 
                                dpar = c("mu"),
                                re_formula = ~ 0)

## quick plot
pred_dat_RMF |> ggplot(aes(x = soil, y = .epred)) +
  facet_wrap(~Treatment) + stat_pointinterval(point_interval = "median_hdci")

## calculate credible intervals
plot_dat_RMF <- pred_dat_RMF |> group_by(soil, Treatment) |> 
  median_hdci(mu, .epred, .width=c(0.95)) |> 
  rename("est" = "mu", "epred.est" = ".epred", 
         "epred.lower" = ".epred.lower", "epred.upper" = ".epred.upper") |> 
  pivot_longer(cols=est:epred.upper) |> 
  #separate_wider_delim(name, delim =".", names=c("stat")) |> 
  mutate(value =  exp(value)) |> 
  pivot_wider(names_from=name, values_from = value) |>
  mutate(soil = case_when(soil == "conspecific" ~ "Conspecific",
                          soil == "heterospecific" ~ "Heterospecific",
                          soil == "control" ~ "Sterile"))

plot_dat_RMF$soil <- factor(plot_dat_RMF$soil,
                            levels = c("Sterile",
                                       "Conspecific",
                                       "Heterospecific"))

Figure4 <- ggplot(filter(plot_dat_RMF, soil != "Sterile"),
                  aes(x = soil, y=est)) +
  geom_point(stat = "identity", position = position_dodge(0.9)) +
  geom_pointrange(mapping = aes( ymin = mu.lower, ymax = mu.upper),
                  position = position_dodge(0.9)) +
  facet_wrap(~Treatment,
             labeller = as_labeller(c("D" = "50% WHC", "W" = "100% WHC"))) +
  labs(x = "Soil Inoculum Source", y = "Root Mass Fraction") +
  theme_bw(12)
Figure4

ggsave(plot = Figure4, filename = "Output/MSFigs/Fig4.png")

### Figure 5 -------------
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


plot_dat_sp2 <- plot_dat_sp |> filter(soil != "sterile")

Figure5 <- ggplot(filter(plot_dat_sp2, vital_rate == "combined"),
                  aes(x = Treatment, y=psf,
                      ymin = .lower, ymax = .upper)) +
  facet_wrap(~Seedling, scales = "free_y",
             labeller = as_labeller(c("COCCMA" ="Coccoloba manzanillensis", 
                                      "EUGENE" =  "Eugenia nesiotica",
                                      "GUAPST" ="Guapira standylana",
                                      "HEISCO" = "Heisteria concinna",
                                      "LACIAG" = "Lacistema aggregatum",
                                      "PSIDFR" =  "Psidium friedrichsthalianum"))) +
  
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_pointrange(data = filter(plot_dat_sp2, vital_rate == "combined",
                                .width == 0.66), linewidth = 1,
                  position = position_dodge2(width=0.2) ) +
  geom_pointrange(data = filter(plot_dat_sp2, vital_rate == "combined", 
                                .width == 0.95), 
                  position = position_dodge2(width=0.2)) +
  labs(y = "Log Response Ratio",
       x = "Soil Moisture Treatment") +
  
  theme_bw() +
  theme(strip.text = element_text(face = "italic")) 
Figure5

ggsave(plot= Figure5, filename = "Output/MSFigs/Fig5.png")
