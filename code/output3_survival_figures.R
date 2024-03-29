

# packages, setup ----
library(tidyverse)
library(ggplot2)
library(RMark)
library(here)
library(cowplot)
options(scipen = 999)


# read in raw data to get observed ranges of predictors 
ana_table <- readRDS(here('data/ana_table_full')) 

ana_table %>% 
  filter(!is.na(sex), !is.na(cr.age)) %>% 
  count(SEASON, sex, OUTCOME) 


# ana_table %>% filter(!is.na(sex), !is.na(cr.age)) %>% group_by(SEASON, sex) %>% summarise_at("cr.mass", list(min = min, mean = mean, max = max, num.chx = length))


# sex_seas_phi_pims <- readRDS(here("fitted_models/sex_seas_phi_pims"))

# total PIM per group should be
# sum(48 - seq(0, 48))
# = 1176
# so total overall PIM = 
# 1176 * 4
# 4704

# define functions for calculating estimates and se, and data management ----
# add sex and season, based on parameter indices. only works for the first cohort, since only using the first cohort below
add_sex_seas <- function(preds) {
preds <- preds %>% 
  left_join(., readRDS(here("data/seas_sex_par_index")))
}



# add season day
add_day <- function(preds) {
preds <- preds %>% 
  mutate(day = rep(seq(50, 97), length.out = nrow(preds)))
}

# extraft daily survival estimates from list into df; for model predictions at different values of a predictor ('covdata')
get_daily_pred <- function(preds) {
  se <- preds$estimates %>% 
    data.frame() %>% 
    add_sex_seas() %>% 
    add_day()
}


# overall survival estimate for model predictions at different values of a predictor ('covdata')
get_overall_pred <- function(preds) {
  se <- preds$estimates %>% 
    data.frame() %>% 
    add_sex_seas()
if("covdata" %in% colnames(se)) {  
se <- se 
} else {
  se <- se %>% 
    rename(covdata = contains("cr."))
}
  se <- se %>%
    group_by(SEASON, sex, covdata) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
    ungroup()
}

#pred <- map_df(flip_preds, get_overall_pred) %>% ungroup()



# standard error for overall survival via delta method
get_overall_se <- function(preds) {
f12_se <- deltamethod.special("prod", preds$estimates$estimate[1 + 17:48], preds$vcv[1 + 17:48, 1 + 17:48])
f13_se <- deltamethod.special("prod", preds$estimates$estimate[49 + 17:96], preds$vcv[49 + 17:96, 49 + 17:96])
m12_se <- deltamethod.special("prod", preds$estimates$estimate[97 + 17:144], preds$vcv[97 + 17:144, 97 + 17:144])
m13_se <- deltamethod.special("prod", preds$estimates$estimate[145 + 17:192], preds$vcv[145 + 17:192, 145 + 17:192])

if("covdata" %in% colnames(preds$estimates)) {
  covdata = preds$estimates %>%
  data_frame() %>% 
    distinct(covdata)
} else {
covdata = preds$estimates %>%
  data_frame() %>% 
  select(contains("cr.")) %>% 
  distinct()
}

pred_se <- rbind(f12_se, f13_se, m12_se, m13_se) %>% 
  data_frame() %>% 
  rename(se = 1) %>% 
  mutate(sex = c("Female", "Female", "Male", "Male"),
         SEASON = c("1213", "1314", "1213", "1314"),
         covdata = covdata[[1]])

}


# standard error for overall survival via delta method
get_overall_se_creched <- function(preds) {
se12 <- deltamethod.special("prod", preds$estimates$estimate[1:38], preds$vcv[1:38, 1:38])
se13 <- deltamethod.special("prod", preds$estimates$estimate[39:76], preds$vcv[39:76, 39:76])

if("covdata" %in% colnames(preds$estimates)) {
  covdata = preds$estimates %>%
  data_frame() %>% 
    distinct(covdata)
} else {
covdata = preds$estimates %>%
  data_frame() %>% 
  select(contains("cr.")) %>% 
  distinct()
}

pred_se <- rbind(se12, se13) %>% 
  data_frame() %>% 
  rename(se = 1) %>% 
  mutate(SEASON = c("2012-13", "2013-14"),
         covdata = covdata[[1]])

}

# se <- map_df(flip_preds, get_overall_se) %>% ungroup()


# pred_se <- full_join(pred, se) %>% 
#  mutate(lci = overall.surv - (1.96 * se),
#         uci = overall.surv + (1.96 * se)) %>% 
#  rename(flipper.length = covdata)


#
# basic summaries
ana_table %>% 
  group_by(SEASON) %>% 
  summarise(min.cr = min(cr.age, na.rm = TRUE), 
            max.cr = max(cr.age, na.rm = TRUE), 
            mean.cr = mean(cr.age, na.rm = TRUE),
            num.chx = n(),
            se.cr = sd(cr.age, na.rm = TRUE) / sqrt(num.chx))

# Figure 1: base survival ----
# read model averaged estimates from CH4_base_survival.R
step2_mod_av_phi <- readRDS(here("fitted_models/survival/step2_mod_av_phi"))

# no cohort effects, all estimates equal for each cohort, so just take estimates for first cohort
step2_mod_av_est <- filter(step2_mod_av_phi$estimates, occ.cohort == 1)


# plot 
mod_av_daily <- step2_mod_av_est %>% 
  mutate(time = as.numeric(as.character(time)),
         SEASON = ifelse(SEASON == 1213, "2012-13", "2013-14")) %>% 
  ggplot() +
  #geom_point(aes(x = time, y = estimate)) +
	geom_line(aes(x = time, y = estimate, group = group)) +
	geom_line(aes(x = time,  y = ucl, group = group), linetype="dotted") +
	geom_line(aes(x = time, y = lcl, group = group), linetype="dotted") +
  facet_grid(sex~SEASON) +
  labs(x = "Day of season",
       y = "Estimated daily survival") +
  theme_bw()

ggsave(here("figures/overall_survival.png"), width = 6, height = 6, dpi = 300)


# adding plot for overall model averaged survival
f12_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[1:48], step2_mod_av_phi$vcv.real[1:48, 1:48])
f13_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[49:96], step2_mod_av_phi$vcv[49:96, 49:96])
m12_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[97:144], step2_mod_av_phi$vcv[97:144, 97:144])
m13_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[145:192], step2_mod_av_phi$vcv[145:192, 145:192])

se <- rbind(f12_se, f13_se, m12_se, m13_se) %>% 
  data.frame() %>% 
  rename(se = 1) %>% 
  rownames_to_column("seas.sex") %>% 
  mutate(seas.sex = case_when(seas.sex == "f12_se" ~ "1213Female",
                              seas.sex == "f13_se" ~ "1314Female",
                              seas.sex == "m12_se" ~ "1213Male",
                              seas.sex == "m13_se" ~ "1314Male"))


# plot 
mod_av_overall <- step2_mod_av_est %>% 
  mutate(time = as.numeric(as.character(time)),
         SEASON = ifelse(SEASON == 1213, "2012-13", "2013-14")) %>% 
  rename(seas.sex = group) %>% 
  group_by(seas.sex) %>% 
  summarise(overall.est = prod(estimate)) %>% 
  full_join(se) %>% 
  mutate(lci = overall.est - (1.96 * se),
         uci = overall.est + (1.96 * se),
         lci = ifelse(lci < 0, 0, lci)) %>% 
  mutate(seas.sex = case_when(seas.sex == "1213Female" ~ "2012\nFemale",
                              seas.sex == "1314Female" ~ "2013\nFemale",
                              seas.sex == "1213Male" ~ "2012\nMale",
                              seas.sex == "1314Male" ~ "2013\nMale")) %>% 
  ggplot() +
	geom_point(aes(x = seas.sex, y = overall.est)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, x = seas.sex, width = 0.25)) +
  labs(x = "",
       y = "Estimated overall survival") +
  theme_bw() 

plot_grid(mod_av_daily, 
          mod_av_overall, rel_widths = c(0.7, 0.3),
          align = "h")

plot_grid(mod_av_daily, 
          mod_av_overall,
          ncol = 1)
#
# Figure 1 alt: base survival from best model ----

step2_best <- readRDS(here("fitted_models/survival/step2"))$phi.year_T_hatch_incr_p.sat

step2_3rd <- readRDS(here("fitted_models/survival/step2"))$phi.year_T.incr_hatch_p.sat

df <- data.frame(matrix(rep(c(0, 1), 49), nrow = 2))

names(df) <- paste("in.cr", 50:98, sep = "")


pred_cr_effect <- function(zdf.ind) {
zday_pred <- covariate.predictions(step2_3rd, data=df[zdf.ind], indices = c(zdf.ind, zdf.ind + 1176), alpha=0.025)$estimates
}

zday_preds = map_df(seq(1:48), pred_cr_effect)


# plot 
best_daily <- zday_preds %>% 
  mutate(SEASON = ifelse(par.index < 1000, "2012-13", "2013-14"), 
         guard.cr = ifelse(covdata == 0, "Guarded", paste("Cr", "\u00E8", "ched", sep = "")),
         time = ifelse(SEASON == "2012-13", model.index + 49, model.index + 1)) %>% 
  filter((SEASON == "2012-13" & covdata == 0 & time <= 74) |
           (SEASON == "2012-13" & covdata == 1 & time >= 60) |
           (SEASON == "2013-14" & covdata == 0 & time <= 77) |
           (SEASON == "2013-14" & covdata == 1 & time >= 63)) %>% 
  ggplot() +
  #geom_point(aes(x = time, y = estimate)) +
	geom_line(aes(x = time, y = estimate, color = guard.cr)) +
  geom_ribbon(aes(x = time, ymin = lcl, ymax = ucl, fill = guard.cr), alpha = 0.2) +
#	geom_line(aes(x = time,  y = ucl, color = guard.cr), linetype="dotted") +
#	geom_line(aes(x = time, y = lcl, color = guard.cr), linetype="dotted") +
  facet_grid(~SEASON) +
  labs(x = "Day of season",
       y = "Estimated daily survival",
       title = "Model Phi(year + hatch date + T * creche)") +
  theme_bw()+ 
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.title = element_blank())

best_daily

ggsave(here("figures/overall_survival_T_X_creche.png"), width = 6, height = 6, dpi = 300)


# adding plot for overall model averaged survival
f12_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[1:48], step2_mod_av_phi$vcv.real[1:48, 1:48])
f13_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[49:96], step2_mod_av_phi$vcv[49:96, 49:96])
m12_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[97:144], step2_mod_av_phi$vcv[97:144, 97:144])
m13_se <- deltamethod.special("prod", step2_mod_av_phi$estimates$estimate[145:192], step2_mod_av_phi$vcv[145:192, 145:192])

se <- rbind(f12_se, f13_se, m12_se, m13_se) %>% 
  data.frame() %>% 
  rename(se = 1) %>% 
  rownames_to_column("seas.sex") %>% 
  mutate(seas.sex = case_when(seas.sex == "f12_se" ~ "1213Female",
                              seas.sex == "f13_se" ~ "1314Female",
                              seas.sex == "m12_se" ~ "1213Male",
                              seas.sex == "m13_se" ~ "1314Male"))


# plot 
mod_av_overall <- step2_mod_av_est %>% 
  mutate(time = as.numeric(as.character(time)),
         SEASON = ifelse(SEASON == 1213, "2012-13", "2013-14")) %>% 
  rename(seas.sex = group) %>% 
  group_by(seas.sex) %>% 
  summarise(overall.est = prod(estimate)) %>% 
  full_join(se) %>% 
  mutate(lci = overall.est - (1.96 * se),
         uci = overall.est + (1.96 * se),
         lci = ifelse(lci < 0, 0, lci)) %>% 
  mutate(seas.sex = case_when(seas.sex == "1213Female" ~ "2012\nFemale",
                              seas.sex == "1314Female" ~ "2013\nFemale",
                              seas.sex == "1213Male" ~ "2012\nMale",
                              seas.sex == "1314Male" ~ "2013\nMale")) %>% 
  ggplot() +
	geom_point(aes(x = seas.sex, y = overall.est)) +
    geom_errorbar(aes(ymin = lci, ymax = uci, x = seas.sex, width = 0.25)) +
  labs(x = "",
       y = "Estimated overall survival") +
  theme_bw() 

plot_grid(mod_av_daily, 
          mod_av_overall, rel_widths = c(0.7, 0.3),
          align = "h")

plot_grid(mod_av_daily, 
          mod_av_overall,
          ncol = 1)
#
# end Figure 1 ----
# Figure 2 overall survival probability plot - flipper model ----

# creching size survival ----

cr_age_size_models <- readRDS(here("fitted_models/survival/step6_cr_timing_surv"))


#model.table(cr_age_size_models, use.lnl = TRUE)

#
# survival estimates based on flipper model ----
#summary(ana_table$cr.flip)
cr_flip_span = seq(55, 155, length.out = 10)

flip_mod <- cr_age_size_models$phi.year_T_crflip_p.sat

PIMS(flip_mod, parameter = "Phi", simplified = FALSE)


make_flip_pred <- function(zcr.flip) {
flip.pred=covariate.predictions(flip_mod, data=data.frame(cr.flip = zcr.flip), indices=c(1:38, 742:779), alpha=0.025)
}

flip_preds <- map(cr_flip_span, make_flip_pred)
  

flip_preds_daily <- map_df(flip_preds, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  rename(cr.flip = covdata) %>% 
  mutate(SEASON = ifelse(between(par.index, 1, 38), "2012-13", "2013-14"),
         day = rep(c(60:97), length.out = 760))

#
# NO RUN daily survival probability plot - flipper model ----

flip_preds_daily %>% 
#  filter(between(par.index, 1, 38) | between(par.index, 1177, 1177 + 47) | between(par.index, 2353, 2353 + 47) | between(par.index, 3529, 3529 + 47)) %>% 
  ggplot(group = cr.flip)+
	geom_line(aes(x = day, y = estimate, color = as.factor(cr.flip)), size=1) +
	#geom_line(aes(x = day, y = ucl, color = flipper.length), linetype = "dotted", size = 1) +
	#geom_line(aes(x = day, y = lcl, color = flipper.length), linetype = "dotted", size = 1) +
	scale_colour_discrete(name="Creching\nflipper (mm)")+
	#geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(~SEASON)+
	labs(y = "Daily survival probability",
	     x = "Day of season",
	     title = "Creched chicks")+
	ylim(0.65,1)+
	theme_bw()+
				theme(plot.background = element_blank(),
				      panel.grid.minor = element_blank(),
				      panel.grid.major = element_blank(),
				      strip.background = element_blank(),
				      axis.ticks = element_blank(),
				      axis.text = element_text(size = 15),
				      axis.title = element_text(size = 20),
				      strip.text = element_text(size = 15)
					)  + 
      guides(color = guide_legend(reverse = TRUE))

ggsave(here("../MSdocs/daily_survival_by_flipper_cr_chicks.png"), width = 6) 
  
# overall survival probability plot - flipper model ----

flip_overall_pred <- flip_preds_daily %>% 
    add_sex_seas() %>%
  filter(day > 68) %>% 
    group_by(SEASON, cr.flip) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON))
  

flip_overall_se <- map_df(flip_preds, get_overall_se_creched) %>% 
  ungroup() %>% 
  rename(cr.flip = covdata)


flip_overall_pred_se <- full_join(flip_overall_pred, flip_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))

# change title to reflect which results
flip_surv_plot <- flip_overall_pred_se %>% 
    mutate(lci = ifelse(lci < 0, 0, lci),
           uci = ifelse(uci > 1, 1, uci)) %>% 
    ggplot() +
    geom_line(aes(cr.flip, overall.surv)) +
    geom_ribbon(aes(x = cr.flip, ymin = lci, ymax = uci), alpha = 0.2) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
    labs(x = "Creching flipper length (mm)",
         y = "") +
    theme_bw() +
    facet_grid(~SEASON)
  
 #ggsave(here("../MSdocs/survival_by_flipper_cr_chicks.png"), width = 6) 
  
  
# make mass model plot ----
#  summary(ana_table$cr.mass)
  zcr.mass <- seq(300, 2900, length.out = 10)

  mass_mod <- cr_age_size_models$phi.year_T_crmass_p.sat

make_mass_pred <- function(zcr.mass) {
  mass.pred=covariate.predictions(mass_mod, data=data.frame(cr.mass=zcr.mass), indices=c(1:38, 742:779), alpha=0.025)
}

mass_preds <- map(zcr.mass, make_mass_pred)

mass_preds_daily <-  map_df(mass_preds, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  rename(cr.mass = covdata) %>% 
  mutate(SEASON = ifelse(between(par.index, 1, 38), "2012-13", "2013-14"),
         day = rep(c(60:97), length.out = 760))

#
# NO RUN daily survival probability plot - mass model ----

mass_preds_daily %>% 
ggplot()+
	geom_line(aes(x = day, y = estimate, color = as.factor(cr.mass)), size=1) +
	#geom_line(aes(x = day, y = ucl, color = pred), linetype = "dotted", size = 1) +
	#geom_line(aes(x = day, y = lcl, color = pred), linetype = "dotted", size = 1) +
	scale_colour_discrete(name="Creching\nmass (g)")+
	#geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(~SEASON)+
	labs(y = "Daily survival probability",
	     x = "Days past hatching",
	     title = "Creched chicks")+
	ylim(0.65,1)+
	theme_bw()+
				theme(plot.background = element_blank(),
				      panel.grid.minor = element_blank(),
				      panel.grid.major = element_blank(),
				      strip.background = element_blank(),
				      axis.ticks = element_blank(),
				      axis.text = element_text(size = 15),
				      axis.title = element_text(size = 20),
				      strip.text = element_text(size = 15)
					)  + 
      guides(color = guide_legend(reverse = TRUE))

ggsave(here("../MSdocs/daily_survival_by_mass_cr_chicks.png")) 

# overall survival probability plot - mass model ----

mass_overall_pred <- mass_preds_daily %>% 
    add_sex_seas() %>%
  filter(day > 68) %>% 
    group_by(SEASON, sex, cr.mass) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON))
  

mass_overall_se <- map_df(mass_preds, get_overall_se_creched) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata)


mass_overall_pred_se <- full_join(mass_overall_pred, mass_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))

# save to have for discussion
mini_mass_preds <- map(c(1350, 1652), make_mass_pred)

mini_mass_preds_daily <-  map_df(mini_mass_preds, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  rename(cr.mass = covdata) %>% 
  mutate(SEASON = ifelse(between(par.index, 1, 38), "2012-13", "2013-14"),
         day = rep(c(60:97), length.out = 152))

mini_mass_overall_pred <- mini_mass_preds_daily %>% 
    add_sex_seas() %>%
  filter(day > 68) %>% 
    group_by(SEASON, sex, cr.mass) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON))
  

mini_mass_overall_se <- map_df(mini_mass_preds, get_overall_se_creched) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata)


mini_mass_overall_pred_se <- full_join(mini_mass_overall_pred, mini_mass_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))

# ----
saveRDS(mass_overall_pred_se, here("fitted_models/survival/mass_overall_pred_se"))

mass_surv_plot <- mass_overall_pred_se %>% 
    mutate(lci = ifelse(lci < 0, 0, lci),
           uci = ifelse(uci > 1, 1, uci),
           cr.mass = cr.mass/1000) %>% 
    ggplot() +
    geom_line(aes(cr.mass, overall.surv)) +
    geom_ribbon(aes(x = cr.mass, ymin = lci, ymax = uci), alpha = 0.2) +
    scale_x_continuous(breaks = seq(0.5, 3, by = 0.5), labels = seq(0.5, 3, by = 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
    labs(x = "Creching mass (kg)",
         y = "") +
    theme_bw() +
    facet_grid(~SEASON)



 # make tib model plot ----
#  summary(ana_table$cr.tib)
  zcr.tib <- seq(75, 140, length.out = 10)

  tib_mod <- cr_age_size_models$phi.year_T_crtib_p.sat

make_tib_pred <- function(zcr.tib) {
  tib.pred=covariate.predictions(tib_mod, data=data.frame(cr.tib=zcr.tib), indices=c(1:38, 742:779), alpha=0.025)
}

tib_preds <- map(zcr.tib, make_tib_pred)

tib_preds_daily <- map_df(tib_preds, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  rename(cr.tib = covdata) %>% 
  mutate(SEASON = ifelse(between(par.index, 1, 38), "2012-13", "2013-14"),
         day = rep(c(60:97), length.out = 760))

#
# daily survival probability plot - tib model ----

tib_preds_daily %>% 
ggplot()+
	geom_line(aes(x = day, y = estimate, color = as.factor(cr.tib)), size=1) +
	#geom_line(aes(x = day, y = ucl, color = pred), linetype = "dotted", size = 1) +
	#geom_line(aes(x = day, y = lcl, color = pred), linetype = "dotted", size = 1) +
	scale_colour_discrete(name="Creching\ntib (mm)")+
	#geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(~SEASON)+
	labs(y = "Daily survival probability",
	     x = "Days past hatching",
	     title = "Creched chicks")+
	ylim(0.65,1)+
	theme_bw()+
				theme(plot.background = element_blank(),
				      panel.grid.minor = element_blank(),
				      panel.grid.major = element_blank(),
				      strip.background = element_blank(),
				      axis.ticks = element_blank(),
				      axis.text = element_text(size = 15),
				      axis.title = element_text(size = 20),
				      strip.text = element_text(size = 15)
					)  + 
      guides(color = guide_legend(reverse = TRUE))

ggsave(here("../MSdocs/daily_survival_by_tib_cr_chicks.png")) 

# overall survival probability plot - tib model ----

tib_overall_pred <-  tib_preds_daily %>% 
    add_sex_seas() %>%
  filter(day > 68) %>% 
    group_by(SEASON, sex, cr.tib) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON),
         sex = as.character(sex))
  

tib_overall_se <- map_df(tib_preds, get_overall_se_creched) %>% 
  ungroup() %>% 
  rename(cr.tib = covdata)


tib_overall_pred_se <- full_join(tib_overall_pred, tib_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))


tib_surv_plot <- tib_overall_pred_se %>% 
    mutate(lci = ifelse(lci < 0, 0, lci),
           uci = ifelse(uci > 1, 1, uci)) %>% 
    ggplot() +
    geom_line(aes(cr.tib, overall.surv)) +
    geom_ribbon(aes(x = cr.tib, ymin = lci, ymax = uci), alpha = 0.2) +
    scale_x_continuous(breaks = seq(80, 140, by = 20), labels = seq(80, 140, by = 20)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
    labs(x = "Creching tibiotarsus length (mm)",
         y = "") +
    theme_bw() +
    facet_grid(~SEASON)
  
  
  
 # ggsave(here("../MSdocs/survival_by_tib_cr_chicks.png")) 
  



# combine all together ----
title <- ggdraw() + 
  draw_label(
    "Probability of surviving mean creche day to flegde",
    angle = 90,
    size = 12
  )

share.marg = c(0.1, 0.5, 0.4, -.1)

plot_grid(title, plot_grid(flip_surv_plot +
  theme(plot.margin=unit(share.marg,"cm")), 
  mass_surv_plot +
  theme(plot.margin=unit(share.marg,"cm")),
  tib_surv_plot +
  theme(plot.margin=unit(share.marg,"cm")), 
          ncol = 1), ncol = 2, rel_widths = c(0.05, 0.95))
 

ggsave(here("figures/cr_size_survival.png"), width = 4, height = 6)


# end Figure 2 ----
# Figure 3: creche age, growth, survival surface plot ----
# mass  
mass_mod <- readRDS(here("fitted_models/survival/step6_cr_timing_surv"))$phi.year_T_crmass_p.sat

dat <- filter(ana_table, weight.slope40 > 0, !is.na(weight.slope40), !is.na(cr.age), !is.na(sex))
cr_mass_mod <- readRDS(here("fitted_models/cr_timing/cr_mass_step2_mods"))$year_hatch_mass_age

mass_cr_age_newdat <- expand.grid(sex = c("F"),
                                  SEASON = factor(c(1213, 1314)),
                                  cr.age = seq(min(dat$cr.age), max(dat$cr.age), length.out = 10),
                                  weight.slope40 = seq(min(dat$weight.slope40), max(dat$weight.slope40), length.out = 10),
                                  resid.hatch = 0)

mean_cr_age_pred = predict(cr_mass_mod, mass_cr_age_newdat, interval = "confidence") %>% 
  data.frame() %>% 
  cbind(., mass_cr_age_newdat)

# functions to get survival predictions at each creching mass
make_mass_pred_12 <- function(zcr.mass) {
  mass.pred=covariate.predictions(mass_mod, data=data.frame(cr.mass=zcr.mass), indices=c(1:38), alpha=0.025)
}
make_mass_pred_13 <- function(zcr.mass) {
  mass.pred=covariate.predictions(mass_mod, data=data.frame(cr.mass=zcr.mass), indices=c(742:779), alpha=0.025)
}


make_mass_pred_12(1652) %>% 
  get_daily_pred()  %>% 
  rename(cr.mass = covdata)%>%
  filter(day > 68) %>% 
    group_by(SEASON, sex, cr.mass) %>% 
    summarise(overall.surv = prod(estimate))
# fit functions and extract predictions
# sticking all this in a wrapper function to keep env cleaner
make_massgr_crage_overall_surv <- function() {
zcr.mass_12 <- filter(mean_cr_age_pred, SEASON == 1213) %>% select(fit) 
mass_pred_12 <- map(zcr.mass_12$fit, make_mass_pred_12)
mass_pred_12_daily <- map_df(mass_pred_12, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata) %>% 
  mutate(SEASON = "1213",
         day = rep(c(60:97), length.out = 3800))

zcr.mass_13 <- filter(mean_cr_age_pred, SEASON == 1314) %>% select(fit) 
mass_pred_13 <- map(zcr.mass_13$fit, make_mass_pred_13)
mass_pred_13_daily <- map_df(mass_pred_13, function(zpred) {zpred$estimates %>% data.frame()}) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata) %>% 
  mutate(SEASON = "1314",
         day = rep(c(60:97), length.out = 3800))

massgr_crage_overall_surv <- rbind(mass_pred_12_daily, mass_pred_13_daily) %>%
  filter(day > 68) %>% 
    group_by(SEASON, cr.mass) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  full_join(mean_cr_age_pred %>% rename(cr.mass = fit))
}

massgr_crage_overall_surv <- make_massgr_crage_overall_surv()

mean_crage_massgr <- dat %>% 
  group_by(SEASON) %>% 
  summarise(mean.cr.age = round(mean(cr.age), 0),
            mean.mass.gr = round(mean(weight.slope40), 0)) %>% 
  mutate(SEASON = ifelse(SEASON == "1213", "2012-13", "2013-14"))

massgr_crage_overall_surv_plot <- massgr_crage_overall_surv %>% 
  mutate(SEASON = ifelse(SEASON == "1213", "2012-13", "2013-14")) %>% 
ggplot(aes(weight.slope40, cr.age, z = overall.surv)) +
  geom_contour_filled() +
  scale_y_continuous(breaks = seq(11, 26), labels = seq(11, 26)) +
  scale_x_continuous(breaks = seq(10, 130, by = 20), labels =  seq(10, 130, by = 20)) +
  labs(x = "Mass growth (g/d)",
       y = "Creching age") +
  scale_fill_brewer(palette = "Spectral", name = "Creching\nmass to\nachieve\nsurvival\nprobability =") +
  theme_bw() +
  geom_hline(data = mean_crage_massgr, aes(yintercept = mean.cr.age)) +
  geom_vline(data = mean_crage_massgr, aes(xintercept = mean.mass.gr)) +
  facet_wrap(~SEASON)


ggsave(here("figures/massgr_crage_survival.png"), massgr_crage_overall_surv_plot, width = 6, height = 4)




