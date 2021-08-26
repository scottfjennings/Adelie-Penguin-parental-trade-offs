

# packages, setup ----
library(tidyverse)
library(ggplot2)
library(RMark)
library(here)
options(scipen = 999)

# setting working directory only to control location for all the Mark output files
setwd(here("mark_output"))

# read in raw data to get observed ranges of predictors 
ana_table <- readRDS(here('data/ana_table_full')) 

ana_table %>% 
  filter(!is.na(sex), did.cr == 1) %>% 
  count(SEASON, sex, OUTCOME) 


ana_table %>% 
  filter(!is.na(sex), did.cr == 1) %>%
  group_by(SEASON, sex) %>% 
  summarise_at("cr.mass", list(min = min, mean = mean, max = max, num.chx = length))



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
f12_se <- deltamethod.special("prod", preds$estimates$estimate[1:48], preds$vcv[1:48, 1:48])
f13_se <- deltamethod.special("prod", preds$estimates$estimate[49:96], preds$vcv[49:96, 49:96])
m12_se <- deltamethod.special("prod", preds$estimates$estimate[97:144], preds$vcv[97:144, 97:144])
m13_se <- deltamethod.special("prod", preds$estimates$estimate[145:192], preds$vcv[145:192, 145:192])

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

# se <- map_df(flip_preds, get_overall_se) %>% ungroup()


# pred_se <- full_join(pred, se) %>% 
#  mutate(lci = overall.surv - (1.96 * se),
#         uci = overall.surv + (1.96 * se)) %>% 
#  rename(flipper.length = covdata)


#
# base survival ----
phi_step2 <- readRDS(here("fitted_models/phi_step1"))

##--model averaging for step 2 models- this does real estimates
step2_mod_av = model.average(phi_step2, "Phi", vcv=TRUE, drop = FALSE)



step2_mod_av <- readRDS(here("fitted_models/step2_mod_av"))

# no cohort effects, all estimates equal for each cohort, so just take estimates for first cohort
step2_mod_av_est <- filter(step2_mod_av$estimates, occ.cohort == 1)


# plot 
step2_mod_av_est %>% 
  mutate(time = as.numeric(as.character(time)),
         days.past.hatch = time - 49) %>% 
  ggplot() +
  #geom_point(aes(x = time, y = estimate)) +
	geom_line(aes(x = days.past.hatch, y = estimate, group = group)) +
	geom_line(aes(x = days.past.hatch,  y = ucl, group = group), linetype="dotted") +
	geom_line(aes(x = days.past.hatch, y = lcl, group = group), linetype="dotted") +
  facet_grid(sex~SEASON) +
  labs(x = "Days past hatching",
       y = "Estimated daily survival") +
  theme_bw()

ggsave(here("../MSdocs/fig1_overall_survival.png"), width = 6, height = 6, dpi = 300)


# creching size survival ----

# 2 options here; at final step, models fitted on all data and on just the chicks that reached creche
cr_age_size_models <- readRDS(here("fitted_models/cr_age_size_models"))
cr_age_size_models <- readRDS(here("fitted_models/cr_age_size_models_cronly"))


model.table(cr_age_size_models, use.lnl = TRUE)
# survival estimates based on flipper model ----
summary(ana_table$cr.flip)
cr_flip_span = seq(88, 131, length.out = 10)

flip_mod <- cr_age_size_models$Phi.SEASON.sex_hatch_T_crflip

make_flip_pred <- function(zcr.flip) {
flip.pred=covariate.predictions(flip_mod, data=data.frame(cr.flip = zcr.flip, did.cr = 1), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
}
make_flip_pred_cronly <- function(zcr.flip) {
flip.pred=covariate.predictions(flip_mod, data=data.frame(cr.flip = zcr.flip), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
}

flip_preds <- map(cr_flip_span, make_flip_pred)
  

# daily survival probability plot - flipper model ----
flip_preds_daily <- map_df(flip_preds, get_daily_pred) %>% 
  ungroup() 

# if creched only
flip_preds_daily <- flip_preds_daily %>% 
  rename(cr.flip = covdata)


flip_preds_daily %>% 
  filter(between(par.index, 1, 48) | between(par.index, 1177, 1177 + 47) | between(par.index, 2353, 2353 + 47) | between(par.index, 3529, 3529 + 47)) %>% 
  ggplot(group = cr.flip)+
	geom_line(aes(x = day, y = estimate, color = as.factor(cr.flip)), size=1) +
	#geom_line(aes(x = day, y = ucl, color = flipper.length), linetype = "dotted", size = 1) +
	#geom_line(aes(x = day, y = lcl, color = flipper.length), linetype = "dotted", size = 1) +
	scale_colour_discrete(name="Creching\nflipper (mm)")+
	#geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(sex~SEASON)+
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

ggsave(here("../MSdocs/daily_survival_by_flipper_cr_chicks.png")) 
  

# overall survival probability plot - flipper model ----

flip_overall_pred <- flip_preds_daily %>% 
    add_sex_seas() %>%
#  rename(cr.flip = covdata) %>% # if cr only
    group_by(SEASON, sex, cr.flip) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON),
         sex = as.character(sex))
  

flip_overall_se <- map_df(flip_preds, get_overall_se) %>% 
  ungroup() 
# if creched only
flip_overall_se <- flip_overall_se %>% 
  rename(cr.flip = covdata)


flip_overall_pred_se <- full_join(flip_overall_pred, flip_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))

# change title to reflect which results
  flip_overall_pred_se %>% 
    ggplot() +
    geom_line(aes(cr.flip, overall.surv)) +
    geom_ribbon(aes(x = cr.flip, ymin = lci, ymax = uci), alpha = 0.2) +
    labs(x = "Creching flipper length (mm)",
         y = "Probability of surviving to flegde",
         title = "Creched chicks") +
    theme_bw() +
    facet_grid(sex~SEASON)
  
 ggsave(here("../MSdocs/survival_by_flipper_cr_chicks.png")) 
  
  
# make mass model plot ----
  summary(ana_table$cr.mass)
  zcr.mass <- seq(100, 2900, length.out = 100)

  mass_mod <- cr_age_size_models$Phi.SEASON.sex_hatch_T_crmass

make_mass_pred <- function(zcr.mass) {
  mass.pred=covariate.predictions(mass_mod, data=data.frame(cr.mass=zcr.mass), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
}

mass_preds <- map(zcr.mass, make_mass_pred)

# daily survival probability plot - mass model ----
mass_preds_daily <- map_df(mass_preds, get_daily_pred) %>% 
  ungroup() %>% 
  rename(mass = covdata)

mass_preds_daily %>% 
ggplot()+
	geom_line(aes(x = day, y = estimate, color = as.factor(mass)), size=1) +
	#geom_line(aes(x = day, y = ucl, color = pred), linetype = "dotted", size = 1) +
	#geom_line(aes(x = day, y = lcl, color = pred), linetype = "dotted", size = 1) +
	scale_colour_discrete(name="Creching\nmass (g)")+
	#geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(sex~SEASON)+
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

mass_overall_pred <- map_df(mass_preds, get_overall_pred) %>% 
  ungroup()
  

mass_overall_se <- map_df(mass_preds, get_overall_se) %>% 
  ungroup()


mass_overall_pred_se <- full_join(mass_overall_pred, mass_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se)) %>% 
  rename(mass = covdata)


  mass_overall_pred_se %>% 
    ggplot() +
    geom_line(aes(mass, overall.surv)) +
    geom_ribbon(aes(x = mass, ymin = lci, ymax = uci), alpha = 0.2) +
    labs(x = "Creching mass (g)",
         y = "Probability of surviving to flegde",
         title = "Creched chicks") +
    theme_bw() +
    facet_grid(sex~SEASON)
  
  
  
 ggsave(here("../MSdocs/survival_by_mass_cr_chicks.png")) 
  




