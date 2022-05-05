



library(tidyverse)
library(here)
library(RMark)
library(plotly)

options(scipen = 999)

ana.table <- readRDS(here("data/ana_table_full"))

penguins=convert.inp(here("data/mark_in.inp"), 
					group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), 
					covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56", 
					 "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64", 
					 "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72",
					 "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80",
					 "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", 
					 "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96",
					 "dayold97", "dayold98","in.cr50", "in.cr51", "in.cr52", "in.cr53", "in.cr54", "in.cr55", "in.cr56", 
					 "in.cr57", "in.cr58", "in.cr59", "in.cr60", "in.cr61", "in.cr62", "in.cr63", "in.cr64", 
					 "in.cr65", "in.cr66", "in.cr67", "in.cr68", "in.cr69", "in.cr70", "in.cr71", "in.cr72",
					 "in.cr73", "in.cr74", "in.cr75", "in.cr76", "in.cr77", "in.cr78", "in.cr79", "in.cr80",
					 "in.cr81", "in.cr82", "in.cr83", "in.cr84", "in.cr85", "in.cr86", "in.cr87", "in.cr88", 
					 "in.cr89", "in.cr90", "in.cr91", "in.cr92", "in.cr93", "in.cr94", "in.cr95", "in.cr96",
					 "in.cr97", "in.cr98", "res.htch", "cr.age", "fail.age", "mass.gr", "tib.gr", "flip.gr", "cr.mass", "cr.flip", "cr.tib"), use.comments = TRUE) 

# sample size for creche-timing ----
ana.table %>% 
  filter(!is.na(cr.seasday), !is.na(weight.slope40), !is.na(sex)) %>% 
  count()

ana.table %>% 
  filter(!is.na(cr.seasday), !is.na(weight.slope40), !is.na(sex)) %>% 
  count(SEASON, sex)

# sample size for overall survival ----
penguins %>% 
  data.frame() %>% 
  count()
penguins %>% 
  data.frame() %>% 
  count(SEASON)
penguins %>% 
  data.frame() %>% 
  count(SEASON, sex)

# sample size for growth rate ~ survival ----
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr)) %>% 
  count()
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr)) %>% 
  count(SEASON)
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr)) %>% 
  count(SEASON, sex)


# sample size for creche age/size ~ survival ----
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr), !is.na(cr.age)) %>% 
  count()
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr), !is.na(cr.age)) %>% 
  count(SEASON)
penguins %>% 
  data.frame() %>% 
  filter(!is.na(mass.gr), !is.na(cr.age)) %>% 
  count(SEASON, sex)


# chicks in full ana.table but not in penguins ----

anti_join(ana.table,
          penguins %>% 
            data.frame() %>% 
            rownames_to_column("CHICK_ID") %>% 
            mutate(CHICK_ID = trimws(CHICK_ID)) %>% 
            select(CHICK_ID)) %>% 
  view()


# creche ages by year ----
ana.table %>% 
  filter(!is.na(cr.age)) %>% 
  group_by(SEASON) %>% 
  summarise(min.cr = min(cr.age),
         mean.cr = mean(cr.age),
         max.cr = max(cr.age),
         se.cr = sd(cr.age)/sqrt(n()))
  
  
  
# correlation between hatch date and growth rates ----
cor.test(ana.table$weight.slope40, ana.table$hatch.seasday)
cor.test(ana.table$flipper.slope40, ana.table$hatch.seasday)
cor.test(ana.table$tibiotar.slope35, ana.table$hatch.seasday)

ana.table %>% 
  filter(!is.na(weight.slope40)) %>% 
  select(hatch.seasday, weight.slope40, flipper.slope40, tibiotar.slope35) %>% 
  pivot_longer(cols = c(weight.slope40, flipper.slope40, tibiotar.slope35), names_to = "morph", values_to = "growth") %>% 
  ggplot() +
  geom_point(aes(x = hatch.seasday, y = growth)) +
  facet_wrap(~morph, scales = "free_y")


# correlation between creche age and growth rates ----
cor.test(ana.table$weight.slope40, ana.table$cr.age)
cor.test(ana.table$flipper.slope40, ana.table$cr.age)
cor.test(ana.table$tibiotar.slope35, ana.table$cr.age)

ana.table %>% 
  filter(!is.na(weight.slope40)) %>% 
  select(hatch.seasday, weight.slope40, flipper.slope40, tibiotar.slope35) %>% 
  pivot_longer(cols = c(weight.slope40, flipper.slope40, tibiotar.slope35), names_to = "morph", values_to = "growth") %>% 
  ggplot() +
  geom_point(aes(x = hatch.seasday, y = growth)) +
  facet_wrap(~morph, scales = "free_y")


# plot to help describe subtle differences in overall survival ----

readRDS(here("fitted_models/survival/overall_survival_step2_mod_av_phi"))$estimates %>% 
  filter(occ.cohort == 1) %>% 
  mutate(time = as.numeric(as.character(time)),
         SEASON = ifelse(SEASON == 1213, "2012-13", "2013-14")) %>% 
  ggplot() +
	geom_line(aes(x = time, y = estimate, color = group))


# tib growth coef estimates from best growth ~ survival model ----

best_growth_surv <- readRDS(here("fitted_models/survival/growth_surv_models"))$Phi.SEASON_hatch_TT_tibgr

best_growth_surv_beta <- best_growth_surv$results$beta





# difference in survival between chicks with crech size 1600g (1970's Croz) and 1350 g (current study) ----
# read functions from output3_survival_figures.r
  zcr.mass <- c(1350, 1600)

  mass_mod <- readRDS(here("fitted_models/survival/cr_age_size_models"))$Phi.SEASON_hatch_TT_crmass

make_mass_pred <- function(zcr.mass) {
  mass.pred=covariate.predictions(mass_mod, data=data.frame(cr.mass=zcr.mass), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
}


mass_preds <- map(zcr.mass, make_mass_pred)


mass_preds_daily <- map_df(mass_preds, get_daily_pred) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata)


mass_overall_pred <- mass_preds_daily %>% 
    add_sex_seas() %>%
  filter(day > 68) %>% 
    group_by(SEASON, sex, cr.mass) %>% 
    summarise(overall.surv = prod(estimate)) %>% 
  ungroup() %>% 
  mutate(SEASON = as.character(SEASON),
         sex = as.character(sex))
  

mass_overall_se <- map_df(mass_preds, get_overall_se_creched) %>% 
  ungroup() %>% 
  rename(cr.mass = covdata)


mass_overall_pred_se <- full_join(mass_overall_pred, mass_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se))



mass_overall_pred_se <- full_join(mass_overall_pred, mass_overall_se) %>% 
  mutate(lci = overall.surv - (1.96 * se),
         uci = overall.surv + (1.96 * se)) 

mass_overall_pred_se %>% 
  filter(sex == "Male", SEASON == "1213") %>% 
  view()



# plots to help explain creche timing results ----
ana.table %>% 
  ggplot() +
  geom_point(aes(x = resid.hatch, y = cr.seasday, color = weight.slope40))

ana.table %>% 
  ggplot() +
  geom_point(aes(x = resid.hatch, y = cr.age, color = cr.mass, size = cr.mass)) +
  facet_grid(~SEASON)

ana.table %>% 
  select(resid.hatch, cr.age, cr.mass, weight.slope40) %>% 
  plot()

# refit best creche timing models with scales predictors


data = ana.table %>% 
  filter(!is.na(cr.seasday), !is.na(weight.slope40), !is.na(sex)) %>% 
  mutate(SEASON = as.factor(SEASON),
         cr.age.sc = scale(cr.age),
         weight.slope40.sc = scale(weight.slope40),
         flipper.slope40.sc = scale(flipper.slope40),
         tibiotar.slope35.sc = scale(tibiotar.slope35))

mass.age_growth_sex_year.sc	= lm(cr.mass ~	cr.age.sc + weight.slope40.sc + sex + SEASON, data = data)

flip.age_growth_sex_year.sc	= lm(cr.mass ~	cr.age.sc + flipper.slope40.sc + sex + SEASON, data = data)

tib.age_growth_sex_year.sc	= lm(cr.mass ~	cr.age.sc + tibiotar.slope35.sc + sex + SEASON, data = data)


# model averaged overall survival
step2_mod_av <- readRDS(here("fitted_models/survival/overall_survival_step2_mod_av_phi"))

# no cohort effects, all estimates equal for each cohort, so just take estimates for first cohort
step2_mod_av_est <- filter(step2_mod_av$estimates, occ.cohort == 1)

f12_se <- deltamethod.special("prod", step2_mod_av$estimates$estimate[1:48], step2_mod_av$vcv.real[1:48, 1:48])
f13_se <- deltamethod.special("prod", step2_mod_av$estimates$estimate[49:96], step2_mod_av$vcv[49:96, 49:96])
m12_se <- deltamethod.special("prod", step2_mod_av$estimates$estimate[97:144], step2_mod_av$vcv[97:144, 97:144])
m13_se <- deltamethod.special("prod", step2_mod_av$estimates$estimate[145:192], step2_mod_av$vcv[145:192, 145:192])

se <- rbind(f12_se, f13_se, m12_se, m13_se) %>% 
  data.frame() %>% 
  rename(se = 1) %>% 
  rownames_to_column("seas.sex") %>% 
  mutate(seas.sex = case_when(seas.sex == "f12_se" ~ "1213Female",
                              seas.sex == "f13_se" ~ "1314Female",
                              seas.sex == "m12_se" ~ "1213Male",
                              seas.sex == "m13_se" ~ "1314Male"))

# plot 
step2_mod_av_est %>% 
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
  

# what growth rates and creche ages are required to achieve different survival curves? ----
dat <- filter(ana.table, weight.slope40 > 0, !is.na(weight.slope40), !is.na(cr.age), !is.na(sex))
cr_mass_mod <- readRDS(here("fitted_models/cr_timing/cr_mass_step2_mods"))$mass.age_growth_sex_year

mean_cr_age_newdat <- expand.grid(sex = c("M", "F"),
                                  SEASON = factor(c(1213, 1314)),
                                  cr.age = seq(11, 26, length.out = 5),
                                  weight.slope40 = seq(min(dat$weight.slope40), max(dat$weight.slope40), length.out = 10))

mean_cr_age_pred = predict(cr_mass_mod, mean_cr_age_newdat, interval = "confidence") %>% 
  data.frame() %>% 
  cbind(., mean_cr_age_newdat)
  
mean_cr_age_pred %>%  
ggplot() +
  geom_line(aes(x = weight.slope40, y = fit, color = as.character(cr.age))) +
  facet_grid(SEASON ~ sex) 


cr_mass_mod_coef <- coef(cr_mass_mod) %>% 
  data.frame() %>% 
  rename(coef = 1) %>% 
  rownames_to_column("varb") %>% 
  mutate(varb = gsub("weight.slope40|flipper.slope40|tibiotar.slope35", "growth", varb)) %>% 
  pivot_wider(names_from = varb, values_from = coef) %>% 
  rename(int = "(Intercept)")


# 1350 = -306.717542 + (20 * 60.165085) + (x * 9.203708) + (1 * 69.223371) + (0 * -256.956913)
# 1350 = -306.717542 + 1203.302 + (x * 9.203708) + 69.223371
# 1350 = 965.8078 + (x * 9.203708)
# 384.1922 = x * 9.203708
# 41.7432 = x

predict_growth <- function(zmod_coef, zcr.size, zcr.age, zsex, zseas) {
  # predict growth rate needed to achieve a given creche size (zcr.size) at mean creche age (20)
(zcr.size - (zmod_coef$int + (zcr.age * zmod_coef$cr.age) + (zsex * zmod_coef$sexM) + (zseas * zmod_coef$SEASON1314)))/zmod_coef$growth
}

predict_growth(cr_mass_mod_coef, 1800, 21, 1, 0)

growth_wrapper <- function(zcr.size) {
  m12 <- predict_growth(cr_mass_mod_coef, zcr.size, 1, 0) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "M_1213")
  f12 <- predict_growth(cr_mass_mod_coef, zcr.size, 0, 0) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "F_1213")  
  m13 <- predict_growth(cr_mass_mod_coef, zcr.size, 1, 1) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "M_1314")
  f13 <- predict_growth(cr_mass_mod_coef, zcr.size, 0, 1) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "F_1314")
  out.growth = rbind(m12, f12, m13, f13) %>% 
    mutate(cr.size = zcr.size)
  }

zzz <- map_df(zcr.mass, growth_wrapper)

# 1350 = -306.717542 + (x * 60.165085) + (60 * 9.203708) + (1 * 69.223371) + (0 * -256.956913)
# 1350 = -306.717542 + (x * 60.165085) + 552.2225 + 69.223371
# 1350 = 314.7283 + (x * 9.203708)
# 1035.272 = x * 60.165085
# 41.7432 = x
predict_cr_age <- function(zmod_coef, zcr.size, zgrowth, zsex, zseas) {
  # predict growth rate needed to achieve a given creche size (zcr.size) at mean creche age (20)
(zcr.size - (zmod_coef$int + (zgrowth * zmod_coef$growth) + (zsex * zmod_coef$sexM) + (zseas * zmod_coef$SEASON1314)))/zmod_coef$cr.age
}

predict_cr_age(cr_mass_mod_coef, 1800, 82, 1, 0)


-306.717542 + (21 * 60.165085) + (68 * 9.203708) + (1 * 69.223371) + (0 * -256.956913)


cr_age_wrapper <- function(zcr.size) {
  m12 <- predict_cr_age(cr_mass_mod_coef, zcr.size, 1, 0) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "M_1213")
  f12 <- predict_cr_age(cr_mass_mod_coef, zcr.size, 0, 0) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "F_1213")  
  m13 <- predict_cr_age(cr_mass_mod_coef, zcr.size, 1, 1) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "M_1314")
  f13 <- predict_cr_age(cr_mass_mod_coef, zcr.size, 0, 1) %>% 
    data.frame() %>% 
    rename(growth = 1) %>% 
    mutate(sex.seas = "F_1314")
  out.cr.age = rbind(m12, f12, m13, f13) %>% 
    mutate(cr.size = zcr.size)
  }


# want to get estimates for guarded survival up until the last chick creched each year, and only get creched survival estimates from the first creched chick onward
ana_table %>% 
  group_by(SEASON) %>% 
  summarise(min.cr = min(cr.seasday, na.rm = TRUE), 
            max.cr = max(cr.seasday, na.rm = TRUE), 
            mean.cr = mean(cr.seasday, na.rm = TRUE))


cr_mod <- readRDS(here("fitted_models/survival/step2"))$phi.year_T_hatch_incr_p.sat

cr_mod$results$beta[1:5,]



# check which PIMS correspond to the desired days
PIMS(cr_mod, parameter = "Phi", simplified = FALSE)

# mean creche day indices are 19 and 1175

guard_data12 = data.frame(in.cr50 = 0:1, 
                          in.cr51 = 0:1, 
                          in.cr52 = 0:1, 
                          in.cr53 = 0:1, 
                          in.cr54 = 0:1, 
                          in.cr55 = 0:1, 
                          in.cr56 = 0:1, 
                          in.cr57 = 0:1, 
                          in.cr58 = 0:1, 
                          in.cr59 = 0:1, 
                          in.cr60 = 0:1, 
                          in.cr61 = 0:1, 
                          in.cr62 = 0:1, 
                          in.cr63 = 0:1, 
                          in.cr64 = 0:1, 
                          in.cr65 = 0:1, 
                          in.cr66 = 0:1, 
                          in.cr67 = 0:1, 
                          in.cr68 = 0:1, 
                          in.cr69 = 0:1, 
                          in.cr70 = 0:1, 
                          in.cr71 = 0:1, 
                          in.cr72 = 0:1, 
                          in.cr73 = 0:1, 
                          in.cr74 = 0:1, 
                          in.cr75 = 0:1, 
                          in.cr76 = 0:1, 
                          in.cr77 = 0:1, 
                          in.cr78 = 0:1, 
                          in.cr79 = 0:1, 
                          in.cr80 = 0:1, 
                          in.cr81 = 0:1, 
                          in.cr82 = 0:1, 
                          in.cr83 = 0:1, 
                          in.cr84 = 0:1, 
                          in.cr85 = 0:1, 
                          in.cr86 = 0:1, 
                          in.cr87 = 0:1, 
                          in.cr88 = 0:1, 
                          in.cr89 = 0:1, 
                          in.cr90 = 0:1, 
                          in.cr91 = 0:1, 
                          in.cr92 = 0:1, 
                          in.cr93 = 0:1, 
                          in.cr94 = 0:1, 
                          in.cr95 = 0:1, 
                          in.cr96 = 0:1, 
                          in.cr97 = 0:1)




gc12 =covariate.predictions(cr_mod, data=guard_data12, indices = c(1:48), alpha=0.025)


gc13 =covariate.predictions(cr_mod, data=guard_data12, indices = c(1177:1224), alpha=0.025)


gc12_est <- gc12$estimates %>% 
  data.frame() %>% 
  select(-contains("in.cr"), in.cr = in.cr50) %>% 
  mutate(year = "2012", 
         day = par.index + 49) %>% 
  filter((in.cr == 0 & day < 75) | (in.cr == 1 & day > 59))

filter(gc12_est, day == 68) %>% 
  select(estimate, se, lcl, ucl, in.cr, year) %>% 
  mutate(across(c(estimate, se, lcl, ucl), ~round(., 3))) %>% 
  view()

gc13_est <- gc13$estimates %>% 
  data.frame() %>% 
  select(-contains("in.cr"), in.cr = in.cr50) %>% 
  mutate(year = "2013", 
         day = par.index - 1127) %>% 
  filter((in.cr == 0 & day < 76) | (in.cr == 1 & day > 62))


gc_est <- bind_rows(gc12_est, gc13_est) %>% 
    mutate(lcl = ifelse(lcl < 0, 0, lcl),
           ucl = ifelse(ucl > 1, 1, ucl),
           in.cr = ifelse(in.cr == 0, "guarded", "creched"))

gc_est %>% 
    ggplot(group = as.factor(in.cr)) +
    geom_line(aes(day, estimate, color = as.factor(in.cr))) +
    geom_ribbon(aes(x = day, ymin = lcl, ymax = ucl, fill = as.factor(in.cr)), alpha = 0.2) +
    theme_bw() +
    facet_grid(~year) +
  labs(color = "",
       fill = "")


# 
flip_cr_mod <- readRDS(here("fitted_models/survival/step4_growth_surv"))$phi.year_T_hatch_tibgr_p.sat

flip_cr_mod$results$beta[1:5,]

