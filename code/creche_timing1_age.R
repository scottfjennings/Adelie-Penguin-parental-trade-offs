

library(tidyverse)
library(here)
library(lme4)
library(AICcmodavg)
library(grid)
options(scipen = 999)

source(here("code/utilities.R"))


# read data ----
#  ana.table.full.csv is created by the sequence of numbered code files in THESIS\thesis_data_work\code_files\analysis_data_prep
#ana.table=read.csv("data/ana.table.full.csv")
ana.table <- readRDS(here("data/ana_table_full"))
#ana.table.old=read.csv("data/ana.table.full.old.csv")%>%  rename(CHICK_ID = id)
# data management ----

data = ana.table %>% 
  filter(!is.na(cr.seasday), !is.na(weight.slope40), !is.na(sex)) %>% 
  mutate(SEASON = as.factor(SEASON))


# analysis ----
# first creching age as response variable ----
# creching age model selection step 1 ----

cr_age_step1_mods <- list(
sex = lm(cr.age ~	sex, data = data),
year = lm(cr.age ~ SEASON, data = data),
hatch = lm(cr.age ~	resid.hatch, data = data),
sex_year = lm(cr.age~	sex + SEASON, data = data),
sex_hatch = lm(cr.age ~	sex + resid.hatch, data = data),
year_hatch = lm(cr.age ~	SEASON + resid.hatch, data = data),
sex_year_hatch = lm(cr.age ~ sex + SEASON + resid.hatch, data = data),
sex.year = lm(cr.age~	sex * SEASON, data = data),
sex.year_hatch = lm(cr.age ~ sex * SEASON + resid.hatch, data = data),
int = lm(cr.age ~ 1, data = data))

cr_age_step1_mods$model.table <- aictab(cr_age_step1_mods, names(cr_age_step1_mods))



saveRDS(cr_age_step1_mods, here("fitted_models/cr_timing/cr_age_step1_mods"))

# uninformative parms in supported step 1 models? ----

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

cr_age_step1_mods <- readRDS(here("fitted_models/cr_timing/cr_age_step1_mods"))

cr_age_mods <- cr_age_step1_mods[-length(cr_age_step1_mods)]

cr_age_informative <- map2_df(cr_age_mods, names(cr_age_mods), lm_parm_informative)

cr_age_informative <- cr_age_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm)) %>% 
  full_join(., cr_age_step1_mods$model.table %>%
              data.frame() %>%  
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_age_informative_wide <- cr_age_informative %>%
  pivot_wider(id_cols = c(mod.name, mod.call), names_from = parm, values_from = informative85) %>% 
  mutate(across(c("SEASON", "hatch", "sex", "sex.SEASON"), ~replace_na(., TRUE))) %>% 
  mutate(mod.call = sub("resid\\.", "", mod.call),
         mod.call = ifelse(sex.SEASON == FALSE, str_replace(mod.call, "sex \\* SEASON", paste("(sex \\* SEASON)", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(sex == FALSE, str_replace(mod.call, "sex", paste("sex", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(hatch == FALSE, str_replace(mod.call, "hatch", paste("hatch", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(SEASON == FALSE, str_replace(mod.call, "SEASON", paste("SEASON", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(mod.name == "int", "Intercept only", mod.call),
         mod.call = mod_call_to_structure(mod.call)) %>% 
  select(Modnames = mod.name, model.structure = mod.call)

# append this column with model.structure to the step 1 model table

cr_age_step1_mods$model.table <- cr_age_step1_mods$model.table %>%
  data.frame() %>% 
  full_join(cr_age_informative_wide)

saveRDS(cr_age_step1_mods, here("fitted_models/cr_timing/cr_age_step1_mods"))

# creching age model selection step 2 ----


cr_age_step2_mods <- list(
  year_hatch = cr_age_step1_mods$year_hatch,
  year_hatch_mass = lm(cr.age ~	SEASON + resid.hatch + weight.slope40, data = data),
  year_hatch_flip = lm(cr.age ~	SEASON + resid.hatch + flipper.slope40, data = data),
  year_hatch_tib = lm(cr.age ~	SEASON + resid.hatch + tibiotar.slope35, data = data))

cr_age_step2_mods$model.table = aictab(cr_age_step2_mods, names(cr_age_step2_mods))


cr_age2_mods <- cr_age_step2_mods[-length(cr_age_step2_mods)]

cr_age2_informative <- map2_df(cr_age2_mods, names(cr_age2_mods), lm_parm_informative)


cr_age2_informative <- cr_age2_informative %>% 
  mutate(mod.name = gsub("cr_age_step2_mods\\$", "", mod.name),
         parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         mod.call = sub("resid.", "", mod.call)) %>% 
  full_join(., cr_age_step2_mods$model.table %>%
              data.frame() %>% 
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_age2_informative_wide <- cr_age2_informative %>%
  pivot_wider(id_cols = c(mod.name, mod.call), names_from = parm, values_from = informative85) %>% 
  mutate(across(c("SEASON", "hatch", "weight.slope40", "flipper.slope40", "tibiotar.slope35"), ~replace_na(., TRUE))) %>% 
  mutate(mod.call = ifelse(SEASON == FALSE, str_replace(mod.call, "SEASON", paste("SEASON", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(hatch == FALSE, str_replace(mod.call, "hatch", paste("hatch", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(weight.slope40 == FALSE, str_replace(mod.call, "weight.slope40", paste("weight.slope40", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(flipper.slope40 == FALSE, str_replace(mod.call, "flipper.slope40", paste("flipper.slope40", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(tibiotar.slope35 == FALSE, str_replace(mod.call, "tibiotar.slope35", paste("tibiotar.slope35", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(mod.name == "int", "Intercept only", mod.call),
         mod.call = mod_call_to_structure(mod.call)) %>% 
  select(Modnames = mod.name, model.structure = mod.call)



cr_age_step2_mods$model.table <- cr_age_step2_mods$model.table %>%
  data.frame() %>% 
  full_join(cr_age2_informative_wide)


saveRDS(cr_age_step2_mods, here("fitted_models/cr_timing/cr_age_step2_mods"))


# all models with resid.hatch are competitive (largest dAICc = 2.14), then a break in dAICc to 43.18 for 5th best.


zzz <- readRDS(here("fitted_models/cr_timing/cr_age_step1_mods")) %>%
    cr_timing_in_table() %>%
    full_join(., readRDS(here("fitted_models/cr_timing/cr_age_step1_uninformative")) %>% rename(Modnames = mod.name)) %>% 
    mutate(uninformative = ifelse(is.na(uninformative), 0, 1))


