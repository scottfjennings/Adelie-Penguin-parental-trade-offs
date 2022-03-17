
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



# creching mass model selection step 1 ----


cr_mass_step1_mods <- list(
sex = lm(cr.mass ~	sex, data = data),
year = lm(cr.mass ~ SEASON, data = data),
hatch = lm(cr.mass ~	resid.hatch, data = data),
sex_year = lm(cr.mass~	sex + SEASON, data = data),
sex_hatch = lm(cr.mass ~	sex + resid.hatch, data = data),
year_hatch = lm(cr.mass ~	SEASON + resid.hatch, data = data),
sex_year_hatch = lm(cr.mass ~ sex + SEASON + resid.hatch, data = data),
sex.year = lm(cr.mass~	sex * SEASON, data = data),
sex.year_hatch = lm(cr.mass ~ sex * SEASON + resid.hatch, data = data),
int = lm(cr.mass ~ 1, data = data)
)

cr_mass_step1_mods$model.table <- aictab(cr_mass_step1_mods, names(cr_mass_step1_mods))



saveRDS(cr_mass_step1_mods, here("fitted_models/cr_timing/cr_mass_step1_mods"))


# mass.sex_year best supported
# uninformative parms in supported step 1 models? ----

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

cr_mass_step1_mods <- readRDS(here("fitted_models/cr_timing/cr_mass_step1_mods"))

cr_mass_mods <- cr_mass_step1_mods[-length(cr_mass_step1_mods)]

cr_mass_informative <- map2_df(cr_mass_mods, names(cr_mass_mods), lm_parm_informative)


cr_mass_informative <- cr_mass_informative %>% 
  mutate(mod.name = gsub("cr_mass_step1_mods\\$", "", mod.name),
         parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         mod.call = sub("resid.", "", mod.call)) %>% 
  full_join(., cr_mass_step1_mods$model.table %>%
              data.frame() %>% 
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_mass_informative_wide <- cr_mass_informative %>%
  mutate(parm = ifelse(mod.name == "int", "intercept", parm)) %>% 
  pivot_wider(id_cols = c(mod.name, mod.call), names_from = parm, values_from = informative85) %>% 
  mutate(across(c("SEASON", "hatch", "sex", "sex.SEASON"), ~replace_na(., TRUE))) %>% 
  mutate(mod.call = ifelse(sex.SEASON == FALSE, str_replace(mod.call, "sex \\* SEASON", paste("(sex \\* SEASON)", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(sex == FALSE, str_replace(mod.call, "sex", paste("sex", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(hatch == FALSE, str_replace(mod.call, "hatch", paste("hatch", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(SEASON == FALSE, str_replace(mod.call, "SEASON", paste("SEASON", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(mod.name == "int", "Intercept only", mod.call),
         mod.call = mod_call_to_structure(mod.call)) %>% 
  select(Modnames = mod.name, model.structure = mod.call)


cr_mass_step1_mods$model.table <- cr_mass_step1_mods$model.table %>% 
  data.frame() %>% 
  full_join(., cr_mass_informative_wide) 


saveRDS(cr_mass_step1_mods, here("fitted_models/cr_timing/cr_mass_step1_mods"))

# creching mass model selection step 2 ----
cr_mass_step1_mods <- readRDS(here("fitted_models/cr_timing/cr_mass_step1_mods"))


cr_mass_step2_mods <- list(
  sex_year_hatch = cr_mass_step1_mods$sex_year_hatch,
  year_hatch = cr_mass_step1_mods$year_hatch,
  #
  year_hatch_mass = lm(cr.mass ~	SEASON + resid.hatch + weight.slope40, data = data),
  sex_year_hatch_mass = lm(cr.mass ~ sex + SEASON + resid.hatch + weight.slope40, data = data),
  #
  year_hatch_age = lm(cr.mass ~	SEASON + resid.hatch + cr.age, data = data),
  sex_year_hatch_age = lm(cr.mass ~ sex + SEASON + resid.hatch + cr.age, data = data),
  #
  year_hatch_mass_age = lm(cr.mass ~	SEASON + resid.hatch + weight.slope40 + cr.age, data = data),
  sex_year_hatch_mass_age = lm(cr.mass ~ sex + SEASON + resid.hatch + weight.slope40 + cr.age, data = data))


cr_mass_step2_mods$model.table = aictab(cr_mass_step2_mods, names(cr_mass_step2_mods))

saveRDS(cr_mass_step2_mods, here("fitted_models/cr_timing/cr_mass_step2_mods"))

# uninformative parms in step 2? ----
cr_mass_step2_mods <- readRDS(here("fitted_models/cr_timing/cr_mass_step2_mods"))


cr_mass2_mods <- cr_mass_step2_mods[-length(cr_mass_step2_mods)]

cr_mass2_informative <- map2_df(cr_mass2_mods, names(cr_mass2_mods), lm_parm_informative)


cr_mass2_informative <- cr_mass2_informative %>% 
  mutate(mod.name = gsub("cr_mass_step2_mods\\$", "", mod.name),
         parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         mod.call = sub("resid.", "", mod.call)) %>% 
  full_join(., cr_mass_step2_mods$model.table %>%
              data.frame() %>% 
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_mass2_informative_wide <- cr_mass2_informative %>%
  mutate(parm = ifelse(mod.name == "int", "intercept", parm)) %>% 
  pivot_wider(id_cols = c(mod.name, mod.call), names_from = parm, values_from = informative85) %>% 
  mutate(across(c("SEASON", "hatch", "sex", "weight.slope40", "cr.age", "sex"), ~replace_na(., TRUE))) %>% 
  mutate(mod.call = ifelse(sex == FALSE, str_replace(mod.call, "sex", paste("sex", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(hatch == FALSE, str_replace(mod.call, "hatch", paste("hatch", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(SEASON == FALSE, str_replace(mod.call, "SEASON", paste("SEASON", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(weight.slope40 == FALSE, str_replace(mod.call, "weight.slope40", paste("weight.slope40", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(cr.age == FALSE, str_replace(mod.call, "cr.age", paste("cr.age", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(mod.name == "int", "Intercept only", mod.call),
         mod.call = mod_call_to_structure(mod.call)) %>% 
  select(Modnames = mod.name, model.structure = mod.call)


cr_mass_step2_mods$model.table <- cr_mass_step2_mods$model.table %>% 
  data.frame() %>% 
  full_join(., cr_mass2_informative_wide) 

# cr_mass_step2_mods$model.table %>% view()



saveRDS(cr_mass_step2_mods, here("fitted_models/cr_timing/cr_mass_step2_mods"))

