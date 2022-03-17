
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


cr_flip_step1_mods <- list(
sex = lm(cr.flip ~	sex, data = data),
year = lm(cr.flip ~ SEASON, data = data),
hatch = lm(cr.flip ~	resid.hatch, data = data),
sex_year = lm(cr.flip~	sex + SEASON, data = data),
sex_hatch = lm(cr.flip ~	sex + resid.hatch, data = data),
year_hatch = lm(cr.flip ~	SEASON + resid.hatch, data = data),
sex_year_hatch = lm(cr.flip ~ sex + SEASON + resid.hatch, data = data),
sex.year = lm(cr.flip~	sex * SEASON, data = data),
sex.year_hatch = lm(cr.flip ~ sex * SEASON + resid.hatch, data = data),
int = lm(cr.flip ~ 1, data = data)
)

cr_flip_step1_mods$model.table <- aictab(cr_flip_step1_mods, names(cr_flip_step1_mods))



saveRDS(cr_flip_step1_mods, here("fitted_models/cr_timing/cr_flip_step1_mods"))


# mass.sex_year best supported
# uninformative parms in supported step 1 models? ----

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

cr_flip_step1_mods <- readRDS(here("fitted_models/cr_timing/cr_flip_step1_mods"))


cr_flip_mods <- cr_flip_step1_mods[-length(cr_flip_step1_mods)]

cr_flip_informative <- map2_df(cr_flip_mods, names(cr_flip_mods), lm_parm_informative)


cr_flip_informative <- cr_flip_informative %>% 
  mutate(mod.name = gsub("cr_flip_step1_mods\\$", "", mod.name),
         parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         mod.call = sub("resid.", "", mod.call)) %>% 
  full_join(., cr_flip_step1_mods$model.table %>%
              data.frame() %>% 
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_flip_informative_wide <- cr_flip_informative %>%
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


cr_flip_step1_mods$model.table <- cr_flip_step1_mods$model.table %>% 
  data.frame() %>% 
  full_join(., cr_flip_informative_wide) 

# cr_flip_step1_mods$model.table %>% view()

saveRDS(cr_flip_step1_mods, here("fitted_models/cr_timing/cr_flip_step1_mods"))
         
# creching mass model selection step 2 ----

cr_flip_step2_mods <- list(
  year_hatch = cr_flip_step1_mods$year_hatch,
  year_hatch_flip = lm(cr.flip ~	SEASON + resid.hatch + flipper.slope40, data = data),
  year_hatch_age = lm(cr.flip ~	SEASON + resid.hatch + cr.age, data = data),
  year_hatch_flip_age = lm(cr.flip ~	SEASON + resid.hatch + flipper.slope40 + cr.age, data = data))


cr_flip_step2_mods$model.table = aictab(cr_flip_step2_mods, names(cr_flip_step2_mods))

saveRDS(cr_flip_step2_mods, here("fitted_models/cr_timing/cr_flip_step2_mods"))


# check for uninformative parms in step 2 ----
cr_flip_step2_mods <- readRDS(here("fitted_models/cr_timing/cr_flip_step2_mods"))



cr_flip2_mods <- cr_flip_step2_mods[-length(cr_flip_step2_mods)]

cr_flip2_informative <- map2_df(cr_flip2_mods, names(cr_flip2_mods), lm_parm_informative)

cr_flip2_informative <- cr_flip2_informative %>% 
  mutate(mod.name = gsub("cr_flip_step2_mods\\$", "", mod.name),
         parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexM", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         mod.call = sub("resid.", "", mod.call)) %>% 
  full_join(., cr_flip_step2_mods$model.table %>%
              data.frame() %>% 
              select(mod.name = Modnames, Delta_AICc)) %>% 
  mutate(informative85 = ifelse(Delta_AICc > 5, TRUE, informative85))

cr_flip2_informative_wide <- cr_flip2_informative %>%
  mutate(parm = ifelse(mod.name == "int", "intercept", parm)) %>% 
  pivot_wider(id_cols = c(mod.name, mod.call), names_from = parm, values_from = informative85) %>% 
  mutate(across(c("SEASON", "hatch", "flipper.slope40", "cr.age"), ~replace_na(., TRUE))) %>% 
  mutate(mod.call = ifelse(hatch == FALSE, str_replace(mod.call, "hatch", paste("hatch", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(SEASON == FALSE, str_replace(mod.call, "SEASON", paste("SEASON", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(flipper.slope40 == FALSE, str_replace(mod.call, "flipper.slope40", paste("flipper.slope40", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(cr.age == FALSE, str_replace(mod.call, "cr.age", paste("cr.age", "\u2020", sep = "")), mod.call),
         mod.call = ifelse(mod.name == "int", "Intercept only", mod.call),
         mod.call = mod_call_to_structure(mod.call)) %>% 
  select(Modnames = mod.name, model.structure = mod.call)


cr_flip_step2_mods$model.table <- cr_flip_step2_mods$model.table %>% 
  data.frame() %>% 
  full_join(., cr_flip2_informative_wide) 

# cr_flip_step2_mods$model.table %>% view()

saveRDS(cr_flip_step2_mods, here("fitted_models/cr_timing/cr_flip_step2_mods"))

