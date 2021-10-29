
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
# creching mass model selection step 2 ----
# uninformative parms in supported step 1 models?
cr_flip_step1_mods <- readRDS(here("fitted_models/cr_timing/cr_flip_step1_mods"))

cr_flip_step1_mods$model.table %>%
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(Delta_AICc <= 5) %>% 
  mutate(uninform.assign = paste("lm_parm_informative(cr_flip_step1_mods$", Modnames, ")", sep = "")) %>% 
  summarise(uninform.assign = paste(uninform.assign, collapse = ", ")) %>% 
  mutate(uninform.assign = paste("cr_flip_informative <- rbind(", uninform.assign, ")", sep = ""))

cr_flip_informative <- rbind(lm_parm_informative(cr_flip_step1_mods$year_hatch), 
                             lm_parm_informative(cr_flip_step1_mods$sex_year_hatch), 
                             lm_parm_informative(cr_flip_step1_mods$sex.year_hatch))

cr_flip_informative <- cr_flip_informative %>% 
  mutate(mod.name = gsub("cr_flip_step1_mods\\$", "", mod.name))

cr_flip_informative_wide <- cr_flip_informative %>% 
  pivot_wider(id_cols = mod.name, names_from = parm, values_from = informative85) %>% 
  mutate(uninformative = "")


cr_flip_informative_wide <- edit(cr_flip_informative_wide)

cr_flip_informative_wide %>% 
  filter(uninformative == TRUE) %>%  
  select(mod.name, uninformative) %>% 
  saveRDS(here("fitted_models/cr_timing/cr_flip_step1_uninformative"))


cr_flip_informative_wide %>% 
  filter(uninformative != TRUE) %>%
  mutate(mod.assign = paste(mod.name, " = cr_flip_step1_mods$", mod.name, sep = "")) %>% 
  select(mod.assign)
         

cr_flip_step2_mods <- list(
  year_hatch = cr_flip_step1_mods$year_hatch,
  #
  year_hatch_flip = lm(cr.flip ~	SEASON + resid.hatch + flipper.slope40, data = data),
  #
  year_hatch_age = lm(cr.flip ~	SEASON + resid.hatch + cr.age, data = data),
  #
  year_hatch_flip_age = lm(cr.flip ~	SEASON + resid.hatch + flipper.slope40 + cr.age, data = data))


cr_flip_step2_mods$model.table = aictab(cr_flip_step2_mods, names(cr_flip_step2_mods))


saveRDS(cr_flip_step2_mods, here("fitted_models/cr_timing/cr_flip_step2_mods"))

