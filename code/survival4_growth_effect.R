


#ACR computer
MarkPath="C:/Program Files (x86)/MARK"
library(RMark)
library(lme4) 
library(AICcmodavg)
library(grid)
library(ggplot2)
library(tidyverse)
library(R2ucare)
library(here)
options(scipen = 999)

source(here("code/utilities.R"))
# setting working directory only to control location for all the Mark output files
setwd(here("mark_output"))


# this part of the analysis needs to be separate from step 3 because it uses the subset of chicks for which I could extimate growth rates



# Nov 1 is season day 1 both years of this study. so begin.time = 50 makes analysis start on
# as.Date("2012-11-01") + 49
# "2012-12-20"
# as.Date("2013-11-01") + 49
# "2013-12-20"

# mean hatch day across both years is Dec 19, so start models on mean 1 day old

# the following script is adapted for my data from appendix C of the Mark book

# CJS analysis of penguin chick survival

# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/

## !!!NOTE 8/26/21- inp file now almost entirely created by code, including in.cr fields. See survival1_make_inp.R for code and notes 

# data ----
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



# filter out chicks that didn't survive long enough for good growth rate estimates.
penguins <- penguins %>% 
  filter(!is.na(mass.gr), !is.na(flip.gr), !is.na(tib.gr))
	

penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 50) 
penguins.ddl = make.design.data(penguins.process)

# growth_effect_pims <- penguins.ddl$Phi



run.models<-function(phi.stru, p.stru = "SEASON * sex + time") {
  phi.stru.list <- list(formula = formula(paste("~", phi.stru)))
  p.stru.list <- list(formula = formula(paste("~", p.stru)))
zmod <- mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru.list, p = p.stru.list), output = FALSE, chat = 1.25)
return(zmod)
}




# evaluating effect of growth rates on survival

# helper code to generate actual model fitting code - NO RUN ----
step1.4mods <- readRDS(here("fitted_models/survival/step1.4"))
step1.4mods$model.table %>% view()
step4_cand_set_base <- step1.4mods %>% 
  model.table() %>% 
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(DeltaQAICc <= 5) %>% 
  left_join(., names(step1.4mods) %>%
              data.frame() %>%
              rownames_to_column("mod.num") %>%
              rename(mod.name = 2)) %>% 
  select(mod.name, Phi, DeltaQAICc) %>% 
  mutate(Phi = gsub("~", "", Phi))

# best from step 1.3 - need to refit these models to the new data (can't just use the step 1.3 models)
step4_cand_set_base %>% 
  mutate(mod.assign = paste(mod.name, " <- run.models(phi.stru = \"", Phi, "\")", sep = "")) %>%  
  select(mod.assign)
# now adding in growth variables.
# best from step 1.3 plus additive mass growth rate
step4_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_massgr_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + mass.gr\")", sep = "")) %>% 
  select(mod.assign) 
# best from step 1.3 plus additive flipper growth
step4_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_flipgr_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + flip.gr\")", sep = "")) %>% 
  select(mod.assign)
# best from step 1.3 plus additive tibiotarsus growth
step4_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_tibgr_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + tib.gr\")", sep = "")) %>% 
  select(mod.assign)

# make code for collect.models call
rbind(step4_cand_set_base %>% 
  select(mod.name),
# now adding in time varying variables. these always need to have time, either time+ or time:
# best from step 1.3 plus additive time and time varying creched.	#same intercept, diff slope for in.cr
data.frame(mod.name = "#massgr"),
step4_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_massgr_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus additive time and time varying age. same intercept, diff slope for dayold
data.frame(mod.name = "#flipgr"),
step4_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_flipgr_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus interaction of time and time varying creched.	#same slope, diff intercept for in.cr
data.frame(mod.name = "#tibgr"),
step4_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_tibgr_p.", mod.name)) %>% 
  select(mod.name)) %>% 
  summarise(collect.call = paste(mod.name, collapse = "\", \""))

# now fitting models ----
# best step1.3
phi.year_TT_hatch_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch")
phi.year_T_hatch_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch")
phi.year_hatch_p.sat <- run.models(phi.stru = "SEASON + res.htch")
# best step1.3 plus mass growth
phi.year_TT_hatch_massgr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + mass.gr")
phi.year_T_hatch_massgr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + mass.gr")
phi.year_hatch_massgr_p.sat <- run.models(phi.stru = "SEASON + res.htch + mass.gr")
# best step1.3 plus flipper growth
phi.year_TT_hatch_flipgr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + flip.gr")
phi.year_T_hatch_flipgr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + flip.gr")
phi.year_hatch_flipgr_p.sat <- run.models(phi.stru = "SEASON + res.htch + flip.gr")
# best step1.3 plus tibiotarsus growth
phi.year_TT_hatch_tibgr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + tib.gr")
phi.year_T_hatch_tibgr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + tib.gr")
phi.year_hatch_tibgr_p.sat <- run.models(phi.stru = "SEASON + res.htch + tib.gr")


step4_growth_surv <- collect.models(lx = c("phi.year_TT_hatch_p.sat", "phi.year_T_hatch_p.sat", "phi.year_hatch_p.sat", 
                                           #massgr
                                           "phi.year_TT_hatch_massgr_p.sat", "phi.year_T_hatch_massgr_p.sat", "phi.year_hatch_massgr_p.sat",
                                           #flipgr
                                           "phi.year_TT_hatch_flipgr_p.sat", "phi.year_T_hatch_flipgr_p.sat", "phi.year_hatch_flipgr_p.sat",
                                           #tibgr
                                           "phi.year_TT_hatch_tibgr_p.sat", "phi.year_T_hatch_tibgr_p.sat", "phi.year_hatch_tibgr_p.sat"))


saveRDS(step4_growth_surv, here("fitted_models/survival/step4_growth_surv"))

# check for uninformative parms in step 4 ----

step4 <- readRDS(here("fitted_models/survival/step4_growth_surv"))

step4_mods <- step4[-length(step4)]

step4_informative <- map2_df(step4_mods, "Phi", cjs_parm_informative)

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

# p_informative has a row for each parameter in each model, but note it isn't meaningful to check each time varying parm for un vs. informative so these are excluded
step4_informative <- step4_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexMale", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         parm = sub("log\\(Time \\+ 1\\)", "lnT", parm),
         parm = sub("I\\(Time\\^2\\)", "TT", parm),
         parm = sub("Time", "T", parm)) %>%
  full_join(., step4$model.table %>% 
              data.frame() %>% 
              select(model.name = model, DeltaQAICc)) %>% 
  mutate(informative85 = ifelse(DeltaQAICc > 5, TRUE, informative85),
         parm = ifelse(grepl("~1", model.name), "p.intercept", parm),
         parm = ifelse(grepl("p\\(\\~time:in.cr\\)", model.name), "p.cr", parm))


         

step4_informative_wide <- step4_informative %>% 
           filter(!is.na(parm)) %>% 
  pivot_wider(id_cols = model.name, names_from = parm, values_from = informative85) %>% 
  mutate(across(c("Phi.T", "Phi.TT", "Phi.res.htch", "Phi.mass.gr",  "Phi.flip.gr",  "Phi.tib.gr"), ~replace_na(., TRUE))) %>% 
  mutate(phi.stru = gsub("\\)p\\(\\~SEASON \\* sex \\+ time\\)", "", model.name),
         phi.stru = gsub(":", ".", phi.stru),
         phi.stru = sub("resid.", "", phi.stru),
         phi.stru = sub("sexMale", "sex", phi.stru),
         phi.stru = sub("SEASON1314", "SEASON", phi.stru),
         phi.stru = sub("Time \\+ I\\(Time\\^2\\)", "TT", phi.stru),
         phi.stru = sub("Time", "T", phi.stru),
         phi.stru = sub("time", "t", phi.stru),
         phi.stru = sub("res.htch", "Hatch date", phi.stru),
         phi.stru = sub("mass.gr", "Mass", phi.stru),
         phi.stru = sub("flip.gr", "Flipper", phi.stru),
         phi.stru = sub("tib.gr", "Tibio", phi.stru),
         phi.stru = sub("Phi\\(\\~", "", phi.stru),
         phi.stru = ifelse(Phi.SEASON == FALSE, str_replace(phi.stru, "SEASON", paste("SEASON", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.T == FALSE & !grepl("TT", phi.stru), str_replace(phi.stru, "T", paste("T", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.TT == FALSE, str_replace(phi.stru, "TT", paste("TT", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.mass.gr == FALSE, str_replace(phi.stru, "Mass", paste("Mass", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.flip.gr == FALSE, str_replace(phi.stru, "Flipper", paste("Flipper", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.tib.gr == FALSE, str_replace(phi.stru, "Tibio", paste("Tibio", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(phi.stru == "1", "Intercept only", phi.stru),
         phi.stru = mod_call_to_structure(phi.stru)) %>% 
  select(model = model.name, phi.stru)

step4$model.table <- step4$model.table %>% 
  data.frame() %>% 
  full_join(step4_informative_wide)

step4$model.table %>% view()


saveRDS(step4, here("fitted_models/survival/step4"))


