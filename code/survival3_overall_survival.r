


#need to tell R that Mark is not stored in C:\Program Files

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



# Nov 1 is season day 1 both years of this study. so begin.time = 50 makes analysis start on
# as.Date("2012-11-01") + 49
# "2012-12-20"
# as.Date("2013-11-01") + 49
# "2013-12-20"

# mean hatch day across both years is Dec 19, so start models on mean 1 day old

# the following script is adapted for my data from appendix C of the Mark book

# CJS analysis of penguin chick survival

# data ----
# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/

## !!!NOTE 8/26/21- inp file now almost entirely created by code, including in.cr fields. See CH4_make_inp.R for code and notes 


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



# Next create the processed dataframe and the design data. We’ll use a group 
# variable for colony so it can be used in the set of models for Phi. Factor 
# variables (covariates with a small finite set of values) are best handled by using 
# them to define groups in the data. 
# need to set begin.time to 51 because RMark considers the first dayold to be the first release dayold, but really I want to consider recaptures beginning on dayold 52.
penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 50) 
penguins.ddl = make.design.data(penguins.process)




# fitting models ----

# strategy is to find the best structure on p and phi variables we want to account for first, then move on to main variables of interest.
# first find the best p structure while holding phi at most saturated (step 1.1), 
# then find best basic phi structure while holding p at most saturated (step 1.2)
# then exclude and uninformative parameters in competitive p and phi structures (step 1.3)
# and recollect the models with only informative parms (step 1.4)
# then finally fit models with parameters of interest (step 2)
# step 1.1 Determine p structure ----

phi.sat = "SEASON * sex + time"

phi.sat_p.year.sex <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex")
phi.sat_p.year_sex <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex")
phi.sat_p.sex <- run.models(phi.stru = phi.sat, p.stru = "sex")
phi.sat_p.year <- run.models(phi.stru = phi.sat, p.stru = "SEASON")
# add t
phi.sat_p.year.sex_t <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + time")
phi.sat_p.year_sex_t <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + time")
phi.sat_p.sex_t <- run.models(phi.stru = phi.sat, p.stru = "sex + time")
phi.sat_p.year_t <- run.models(phi.stru = phi.sat, p.stru = "SEASON + time")
# add T
phi.sat_p.year.sex_T <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + Time")
phi.sat_p.year_sex_T <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + Time")
phi.sat_p.sex_T <- run.models(phi.stru = phi.sat, p.stru = "sex + Time")
phi.sat_p.year_T <- run.models(phi.stru = phi.sat, p.stru = "SEASON + Time")
# add TT
phi.sat_p.year.sex_TT <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + Time + I(Time^2)")
phi.sat_p.year_sex_TT <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + Time + I(Time^2)")
phi.sat_p.sex_TT <- run.models(phi.stru = phi.sat, p.stru = "sex + Time + I(Time^2)")
phi.sat_p.year_TT <- run.models(phi.stru = phi.sat, p.stru = "SEASON + Time + I(Time^2)")
#add lnT
phi.sat_p.year.sex_lnT <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + log(Time+1)")
phi.sat_p.year_sex_lnT <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + log(Time+1)")
phi.sat_p.sex_lnT <- run.models(phi.stru = phi.sat, p.stru = "sex + log(Time+1)")
phi.sat_p.year_lnT <- run.models(phi.stru = phi.sat, p.stru = "SEASON + log(Time+1)")
# add creched _ I don't think this is the right way to include individual time-varying covs: http://www.phidot.org/forum/viewtopic.php?f=21&t=3775
# phi.sat_p.year.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + time:in.cr")
# phi.sat_p.year_sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + time:in.cr")
# phi.sat_p.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "sex + time:in.cr")
# phi.sat_p.year_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + time:in.cr")
# add creched - I think this is the right way to include individual time-varying covs: http://www.phidot.org/forum/viewtopic.php?f=21&t=3775
# pretty sure we should just end up with 1 beta for the time-varying indiv cov. using incr alone yields 1 beta, time:in.cr yields many betas
phi.sat_p.year.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + in.cr")
phi.sat_p.year_sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + in.cr")
phi.sat_p.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "sex + in.cr")
phi.sat_p.year_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + in.cr")

phi.sat_p.cr <- run.models(phi.stru = phi.sat, p.stru = "in.cr")
phi.sat_p.int <- run.models(phi.stru = phi.sat, p.stru = "1")




##  now collect models into one object to make a model comparison table
step1.1_p = collect.models(lx = c("phi.sat_p.year.sex", "phi.sat_p.year_sex", "phi.sat_p.sex", "phi.sat_p.year", 
                                  "phi.sat_p.year.sex_t", "phi.sat_p.year_sex_t", "phi.sat_p.sex_t", "phi.sat_p.year_t", 
                                  "phi.sat_p.year.sex_T", "phi.sat_p.year_sex_T", "phi.sat_p.sex_T", "phi.sat_p.year_T", 
                                  "phi.sat_p.year.sex_TT", "phi.sat_p.year_sex_TT", "phi.sat_p.sex_TT", "phi.sat_p.year_TT", 
                                  "phi.sat_p.year.sex_lnT", "phi.sat_p.year_sex_lnT", "phi.sat_p.sex_lnT", "phi.sat_p.year_lnT", 
                                  "phi.sat_p.year.sex_cr", "phi.sat_p.year_sex_cr", "phi.sat_p.sex_cr", "phi.sat_p.year_cr", 
                                  "phi.sat_p.cr", "phi.sat_p.int"))

# model.table(step1.1_p) %>% view()


step1.1_p$phi.sat_p.year_cr$results$beta %>% view()


saveRDS(step1.1_p, here("fitted_models/survival/step1.1_p"))


# step 1.2 Determine phi structure ----
###	now moving on to modeling Phi, with the best p structure from above, which is SEASON+time

# for my thesis and first round of publication draft, Phi structure began with 2 steps: first just SEASON and Sex structure (table 4.4 in thesis; table 3b in pub ms), then relative hatch date (tables 4.5 and 3c).
# in editing ms version we decided to combine these 2 steps for a cleaner and more readable analysis; this required fitting more models than the sum of the 2 original candidate sets.

p.sat = "SEASON * sex + time"

phi.year.sex_p.sat <- run.models(phi.stru = "SEASON * sex", p.stru = p.sat)
phi.year_sex_p.sat <- run.models(phi.stru = "SEASON + sex", p.stru = p.sat)
phi.sex_p.sat <- run.models(phi.stru = "sex", p.stru = p.sat)
phi.year_p.sat <- run.models(phi.stru = "SEASON", p.stru = p.sat)
# add t
phi.year.sex_t_p.sat <- run.models(phi.stru = "SEASON * sex + time", p.stru = p.sat)
phi.year_sex_t_p.sat <- run.models(phi.stru = "SEASON + sex + time", p.stru = p.sat)
phi.sex_t_p.sat <- run.models(phi.stru = "sex + time", p.stru = p.sat)
phi.year_t_p.sat <- run.models(phi.stru = "SEASON + time", p.stru = p.sat)
# add T
phi.year.sex_T_p.sat <- run.models(phi.stru = "SEASON * sex + Time", p.stru = p.sat)
phi.year_sex_T_p.sat <- run.models(phi.stru = "SEASON + sex + Time", p.stru = p.sat)
phi.sex_T_p.sat <- run.models(phi.stru = "sex + Time", p.stru = p.sat)
phi.year_T_p.sat <- run.models(phi.stru = "SEASON + Time", p.stru = p.sat)
# add TT
phi.year.sex_TT_p.sat <- run.models(phi.stru = "SEASON * sex + Time + I(Time^2)", p.stru = p.sat)
phi.year_sex_TT_p.sat <- run.models(phi.stru = "SEASON + sex + Time + I(Time^2)", p.stru = p.sat)
phi.sex_TT_p.sat <- run.models(phi.stru = "sex + Time + I(Time^2)", p.stru = p.sat)
phi.year_TT_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2)", p.stru = p.sat)
#add lnT
phi.year.sex_lnT_p.sat <- run.models(phi.stru = "SEASON * sex + log(Time+1)", p.stru = p.sat)
phi.year_sex_lnT_p.sat <- run.models(phi.stru = "SEASON + sex + log(Time+1)", p.stru = p.sat)
phi.sex_lnT_p.sat <- run.models(phi.stru = "sex + log(Time+1)", p.stru = p.sat)
phi.year_lnT_p.sat <- run.models(phi.stru = "SEASON + log(Time+1)", p.stru = p.sat)
# add hatch date
phi.year.sex_hatch_p.sat <- run.models(phi.stru = "SEASON * sex + res.htch", p.stru = p.sat)
phi.year_sex_hatch_p.sat <- run.models(phi.stru = "SEASON + sex + res.htch", p.stru = p.sat)
phi.sex_hatch_p.sat <- run.models(phi.stru = "sex + res.htch", p.stru = p.sat)
phi.year_hatch_p.sat <- run.models(phi.stru = "SEASON + res.htch", p.stru = p.sat)
# add t and hatch date
phi.year.sex_t_hatch_p.sat <- run.models(phi.stru = "SEASON * sex + time + res.htch", p.stru = p.sat)
phi.year_sex_t_hatch_p.sat <- run.models(phi.stru = "SEASON + sex + time + res.htch", p.stru = p.sat)
phi.sex_t_hatch_p.sat <- run.models(phi.stru = "sex + time + res.htch", p.stru = p.sat)
phi.year_t_hatch_p.sat <- run.models(phi.stru = "SEASON + time + res.htch", p.stru = p.sat)
# add T and hatch date
phi.year.sex_T_hatch_p.sat <- run.models(phi.stru = "SEASON * sex + Time + res.htch", p.stru = p.sat)
phi.year_sex_T_hatch_p.sat <- run.models(phi.stru = "SEASON + sex + Time + res.htch", p.stru = p.sat)
phi.sex_T_hatch_p.sat <- run.models(phi.stru = "sex + Time + res.htch", p.stru = p.sat)
phi.year_T_hatch_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch", p.stru = p.sat)
# add TT and hatch date
phi.year.sex_TT_hatch_p.sat <- run.models(phi.stru = "SEASON * sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.year_sex_TT_hatch_p.sat <- run.models(phi.stru = "SEASON + sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.sex_TT_hatch_p.sat <- run.models(phi.stru = "sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.year_TT_hatch_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch", p.stru = p.sat)
#add lnT and hatch date
phi.year.sex_lnT_hatch_p.sat <- run.models(phi.stru = "SEASON * sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.year_sex_lnT_hatch_p.sat <- run.models(phi.stru = "SEASON + sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.sex_lnT_hatch_p.sat <- run.models(phi.stru = "sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.year_lnT_hatch_p.sat <- run.models(phi.stru = "SEASON + log(Time+1) + res.htch", p.stru = p.sat)

phi.int_p.sat <- run.models(phi.stru = "1", p.stru = p.sat)



step1.2_phi =collect.models(lx=c("phi.year.sex_p.sat", "phi.year_sex_p.sat", "phi.sex_p.sat", "phi.year_p.sat", 
                              "phi.year.sex_t_p.sat", "phi.year_sex_t_p.sat", "phi.sex_t_p.sat", "phi.year_t_p.sat", 
                              "phi.year.sex_T_p.sat", "phi.year_sex_T_p.sat", "phi.sex_T_p.sat", "phi.year_T_p.sat", 
                              "phi.year.sex_TT_p.sat", "phi.year_sex_TT_p.sat", "phi.sex_TT_p.sat", "phi.year_TT_p.sat", 
                              "phi.year.sex_lnT_p.sat", "phi.year_sex_lnT_p.sat", "phi.sex_lnT_p.sat", "phi.year_lnT_p.sat", 
                              "phi.year.sex_hatch_p.sat", "phi.year_sex_hatch_p.sat", "phi.sex_hatch_p.sat", "phi.year_hatch_p.sat", 
                              "phi.year.sex_t_hatch_p.sat", "phi.year_sex_t_hatch_p.sat", "phi.sex_t_hatch_p.sat", "phi.year_t_hatch_p.sat", 
                              "phi.year.sex_T_hatch_p.sat", "phi.year_sex_T_hatch_p.sat", "phi.sex_T_hatch_p.sat", "phi.year_T_hatch_p.sat", 
                              "phi.year.sex_TT_hatch_p.sat", "phi.year_sex_TT_hatch_p.sat", "phi.sex_TT_hatch_p.sat", "phi.year_TT_hatch_p.sat", 
                              "phi.year.sex_lnT_hatch_p.sat", "phi.year_sex_lnT_hatch_p.sat", "phi.sex_lnT_hatch_p.sat", "phi.year_lnT_hatch_p.sat", 
                              "phi.int_p.sat"))

model.table(step1.2_phi)


saveRDS(step1.2_phi, here("fitted_models/survival/step1.2_phi"))

# phi_step1_2 <- readRDS(here("fitted_models/survival/phi_step1_models"))


#
# step 1.3: first check for uninformative parms where DQAICC < 5 in step 5.1  ----
# p ----
step1.1_p <- readRDS(here("fitted_models/survival/step1.1_p"))

p_mods <- step1.1_p[-length(step1.1_p)]

p_informative <- map2_df(p_mods, "p", cjs_parm_informative)

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

# p_informative has a row for each parameter in each model, but note it isn't meaningful to check each time varying parm for un vs. informative so these are excluded
p_informative <- p_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexMale", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         parm = sub("log\\(Time \\+ 1\\)", "lnT", parm),
         parm = sub("I\\(Time\\^2\\)", "TT", parm),
         parm = sub("Time", "T", parm)) %>%
  full_join(., step1.1_p$model.table %>% 
              data.frame() %>% 
              select(model.name = model, DeltaQAICc)) %>% 
  mutate(informative85 = ifelse(DeltaQAICc > 5, TRUE, informative85),
         parm = ifelse(grepl("~1", model.name), "p.intercept", parm),
         parm = ifelse(grepl("p\\(\\~time:in.cr\\)", model.name), "p.cr", parm))


         

p_informative_wide <- p_informative %>% 
           filter(!is.na(parm)) %>% 
  pivot_wider(id_cols = model.name, names_from = parm, values_from = informative85) %>% 
  mutate(across(c("p.SEASON", "p.sex", "p.SEASON.sex", "p.T", "p.TT", "p.lnT", "p.intercept"), ~replace_na(., TRUE))) %>% 
  mutate(p.stru = gsub("Phi\\(\\~SEASON \\* sex \\+ time\\)p\\(\\~", "", model.name),
         p.stru = gsub(":", ".", p.stru),
         p.stru = sub("resid.", "", p.stru),
         p.stru = sub("sexMale", "sex", p.stru),
         p.stru = sub("SEASON1314", "SEASON", p.stru),
         p.stru = sub("log\\(Time \\+ 1\\)", "lnT", p.stru),
         p.stru = sub("Time \\+ I\\(Time\\^2\\)", "TT", p.stru),
         p.stru = sub("Time", "T", p.stru),
         p.stru = sub("in.cr", paste("Cr", "\u00E8", "che", sep = ""), p.stru),
         p.stru = sub("time", "t", p.stru),
         p.stru = sub("\\)", "", p.stru),
         p.stru = ifelse(p.SEASON.sex == FALSE, str_replace(p.stru, "SEASON \\* sex", paste("(SEASON \\* sex)", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.SEASON == FALSE, str_replace(p.stru, "SEASON", paste("SEASON", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.sex == FALSE, str_replace(p.sex, "sex", paste("sex", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.T == FALSE, str_replace(p.T, "T", paste("T", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.TT == FALSE, str_replace(p.TT, "TT", paste("TT", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.lnT == FALSE, str_replace(p.lnT, "lnT", paste("lnT", "\u2020", sep = "")), p.stru),
         p.stru = ifelse(p.stru == "1", "Intercept only", p.stru),
         p.stru = mod_call_to_structure(p.stru)) %>% 
  select(model = model.name, p.stru)

step1.1_p$model.table <- step1.1_p$model.table %>% 
  data.frame() %>% 
  full_join(p_informative_wide)

step1.1_p$model.table %>% view()

saveRDS(step1.1_p, here("fitted_models/survival/step1.1_p"))


# step1.1_p$model.table %>% view()

# sex:SEASON parm is uninformative at P = 0.157 (85% CI) in phi.sat_p.year.sex_t (DAICc = 2.19)
# so just carrying phi.sat_p.year_sex_t and phi.sat_p.year_t forward

# phi ----
step1.2_phi <- readRDS(here("fitted_models/survival/step1.2_phi"))

phi_mods <- step1.2_phi[-length(step1.2_phi)]

phi_informative <- map2_df(phi_mods, "Phi", cjs_parm_informative)


phi_informative <- phi_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexMale", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         parm = sub("log\\(Time \\+ 1\\)", "lnT", parm),
         parm = sub("I\\(Time\\^2\\)", "TT", parm),
         parm = sub("Time", "T", parm)) %>%
  full_join(., step1.2_phi$model.table %>% 
              data.frame() %>% 
              select(model.name = model, DeltaQAICc))%>% 
  mutate(informative85 = ifelse(DeltaQAICc > 5, TRUE, informative85),
         parm = ifelse(grepl("~1", model.name), "Phi.intercept", parm))


         

phi_informative_wide <- phi_informative %>% 
           filter(!is.na(parm)) %>% 
  pivot_wider(id_cols = model.name, names_from = parm, values_from = informative85) %>% 
  mutate(across(c("Phi.SEASON", "Phi.sex", "Phi.SEASON.sex", "Phi.T", "Phi.TT", "Phi.lnT", "Phi.intercept"), ~replace_na(., TRUE))) %>% 
  mutate(phi.stru = gsub("\\)p\\(\\~SEASON \\* sex \\+ time\\)", "", model.name),
         phi.stru = gsub(":", ".", phi.stru),
         phi.stru = sub("resid.", "", phi.stru),
         phi.stru = sub("sexMale", "sex", phi.stru),
         phi.stru = sub("SEASON1314", "SEASON", phi.stru),
         phi.stru = sub("log\\(Time \\+ 1\\)", "lnT", phi.stru),
         phi.stru = sub("Time \\+ I\\(Time\\^2\\)", "TT", phi.stru),
         phi.stru = sub("Time", "T", phi.stru),
         phi.stru = sub("time.in.cr", "Crèche", phi.stru),
         phi.stru = sub("time", "t", phi.stru),
         phi.stru = sub("res.htch", "Hatch date", phi.stru),
         phi.stru = sub("Phi\\(\\~", "", phi.stru),
         phi.stru = ifelse(Phi.SEASON.sex == FALSE, str_replace(phi.stru, "SEASON \\* sex", paste("(SEASON \\* sex)", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.SEASON == FALSE, str_replace(phi.stru, "SEASON", paste("SEASON", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.sex == FALSE, str_replace(phi.stru, "sex", paste("sex", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.T == FALSE & !grepl("TT", phi.stru), str_replace(phi.stru, "T", paste("T", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.TT == FALSE, str_replace(phi.stru, "TT", paste("TT", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.lnT == FALSE, str_replace(phi.stru, "lnT", paste("lnT", "\u2020", sep = "")), phi.stru),
         phi.stru = ifelse(phi.stru == "1", "Intercept only", phi.stru),
         phi.stru = mod_call_to_structure(phi.stru)) %>% 
  select(model = model.name, phi.stru)

step1.2_phi$model.table <- step1.2_phi$model.table %>% 
  data.frame() %>% 
  full_join(phi_informative_wide)

step1.2_phi$model.table %>% view()

saveRDS(step1.2_phi, here("fitted_models/survival/step1.2_phi"))


# in TT models, T2 parm is informative but T not; keeping TT models
# lnT always uninformative so dropping those mods
# sex always uninformative so dropping those mods
# keeping:
# "phi.year_TT_hatch_p.sat"
# "phi.year_T_hatch_p.sat"
# "phi.year_hatch_p.sat"


# Step 1.4: create new candidate set from the best structures from above ----

step1.1_p <- readRDS(here("fitted_models/survival/step1.1_p"))
step1.2_phi <- readRDS(here("fitted_models/survival/step1.2_phi"))

phi.sat_p.year_sex_t <- step1.1_p$phi.sat_p.year_sex_t
phi.sat_p.year_t <- step1.1_p$phi.sat_p.year_t
phi.year_TT_hatch_p.sat <- step1.2_phi$phi.year_TT_hatch_p.sat
phi.year_T_hatch_p.sat <- step1.2_phi$phi.year_T_hatch_p.sat
phi.year_hatch_p.sat <- step1.2_phi$phi.year_hatch_p.sat

step1.4 = collect.models(lx=c("phi.sat_p.year_sex_t", "phi.sat_p.year_t", 
                              "phi.year_TT_hatch_p.sat", "phi.year_T_hatch_p.sat", "phi.year_hatch_p.sat"))


step1.4$model.table <- step1.4$model.table %>% 
  mutate(model.structure = model,
         model.structure = gsub(":", ".", model.structure),
         model.structure = gsub("res.htch", "Hatch date", model.structure),
         model.structure = gsub("sex", "Sex", model.structure),
         model.structure = gsub("SEASON", "Yr", model.structure),
         model.structure = gsub("log\\(Time \\+ 1\\)", "lnT", model.structure),
         model.structure = gsub("Time \\+ I\\(Time\\^2\\)", "TT", model.structure),
         model.structure = gsub("Time", "T", model.structure),
         model.structure = gsub("in.cr", paste("Cr", "\u00E8", "che", sep = ""), model.structure),
         model.structure = gsub("time", "t", model.structure),
         model.structure = gsub("\\)p\\(", "\\) p\\(", model.structure))


step1.4$model.table %>% view()
saveRDS(step1.4, here("fitted_models/survival/step1.4"))

#
# step 2. main variables of interest for overall survival: Age and creched status ----


step1.4mods <- readRDS(here("fitted_models/survival/step1.4"))

#
# more code to generate code NO RUN ----
# get DQAICc <= 5 models from previous step
step2_cand_set_base <- step1.4mods %>% 
  model.table() %>% 
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(DeltaQAICc <= 5) %>% 
  left_join(., names(step1.4mods) %>%
              data.frame() %>%
              rownames_to_column("mod.num") %>%
              rename(mod.name = 2)) %>% 
  select(mod.name, Phi, p, DeltaQAICc) %>% 
  mutate(Phi = gsub("~", "", Phi),
         p = gsub("~", "", p))

# best from step 1.3
step2_cand_set_base %>% 
  mutate(mod.assign = paste(mod.name, " <- step1.3mods$", mod.name, sep = "")) %>% 
  select(mod.assign)
# now adding in time varying variables. these always need to have time, either time+ or time:
# best from step 1.3 plus additive time and time varying creched.	#same intercept, diff slope for in.cr
step2_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_t_incr_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + time + in.cr\")", sep = "")) %>% 
  select(mod.assign) 
# best from step 1.3 plus additive time and time varying age. same intercept, diff slope for dayold
step2_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_t_age_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + time + dayold\")", sep = "")) %>% 
  select(mod.assign)
# best from step 1.3 plus interaction of time and time varying creched.	#same slope, diff intercept for in.cr
step2_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_t.incr_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + time:in.cr\")", sep = "")) %>% 
  select(mod.assign)
# best from step 1.3 plus additive time and time varying age. same slope, diff intercept for dayold
step2_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_t.age_p.", mod.name), " <- run.models(phi.stru = \"", Phi, " + time:dayold\")", sep = "")) %>% 
  select(mod.assign)

# make code for collect.models call
rbind(step2_cand_set_base %>% 
  select(mod.name),
# now adding in time varying variables. these always need to have time, either time+ or time:
# best from step 1.3 plus additive time and time varying creched.	#same intercept, diff slope for in.cr
data.frame(mod.name = "#t_incr"),
step2_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_t_incr_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus additive time and time varying age. same intercept, diff slope for dayold
data.frame(mod.name = "#t_age"),
step2_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_t_age_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus interaction of time and time varying creched.	#same slope, diff intercept for in.cr
data.frame(mod.name = "#t:incr"),
step2_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_t.incr_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus additive time and time varying age. same slope, diff intercept for dayold
data.frame(mod.name = "#t:age"),
step2_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_t.age_p.", mod.name)) %>% 
  select(mod.name)) %>% 
  summarise(collect.call = paste(mod.name, collapse = "\", \""))


# now assign/run all the models ----
# best from 1.3, don't need to refit these
phi.year_TT_hatch_p.sat <- step1.4mods$phi.year_TT_hatch_p.sat
phi.year_T_hatch_p.sat <- step1.4mods$phi.year_T_hatch_p.sat
phi.year_hatch_p.sat <- step1.4mods$phi.year_hatch_p.sat
# plus time-varying individual covariate for creched status - this is the proper way to code the time-varying indiv covs
phi.year_TT_hatch_incr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + in.cr")
phi.year_T_hatch_incr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + in.cr")
phi.year_hatch_incr_p.sat <- run.models(phi.stru = "SEASON + res.htch + in.cr")
# plus time varying age - this is the proper way to code the time-varying indiv covs
phi.year_TT_hatch_age_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + dayold")
phi.year_T_hatch_age_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + dayold")
phi.year_hatch_age_p.sat <- run.models(phi.stru = "SEASON + res.htch + dayold")


# export some DMs to double check
phi.year_T_hatch_incr_p.sat$design.matrix %>% write.csv(here("fitted_models/phi.year_TT_hatch_t_incr_p.sat.csv"))
phi.year_T_hatch.age_p.sat$design.matrix %>% write.csv(here("fitted_models/phi.year_TT_hatch_t.incr_p.sat.csv"))


step2 <- collect.models(c("phi.year_TT_hatch_p.sat", "phi.year_T_hatch_p.sat", "phi.year_hatch_p.sat", 
                          #incr
                          "phi.year_TT_hatch_incr_p.sat", "phi.year_T_hatch_incr_p.sat", "phi.year_hatch_incr_p.sat", 
                          #age
                          "phi.year_TT_hatch_age_p.sat", "phi.year_T_hatch_age_p.sat", "phi.year_hatch_age_p.sat"))


step2 <- readRDS(here("fitted_models/survival/step2"))

step2$model.table <- step2$model.table %>%
  data.frame() %>%  
  mutate(model.structure = model,
         model.structure = gsub(":", ".", model.structure),
         model.structure = gsub("res.htch", "Hatch date", model.structure),
         model.structure = gsub("sex", "Sex", model.structure),
         model.structure = gsub("SEASON", "Yr", model.structure),
         model.structure = gsub("log\\(Time \\+ 1\\)", "lnT", model.structure),
         model.structure = gsub("Time \\+ I\\(Time\\^2\\)", "TT", model.structure),
         model.structure = gsub("Time", "T", model.structure),
         model.structure = gsub("in.cr", paste("Cr", "\u00E8", "che", sep = ""), model.structure),
         model.structure = gsub("time\\.", "t \\* ", model.structure),
         model.structure = gsub("time", "t", model.structure),
         model.structure = gsub("dayold", "Age", model.structure),
         model.structure = gsub("\\)p\\(", "\\) p\\(", model.structure))


step2$model.table %>% view()

saveRDS(step2, here("fitted_models/survival/step2"))
# model averaging ----




step2 <- readRDS(here("fitted_models/survival/step2"))

step2_mod_av_phi = model.average(step2, "Phi", vcv=TRUE, drop = FALSE)

saveRDS(step2_mod_av_phi, here("fitted_models/survival/step2_mod_av_phi"))

