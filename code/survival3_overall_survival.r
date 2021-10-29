


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
# add creched
phi.sat_p.year.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON * sex + time:in.cr")
phi.sat_p.year_sex_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + sex + time:in.cr")
phi.sat_p.sex_cr <- run.models(phi.stru = phi.sat, p.stru = "sex + time:in.cr")
phi.sat_p.year_cr <- run.models(phi.stru = phi.sat, p.stru = "SEASON + time:in.cr")

phi.sat_p.cr <- run.models(phi.stru = phi.sat, p.stru = "time:in.cr")
phi.sat_p.int <- run.models(phi.stru = phi.sat, p.stru = "1")




##  now collect models into one object to make a model comparison table
step1.1_p = collect.models(lx = c("phi.sat_p.year.sex", "phi.sat_p.year_sex", "phi.sat_p.sex", "phi.sat_p.year", 
                                  "phi.sat_p.year.sex_t", "phi.sat_p.year_sex_t", "phi.sat_p.sex_t", "phi.sat_p.year_t", 
                                  "phi.sat_p.year.sex_T", "phi.sat_p.year_sex_T", "phi.sat_p.sex_T", "phi.sat_p.year_T", 
                                  "phi.sat_p.year.sex_TT", "phi.sat_p.year_sex_TT", "phi.sat_p.sex_TT", "phi.sat_p.year_TT", 
                                  "phi.sat_p.year.sex_lnT", "phi.sat_p.year_sex_lnT", "phi.sat_p.sex_lnT", "phi.sat_p.year_lnT", 
                                  "phi.sat_p.year.sex_cr", "phi.sat_p.year_sex_cr", "phi.sat_p.sex_cr", "phi.sat_p.year_cr", 
                                  "phi.sat_p.cr", "phi.sat_p.int"))

# model.table(step1.1_p)

step1.1_p$phi.sat_p.year_cr$results$beta %>% view()


saveRDS(step1.1_p, here("fitted_models/survival/step1.1_p"))

readRDS(here("fitted_models/survival/step1.1_p"))$model.table
# peng.p2 <- readRDS(here("fitted_models/survival/p_models"))

readRDS(here("fitted_models/survival/step1.1_p"))$phi.sat_p.year_sex_cr$results$beta %>% view()
	
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
# Step 1.3: create new candidate set from the best structures from above ----
# step 1.3.1: first check for uninformative parms where DQAICC < 5 in step 5.1  ----

step1.1_p <- readRDS(here("fitted_models/survival/step1.1_p"))

step1.1_p %>% 
  model.table(model.name = FALSE) %>% 
  filter(DeltaQAICc <= 5) %>% 
  mutate(uninform.assign = paste("cjs_parm_informative(step1.1_p$", model, ", \"p\")", sep = "")) %>% 
  summarise(uninform.assign = paste(uninform.assign, collapse = ", ")) %>% 
  mutate(uninform.assign = paste("p_informative <- rbind(", uninform.assign, ")", sep = ""))

# for some reason map2_df not working for this so need to manually rbind
p_informative <- rbind(cjs_parm_informative(step1.1_p$phi.sat_p.year_sex_t, "p"), 
                             cjs_parm_informative(step1.1_p$phi.sat_p.year_t, "p"), 
                             cjs_parm_informative(step1.1_p$phi.sat_p.year.sex_t, "p"))

p_informative <- p_informative %>% 
  mutate(mod.name = gsub("step1.1_p\\$", "", mod.name))

p_informative_wide <- p_informative %>% 
  pivot_wider(id_cols = mod.name, names_from = parm, values_from = informative85) %>% 
  mutate(uninformative = "")


p_informative_wide <- edit(p_informative_wide)

p_informative_wide %>%
  filter(uninformative == TRUE) %>%  
  select(mod.name, uninformative) %>% 
  saveRDS(here("fitted_models/survival/step1.1_uninformative"))


# sex:SEASON parm is uninformative at P = 0.157 (85% CI) in phi.sat_p.year.sex_t (DAICc = 2.19)
# so just carrying phi.sat_p.year_sex_t and phi.sat_p.year_t forward


step1.2_phi <- readRDS(here("fitted_models/survival/step1.2_phi"))

step1.2_phi %>% 
  model.table(model.name = FALSE) %>% 
  filter(DeltaQAICc <= 5) %>% 
  mutate(uninform.assign = paste("cjs_parm_informative(step1.2_phi$", model, ", \"Phi\")", sep = "")) %>% 
  summarise(uninform.assign = paste(uninform.assign, collapse = ", ")) %>% 
  mutate(uninform.assign = paste("phi_informative <- rbind(", uninform.assign, ")", sep = ""))



phi_informative <- rbind(cjs_parm_informative(step1.2_phi$phi.year_TT_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_T_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year.sex_TT_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_lnT_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year.sex_T_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_sex_TT_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year.sex_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_sex_T_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year.sex_lnT_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_sex_hatch_p.sat, "Phi"), 
                             cjs_parm_informative(step1.2_phi$phi.year_sex_lnT_hatch_p.sat, "Phi"))


phi_informative <- phi_informative %>% 
  mutate(mod.name = gsub("step1.2_phi\\$", "", mod.name))

phi_informative_wide <- phi_informative %>% 
  pivot_wider(id_cols = mod.name, names_from = parm, values_from = informative85) %>% 
  mutate(uninformative = "")

phi_informative_wide <- edit(phi_informative_wide)

phi_informative_wide %>%
  filter(uninformative == TRUE) %>%  
  select(mod.name, uninformative) %>% 
  saveRDS(here("fitted_models/survival/step1.2_uninformative"))


# in TT models, T2 parm is informative but T not; keeping TT models
# lnT always uninformative so dropping those mods
# sex always uninformative so dropping those mods
# keeping:
# "phi.year_TT_hatch_p.sat"
# "phi.year_T_hatch_p.sat"
# "phi.year_hatch_p.sat"



step5.2_uninformative <- phi_informative_wide %>%
  select(mod.name, uninformative) %>% 
  filter(uninformative == TRUE) %>% 
  left_join(., phi_best) %>% 
  select(model, mod.name, uninformative)

saveRDS(step5.2_uninformative, here("fitted_models/survival/step5.2_uninformative"))



step1.3 = collect.models(lx=c("phi.sat_p.year_sex_t", "phi.sat_p.year_t", 
                              # "phi.sat_p.year.sex_t", uninformative
                              "phi.year_TT_hatch_p.sat", "phi.year_T_hatch_p.sat", "phi.year_hatch_p.sat"
                              # "phi.year.sex_TT_hatch_p.sat", "phi.year_lnT_hatch_p.sat", "phi.year.sex_T_hatch_p.sat", "phi.year_sex_TT_hatch_p.sat", "phi.year.sex_hatch_p.sat", "phi.year_sex_T_hatch_p.sat", "phi.year.sex_lnT_hatch_p.sat", "phi.year_sex_hatch_p.sat", "phi.year_sex_lnT_hatch_p.sat" these uninformative
                              ))


model.table(step1.3) %>% 
  select(-model)
saveRDS(step1.3, here("fitted_models/survival/step1.3"))


# step 2. main variables of interest for overall survival: Age and creched status ----


step1.3mods <- readRDS(here("fitted_models/survival/step1.3"))

step2_cand_set_base <- step1.3mods %>% 
  model.table() %>% 
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(DeltaQAICc <= 5) %>% 
  left_join(., names(step1.3mods) %>%
              data.frame() %>%
              rownames_to_column("mod.num") %>%
              rename(mod.name = 2)) %>% 
  select(mod.name, Phi, DeltaQAICc) %>% 
  mutate(Phi = gsub("~", "", Phi))

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


# assign/run all the models
# best from 1.3, don't need to refit these
phi.year_TT_hatch_p.sat <- step1.3mods$phi.year_TT_hatch_p.sat
phi.year_T_hatch_p.sat <- step1.3mods$phi.year_T_hatch_p.sat
phi.year_hatch_p.sat <- step1.3mods$phi.year_hatch_p.sat
# plus additive time and creched
phi.year_TT_hatch_t_incr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + time + in.cr")
phi.year_T_hatch_t_incr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + time + in.cr")
phi.year_hatch_t_incr_p.sat <- run.models(phi.stru = "SEASON + res.htch + time + in.cr")
# plus additive time and age
phi.year_TT_hatch_t_age_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + time + dayold")
phi.year_T_hatch_t_age_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + time + dayold")
phi.year_hatch_t_age_p.sat <- run.models(phi.stru = "SEASON + res.htch + time + dayold")
# plus interaction time and creched
phi.year_TT_hatch_t.incr_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + time:in.cr")
phi.year_T_hatch_t.incr_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + time:in.cr")
phi.year_hatch_t.incr_p.sat <- run.models(phi.stru = "SEASON + res.htch + time:in.cr")
# plus interaction time and age
phi.year_TT_hatch_t.age_p.sat <- run.models(phi.stru = "SEASON + Time + I(Time^2) + res.htch + time:dayold")
phi.year_T_hatch_t.age_p.sat <- run.models(phi.stru = "SEASON + Time + res.htch + time:dayold")
phi.year_hatch_t.age_p.sat <- run.models(phi.stru = "SEASON + res.htch + time:dayold")


# export some DMs to double check
phi.year_TT_hatch_t_incr_p.sat$design.matrix %>% write.csv(here("fitted_models/phi.year_TT_hatch_t_incr_p.sat.csv"))
phi.year_TT_hatch_t.incr_p.sat$design.matrix %>% write.csv(here("fitted_models/phi.year_TT_hatch_t.incr_p.sat.csv"))


step2 <- collect.models(c("phi.year_TT_hatch_p.sat", "phi.year_T_hatch_p.sat", "phi.year_hatch_p.sat", 
                          #t_incr
                          "phi.year_TT_hatch_t_incr_p.sat", "phi.year_T_hatch_t_incr_p.sat", "phi.year_hatch_t_incr_p.sat", 
                          #t_age
                          "phi.year_TT_hatch_t_age_p.sat", "phi.year_T_hatch_t_age_p.sat", "phi.year_hatch_t_age_p.sat", 
                          #t:incr
                          "phi.year_TT_hatch_t.incr_p.sat", "phi.year_T_hatch_t.incr_p.sat", "phi.year_hatch_t.incr_p.sat", 
                          #t:age
                          "phi.year_TT_hatch_t.age_p.sat", "phi.year_T_hatch_t.age_p.sat", "phi.year_hatch_t.age_p.sat"))

step2$model.table %>%
  data.frame() %>%  
  rownames_to_column("mod.num") %>% 
  left_join(names(step2) %>%
              data.frame() %>%
              rename(mod.name = 1) %>%
              rownames_to_column("mod.num")) %>% view()

saveRDS(step2, here("fitted_models/survival/step2"))
# model averaging ----




step2 <- readRDS(here("fitted_models/survival/step2"))

step2_mod_av_phi = model.average(step2, "Phi", vcv=TRUE, drop = FALSE)

saveRDS(step2_mod_av_phi, here("fitted_models/survival/step2_mod_av_phi"))

