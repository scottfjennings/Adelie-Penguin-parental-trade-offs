




#ACR computer
MarkPath="C:/Program Files (x86)/MARK"
library(RMark)
library(AICcmodavg)
library(tidyverse)
library(R2ucare)
library(here)
options(scipen = 999)

source(here("code/utilities.R"))

# setting working directory only to control location for all the Mark output files
setwd(here("mark_output"))


# for this final part of the analysis we're considering just chicks that reached the creche stage. so, go back and do the first stages of model selection (base p and phi) with this subset (mirroring steps done in survival3_...R)


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
penguins=convert.inp(here("data/mark_in_cr.inp"), 
					group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), 
					covariates=c("res.htch", "cr.age", "cr.mass", "cr.flip", "cr.tib"), use.comments = TRUE) 




penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 60) 
penguins.ddl = make.design.data(penguins.process)

# best model for overall survival is Phi(~SEASON + res.htch + Time + I(Time^2))p(~SEASON + time)

# define helper functions ----

# generic function for running models
run_penguin_models<-function(phi.stru, p.stru = "SEASON * sex + time") {
  phi.stru.list <- list(formula = formula(paste("~", phi.stru)))
  p.stru.list <- list(formula = formula(paste("~", p.stru)))
zmod <- mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru.list, p = p.stru.list), output = FALSE, chat = 1.16)
return(zmod)
}


# step 5.1 Determine p structure ----

phi.sat = "SEASON * sex + time"

phi.sat_p.year.sex <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON * sex")
phi.sat_p.year_sex <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + sex")
phi.sat_p.sex <- run_penguin_models(phi.stru = phi.sat, p.stru = "sex")
phi.sat_p.year <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON")
# add t
phi.sat_p.year.sex_t <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON * sex + time")
phi.sat_p.year_sex_t <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + sex + time")
phi.sat_p.sex_t <- run_penguin_models(phi.stru = phi.sat, p.stru = "sex + time")
phi.sat_p.year_t <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + time")
# add T
phi.sat_p.year.sex_T <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON * sex + Time")
phi.sat_p.year_sex_T <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + sex + Time")
phi.sat_p.sex_T <- run_penguin_models(phi.stru = phi.sat, p.stru = "sex + Time")
phi.sat_p.year_T <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + Time")
# add TT
phi.sat_p.year.sex_TT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON * sex + Time + I(Time^2)")
phi.sat_p.year_sex_TT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + sex + Time + I(Time^2)")
phi.sat_p.sex_TT <- run_penguin_models(phi.stru = phi.sat, p.stru = "sex + Time + I(Time^2)")
phi.sat_p.year_TT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + Time + I(Time^2)")
#add lnT
phi.sat_p.year.sex_lnT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON * sex + log(Time+1)")
phi.sat_p.year_sex_lnT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + sex + log(Time+1)")
phi.sat_p.sex_lnT <- run_penguin_models(phi.stru = phi.sat, p.stru = "sex + log(Time+1)")
phi.sat_p.year_lnT <- run_penguin_models(phi.stru = phi.sat, p.stru = "SEASON + log(Time+1)")

phi.sat_p.int <- run_penguin_models(phi.stru = phi.sat, p.stru = "1")




##  now collect models into one object to make a model comparison table
step5.1_p = collect.models(lx = c("phi.sat_p.year.sex", "phi.sat_p.year_sex", "phi.sat_p.sex", "phi.sat_p.year", 
                                  "phi.sat_p.year.sex_t", "phi.sat_p.year_sex_t", "phi.sat_p.sex_t", "phi.sat_p.year_t", 
                                  "phi.sat_p.year.sex_T", "phi.sat_p.year_sex_T", "phi.sat_p.sex_T", "phi.sat_p.year_T", 
                                  "phi.sat_p.year.sex_TT", "phi.sat_p.year_sex_TT", "phi.sat_p.sex_TT", "phi.sat_p.year_TT", 
                                  "phi.sat_p.year.sex_lnT", "phi.sat_p.year_sex_lnT", "phi.sat_p.sex_lnT", "phi.sat_p.year_lnT", 
                                  "phi.sat_p.int"))

# model.table(step5.1_p)

saveRDS(step5.1_p, here("fitted_models/survival/step5.1_creched_p"))

# peng.p2 <- readRDS(here("fitted_models/survival/p_models"))


	
# step 5.2 Determine phi structure ----


p.sat = "SEASON * sex + time"

phi.year.sex_p.sat <- run_penguin_models(phi.stru = "SEASON * sex", p.stru = p.sat)
phi.year_sex_p.sat <- run_penguin_models(phi.stru = "SEASON + sex", p.stru = p.sat)
phi.sex_p.sat <- run_penguin_models(phi.stru = "sex", p.stru = p.sat)
phi.year_p.sat <- run_penguin_models(phi.stru = "SEASON", p.stru = p.sat)
# add t
phi.year.sex_t_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + time", p.stru = p.sat)
phi.year_sex_t_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + time", p.stru = p.sat)
phi.sex_t_p.sat <- run_penguin_models(phi.stru = "sex + time", p.stru = p.sat)
phi.year_t_p.sat <- run_penguin_models(phi.stru = "SEASON + time", p.stru = p.sat)
# add T
phi.year.sex_T_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + Time", p.stru = p.sat)
phi.year_sex_T_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + Time", p.stru = p.sat)
phi.sex_T_p.sat <- run_penguin_models(phi.stru = "sex + Time", p.stru = p.sat)
phi.year_T_p.sat <- run_penguin_models(phi.stru = "SEASON + Time", p.stru = p.sat)
# add TT
phi.year.sex_TT_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + Time + I(Time^2)", p.stru = p.sat)
phi.year_sex_TT_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + Time + I(Time^2)", p.stru = p.sat)
phi.sex_TT_p.sat <- run_penguin_models(phi.stru = "sex + Time + I(Time^2)", p.stru = p.sat)
phi.year_TT_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + I(Time^2)", p.stru = p.sat)
#add lnT
phi.year.sex_lnT_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + log(Time+1)", p.stru = p.sat)
phi.year_sex_lnT_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + log(Time+1)", p.stru = p.sat)
phi.sex_lnT_p.sat <- run_penguin_models(phi.stru = "sex + log(Time+1)", p.stru = p.sat)
phi.year_lnT_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time+1)", p.stru = p.sat)
# add hatch date
phi.year.sex_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + res.htch", p.stru = p.sat)
phi.year_sex_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + res.htch", p.stru = p.sat)
phi.sex_hatch_p.sat <- run_penguin_models(phi.stru = "sex + res.htch", p.stru = p.sat)
phi.year_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + res.htch", p.stru = p.sat)
# add t and hatch date
phi.year.sex_t_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + time + res.htch", p.stru = p.sat)
phi.year_sex_t_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + time + res.htch", p.stru = p.sat)
phi.sex_t_hatch_p.sat <- run_penguin_models(phi.stru = "sex + time + res.htch", p.stru = p.sat)
phi.year_t_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + time + res.htch", p.stru = p.sat)
# add T and hatch date
phi.year.sex_T_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + Time + res.htch", p.stru = p.sat)
phi.year_sex_T_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + Time + res.htch", p.stru = p.sat)
phi.sex_T_hatch_p.sat <- run_penguin_models(phi.stru = "sex + Time + res.htch", p.stru = p.sat)
phi.year_T_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch", p.stru = p.sat)
# add TT and hatch date
phi.year.sex_TT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.year_sex_TT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.sex_TT_hatch_p.sat <- run_penguin_models(phi.stru = "sex + Time + I(Time^2) + res.htch", p.stru = p.sat)
phi.year_TT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + I(Time^2) + res.htch", p.stru = p.sat)
#add lnT and hatch date
phi.year.sex_lnT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON * sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.year_sex_lnT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.sex_lnT_hatch_p.sat <- run_penguin_models(phi.stru = "sex + log(Time+1) + res.htch", p.stru = p.sat)
phi.year_lnT_hatch_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time+1) + res.htch", p.stru = p.sat)

phi.int_p.sat <- run_penguin_models(phi.stru = "1", p.stru = p.sat)



step5.2_phi =collect.models(lx=c("phi.year.sex_p.sat", "phi.year_sex_p.sat", "phi.sex_p.sat", "phi.year_p.sat", 
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

model.table(step5.2_phi) %>% view()


saveRDS(step5.2_phi, here("fitted_models/survival/step5.2_creched_phi"))

# phi_step5_2 <- readRDS(here("fitted_models/survival/phi_step5_models"))


#
# step 5.3.1: first check for uninformative parms where DQAICC < 5 in step 5.1  ----
# p ----
step5.1_p <- readRDS(here("fitted_models/survival/step5.1_creched_p"))

p_mods <- step5.1_p[-length(step5.1_p)]

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
  full_join(., step5.1_p$model.table %>% 
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
         p.stru = ifelse(p.SEASON.sex == FALSE, str_replace(p.stru, "SEASON \\* sex", sprintf(paste("(SEASON \\* sex)", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.SEASON == FALSE, str_replace(p.stru, "SEASON", sprintf(paste("SEASON", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.sex == FALSE, str_replace(p.stru, "sex", sprintf(paste("sex", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.T == FALSE, str_replace(p.stru, "T", sprintf(paste("T", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.TT == FALSE, str_replace(p.stru, "TT", sprintf(paste("TT", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.lnT == FALSE, str_replace(p.stru, "lnT", sprintf(paste("lnT", "^\u2020^", sep = ""))), p.stru),
         p.stru = ifelse(p.stru == "1", "Intercept only", p.stru),
         p.stru = mod_call_to_structure(p.stru)) %>% 
  select(model = model.name, p.stru)


# sprintf(paste("lnT", "^\u2020^", sep = ""))

step5.1_p$model.table <- step5.1_p$model.table %>% 
  data.frame() %>% 
  full_join(p_informative_wide, by = c("model"))

step5.1_p$model.table %>% view()

saveRDS(step5.1_p, here("fitted_models/survival/step5.1_creched_p"))


# step5.1_p$model.table %>% view()

# sex and sex:SEASON parms are uninformative at P = 0.157 (85% CI) in phi.sat_p.year_sex_t (DAICc = 2.29) and phi.sat_p.year.sex_t (DAICc = 3.78)
# so just carrying phi.sat_p.year_t forward

# phi ----
step5.2_phi <- readRDS(here("fitted_models/survival/step5.2_creched_phi"))

phi_mods <- step5.2_phi[-length(step5.2_phi)]

phi_informative <- map2_df(phi_mods, "Phi", cjs_parm_informative)


phi_informative <- phi_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexMale", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         parm = sub("log\\(Time \\+ 1\\)", "lnT", parm),
         parm = sub("I\\(Time\\^2\\)", "TT", parm),
         parm = sub("Time", "T", parm)) %>%
  full_join(., step5.2_phi$model.table %>% 
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
         phi.stru = ifelse(Phi.SEASON.sex == FALSE, str_replace(phi.stru, "SEASON \\* sex", sprintf(paste("(SEASON \\* sex)", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.SEASON == FALSE, str_replace(phi.stru, "SEASON", sprintf(paste("SEASON", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.sex == FALSE, str_replace(phi.stru, "sex", sprintf(paste("sex", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.T == FALSE & !grepl("TT", phi.stru), str_replace(phi.stru, "T", sprintf(paste("T", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.TT == FALSE, str_replace(phi.stru, "TT", sprintf(paste("TT", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.lnT == FALSE, str_replace(phi.stru, "lnT", sprintf(paste("lnT", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(phi.stru == "1", "Intercept only", phi.stru),
         phi.stru = mod_call_to_structure(phi.stru)) %>% 
  select(model = model.name, phi.stru)

step5.2_phi$model.table <- step5.2_phi$model.table %>% 
  data.frame() %>% 
  full_join(phi_informative_wide, by = c("model"))

step5.2_phi$model.table %>% view()

saveRDS(step5.2_phi, here("fitted_models/survival/step5.2_creched_phi"))


# in TT models, T and T2 parms uninformative so dropping TT models
# sex and sex:SEASON always uninformative so dropping those mods
# keeping:
# "phi.year_lnT_hatch_p.sat"
# "phi.year_T_hatch_p.sat"
# "phi.year_lnT_hatch_p.sat"
# "phi.year_hatch_p.sat"



# Step 5.4: create new candidate set from the best structures from above ----
# don't need to refit, just get the appropriate models from step5.1_creched_p and step5.2_creched_phi

step5.1_p <- readRDS(here("fitted_models/survival/step5.1_creched_p"))
step5.2_phi <- readRDS(here("fitted_models/survival/step5.2_creched_phi"))


#code-generating code
rbind(p_informative_wide %>% filter(uninformative != TRUE) %>% select(mod.name) %>% mutate(mod.name = paste(mod.name, " <- step5.1_p$", mod.name, sep = "")),
      phi_informative_wide %>% filter(uninformative != TRUE) %>% select(mod.name) %>% mutate(mod.name = paste(mod.name, " <- step5.2_phi$", mod.name, sep = "")))
#
rbind(p_informative_wide %>% filter(uninformative != TRUE) %>% select(mod.name),
      phi_informative_wide %>% filter(uninformative != TRUE) %>% select(mod.name)) %>% 
  summarise(collect.call = paste(mod.name, collapse = "\", \"")) %>% 
  mutate(collect.call = paste("step5.3 = collect.models(lx=c(\"", collect.call, "\"))", sep = ""))


phi.sat_p.year_t <- step5.1_p$phi.sat_p.year_t
phi.year_T_hatch_p.sat <- step5.2_phi$phi.year_T_hatch_p.sat
phi.year_lnT_hatch_p.sat <- step5.2_phi$phi.year_lnT_hatch_p.sat
phi.year_T_p.sat <- step5.2_phi$phi.year_T_p.sat
phi.year_lnT_p.sat <- step5.2_phi$phi.year_lnT_p.sat

step5.4 = collect.models(lx=c("phi.sat_p.year_t", "phi.year_T_hatch_p.sat", "phi.year_lnT_hatch_p.sat", "phi.year_T_p.sat", "phi.year_lnT_p.sat"))

step5.4$model.table <- step5.4$model.table %>% 
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

step5.4$model.table %>% 
  view()

saveRDS(step5.4, here("fitted_models/survival/step5.4_creched_base"))
# forgot to change Phi in model.name
step5.4 <- readRDS(here("fitted_models/survival/step5.4_creched_base"))
step5.4$model.table <- step5.4$model.table %>% 
  mutate(model.structure = gsub("Phi", "\u03D5", model.structure))
saveRDS(step5.4, here("fitted_models/survival/step5.4_creched_base"))
# step 6 add creching age and size ----
# helper code to build up candidata set code ----

step5.4mods <- readRDS(here("fitted_models/survival/step5.4_creched_base"))

step6_cand_set_base <- step5.4mods %>% 
  model.table() %>% 
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(DeltaQAICc <= 5) %>% 
  left_join(., names(step5.4mods) %>%
              data.frame() %>%
              rownames_to_column("mod.num") %>%
              rename(mod.name = 2)) %>% 
  select(mod.name, Phi, DeltaQAICc) %>% 
  mutate(Phi = gsub("~", "", Phi))

# code-generating code ----
# code for fitting models ----
# best from step 1.3
step6_cand_set_base %>% 
  mutate(mod.assign = paste(mod.name, " <- step5.3mods$", mod.name, sep = "")) %>% 
  select(mod.assign)
# now adding in time varying variables. these always need to have time, either time+ or time:
# best from step 1.3 plus creche age
step6_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_crage_p.", mod.name), " <- run_penguin_models(phi.stru = \"", Phi, " + cr.age\")", sep = "")) %>% 
  select(mod.assign) 
# best from step 1.3 plus creching mass
step6_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_crmass_p.", mod.name), " <- run_penguin_models(phi.stru = \"", Phi, " + cr.mass\")", sep = "")) %>% 
  select(mod.assign)
# best from step 1.3 plus creching tib
step6_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_crtib_p.", mod.name), " <- run_penguin_models(phi.stru = \"", Phi, " + cr.tib\")", sep = "")) %>% 
  select(mod.assign)
# best from step 1.3 plus creching flip
step6_cand_set_base %>% 
  mutate(mod.assign = paste(gsub("_p.", "_crflip_p.", mod.name), " <- run_penguin_models(phi.stru = \"", Phi, " + cr.flip\")", sep = "")) %>% 
  select(mod.assign)

#
# make code for collect.models call ----
rbind(step6_cand_set_base %>% 
  select(mod.name),
# now adding in time varying variables. these always need to have time, either time+ or time:
# best from step 1.3 plus additive time and time varying creched.	#same intercept, diff slope for in.cr
data.frame(mod.name = "#crage"),
step6_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_crage_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus additive time and time varying age. same intercept, diff slope for dayold
data.frame(mod.name = "#crmass"),
step6_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_crmass_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus interaction of time and time varying creched.	#same slope, diff intercept for in.cr
data.frame(mod.name = "#crtib"),
step6_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_crtib_p.", mod.name)) %>% 
  select(mod.name),
# best from step 1.3 plus additive time and time varying age. same slope, diff intercept for dayold
data.frame(mod.name = "#crflip"),
step6_cand_set_base %>% 
  mutate(mod.name = gsub("_p.", "_crflip_p.", mod.name)) %>% 
  select(mod.name)) %>% 
  summarise(collect.call = paste(mod.name, collapse = "\", \"")) %>% 
  mutate(collect.call = paste("step6_cr_timing_surv = collect.models(lx=c(\"", collect.call, "\"))", sep = ""))

# now actually fitting models ----

step5.4mods <- readRDS(here("fitted_models/survival/step5.4_creched_base"))

# best step 5.3 
phi.year_T_hatch_p.sat <- step5.4mods$phi.year_T_hatch_p.sat
phi.year_lnT_hatch_p.sat <- step5.4mods$phi.year_lnT_hatch_p.sat
phi.year_T_p.sat <- step5.4mods$phi.year_T_p.sat
phi.year_lnT_p.sat <- step5.4mods$phi.year_lnT_p.sat
# best step 5.3 + creche age 
phi.year_T_hatch_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.age")
phi.year_lnT_hatch_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.age")
phi.year_T_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.age")
phi.year_lnT_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.age")
# best step 5.3 + creching mass 
phi.year_T_hatch_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.mass")
phi.year_lnT_hatch_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.mass")
phi.year_T_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.mass")
phi.year_lnT_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.mass")
# best step 5.3 + creching tibiotarsus 
phi.year_T_hatch_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.tib")
phi.year_lnT_hatch_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.tib")
phi.year_T_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.tib")
phi.year_lnT_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.tib")
# best step 5.3 + creching flipper 
phi.year_T_hatch_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.flip")
phi.year_lnT_hatch_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.flip")
phi.year_T_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.flip")
phi.year_lnT_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.flip")


step6_cr_timing_surv = collect.models(lx=c("phi.year_T_hatch_p.sat", "phi.year_lnT_hatch_p.sat", "phi.year_T_p.sat", "phi.year_lnT_p.sat", "#crage", "phi.year_T_hatch_crage_p.sat", "phi.year_lnT_hatch_crage_p.sat", "phi.year_T_crage_p.sat", "phi.year_lnT_crage_p.sat", "#crmass", "phi.year_T_hatch_crmass_p.sat", "phi.year_lnT_hatch_crmass_p.sat", "phi.year_T_crmass_p.sat", "phi.year_lnT_crmass_p.sat", "#crtib", "phi.year_T_hatch_crtib_p.sat", "phi.year_lnT_hatch_crtib_p.sat", "phi.year_T_crtib_p.sat", "phi.year_lnT_crtib_p.sat", "#crflip", "phi.year_T_hatch_crflip_p.sat", "phi.year_lnT_hatch_crflip_p.sat", "phi.year_T_crflip_p.sat", "phi.year_lnT_crflip_p.sat"))



step6_cr_timing_surv$model.table %>% 
  view()

saveRDS(step6_cr_timing_surv, here("fitted_models/survival/step6_cr_timing_surv"))

# check for uninformative parms in step 6 ----

step6 <- readRDS(here("fitted_models/survival/step6_cr_timing_surv"))

step6_mods <- step6[-length(step6)]

step6_informative <- map2_df(step6_mods, "Phi", cjs_parm_informative)

# we want to indicate uninformative parms in model selection tables. It ended up being simplest to just make the entire output model names here rather than when creating those tables, so we just run all models through ..._parm_informative (not just DAICc < 5)

# p_informative has a row for each parameter in each model, but note it isn't meaningful to check each time varying parm for un vs. informative so these are excluded
step6_informative <- step6_informative %>% 
  mutate(parm = gsub(":", ".", parm),
         parm = sub("resid.", "", parm),
         parm = sub("sexMale", "sex", parm),
         parm = sub("SEASON1314", "SEASON", parm),
         parm = sub("log\\(Time \\+ 1\\)", "lnT", parm),
         parm = sub("I\\(Time\\^2\\)", "TT", parm),
         parm = sub("Time", "T", parm)) %>%
  full_join(., step6$model.table %>% 
              data.frame() %>% 
              select(model.name = model, DeltaQAICc)) %>% 
  mutate(informative85 = ifelse(DeltaQAICc > 5, TRUE, informative85),
         parm = ifelse(grepl("~1", model.name), "p.intercept", parm),
         parm = ifelse(grepl("p\\(\\~time:in.cr\\)", model.name), "p.cr", parm))


         

step6_informative_wide <- step6_informative %>% 
           filter(!is.na(parm)) %>% 
  pivot_wider(id_cols = model.name, names_from = parm, values_from = informative85) %>% 
  mutate(across(c("Phi.T", "Phi.res.htch", "Phi.lnT", "Phi.cr.age", "Phi.cr.mass", "Phi.cr.tib", "Phi.cr.flip"), ~replace_na(., TRUE))) %>% 
  mutate(phi.stru = gsub("\\)p\\(\\~SEASON \\* sex \\+ time\\)", "", model.name),
         phi.stru = gsub(":", ".", phi.stru),
         phi.stru = sub("resid.", "", phi.stru),
         phi.stru = sub("sexMale", "sex", phi.stru),
         phi.stru = sub("SEASON1314", "SEASON", phi.stru),
         phi.stru = sub("log\\(Time \\+ 1\\)", "lnT", phi.stru),
         phi.stru = sub("Time", "T", phi.stru),
         phi.stru = sub("time", "t", phi.stru),
         phi.stru = sub("res.htch", "Hatch date", phi.stru),
         phi.stru = sub("cr.age", paste("Cr", "\u00E8", "ching age", sep = ""), phi.stru),
         phi.stru = sub("cr.mass", paste("Cr", "\u00E8", "ching mass", sep = ""), phi.stru),
         phi.stru = sub("cr.flip", paste("Cr", "\u00E8", "ching flipper", sep = ""), phi.stru),
         phi.stru = sub("cr.tib", paste("Cr", "\u00E8", "ching tibio", sep = ""), phi.stru),
         phi.stru = sub("Phi\\(\\~", "", phi.stru),
         phi.stru = ifelse(Phi.SEASON == FALSE, str_replace(phi.stru, "SEASON", paste("SEASON", "^\u2020^", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.res.htch == FALSE, str_replace(phi.stru, "Hatch date", paste("Hatch date", "^\u2020^", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.T == FALSE & !grepl("TT", phi.stru), str_replace(phi.stru, "T", paste("T", "^\u2020^", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.lnT == FALSE, str_replace(phi.stru, "lnT", paste("lnT", "^\u2020^", sep = "")), phi.stru),
         phi.stru = ifelse(Phi.cr.age == FALSE, str_replace(phi.stru, paste("Cr", "\u00E8", "ching age", sep = ""), sprintf(paste("Cr", "\u00E8", "ching age", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.cr.mass == FALSE, str_replace(phi.stru, paste("Cr", "\u00E8", "ching mass", sep = ""), sprintf(paste("Cr", "\u00E8", "ching mass", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.cr.flip == FALSE, str_replace(phi.stru, paste("Cr", "\u00E8", "ching flipper", sep = ""), sprintf(paste("Cr", "\u00E8", "ching flipper", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(Phi.cr.tib == FALSE, str_replace(phi.stru, paste("Cr", "\u00E8", "ching tibio", sep = ""), sprintf(paste("Cr", "\u00E8", "ching tibio", "^\u2020^", sep = ""))), phi.stru),
         phi.stru = ifelse(phi.stru == "1", "Intercept only", phi.stru),
         phi.stru = mod_call_to_structure(phi.stru)) %>% 
  select(model = model.name, phi.stru)

step6$model.table <- step6$model.table %>% 
  data.frame() %>% 
  full_join(step6_informative_wide, by = c("model"))

step6$model.table %>% view()


saveRDS(step6, here("fitted_models/survival/step6_cr_timing_surv"))


