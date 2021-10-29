




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

# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/

## !!!NOTE 8/26/21- inp file now almost entirely created by code, including in.cr fields. See CH4_make_inp.R for code and notes 


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
zmod <- mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru.list, p = p.stru.list), output = FALSE, chat = 1.25)
return(zmod)
}

cjs_parm_informative <- function(zmodel, phi.p) {
  zmod.name = deparse(substitute(zmodel))
  
  coef(zmodel) %>% 
  data.frame() %>% 
  rownames_to_column("parm") %>% 
  filter(grepl(paste(phi.p, ":", sep = ""), parm), !grepl("ntercept", parm), !grepl("time", parm)) %>% 
    mutate(lci85 = estimate - (1.44 * se),
           uci85 = estimate + (1.44 * se),
           ci95check = lcl == estimate - (1.96 * se) & ucl == estimate + (1.96 * se),
           informative95 = ifelse((lcl < 0 & ucl < 0) | (lcl > 0 & ucl > 0), TRUE, FALSE),
           informative85 = ifelse((lci85 < 0 & uci85 < 0) | (lci85 > 0 & uci85 > 0), TRUE, FALSE),
           mod.name = zmod.name)
}
#


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
###	now moving on to modeling Phi

# for my thesis and first round of publication draft, Phi structure began with 2 steps: first just SEASON and Sex structure (table 4.4 in thesis; table 3b in pub ms), then relative hatch date (tables 4.5 and 3c).
# in editing ms version we decided to combine these 2 steps for a cleaner and more readable analysis; this required fitting more models than the sum of the 2 original candidate sets.

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

model.table(step5.2_phi)


saveRDS(step5.2_phi, here("fitted_models/survival/step5.2_creched_phi"))

# phi_step5_2 <- readRDS(here("fitted_models/survival/phi_step5_models"))


#
# Step 5.3: create new candidate set from the best structures from above ----
# step 5.3.1: first check for uninformative parms where DQAICC < 5 in step 5.1  ----

step5.1_p <- readRDS(here("fitted_models/survival/step5.1_creched_p"))


step5.1_p %>% 
  model.table(model.name = FALSE) %>% 
  filter(DeltaQAICc <= 5) %>% 
  mutate(uninform.assign = paste("cjs_parm_informative(step5.1_p$", model, ", \"p\")", sep = "")) %>% 
  summarise(uninform.assign = paste(uninform.assign, collapse = ", ")) %>% 
  mutate(uninform.assign = paste("p_informative <- rbind(", uninform.assign, ")", sep = ""))

p_informative <- rbind(cjs_parm_informative(step5.1_p$phi.sat_p.year_t, "p"), 
                       cjs_parm_informative(step5.1_p$phi.sat_p.year_sex_t, "p"), 
                       cjs_parm_informative(step5.1_p$phi.sat_p.year.sex_t, "p"))

p_informative <- p_informative %>% 
  mutate(mod.name = gsub("step5.1_p\\$", "", mod.name))

p_informative_wide <- p_informative %>% 
  pivot_wider(id_cols = mod.name, names_from = parm, values_from = informative85) %>% 
  mutate(uninformative = "")

p_informative_wide <- edit(p_informative_wide)

p_informative_wide %>%
  filter(uninformative == TRUE) %>%  
  select(mod.name, uninformative) %>% 
  saveRDS(here("fitted_models/survival/step5.1_uninformative"))


# year and year:sex parms uninformative. only moving phi.sat_p.year_t forward

# step 5.3.2: then check for uninformative parms where DQAICC < 5 in step 5.2  ----
step5.2_phi <- readRDS(here("fitted_models/survival/step5.2_creched_phi"))


step5.2_phi %>% 
  model.table(model.name = FALSE) %>% 
  filter(DeltaQAICc <= 5) %>% 
  mutate(uninform.assign = paste("cjs_parm_informative(step5.2_phi$", model, ", \"Phi\")", sep = "")) %>% 
  summarise(uninform.assign = paste(uninform.assign, collapse = ", ")) %>% 
  mutate(uninform.assign = paste("phi_informative <- rbind(", uninform.assign, ")", sep = ""))



# check for uninformative parms in phi structure models
phi_informative <- rbind(
  cjs_parm_informative(step5.2_phi$phi.year_T_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_lnT_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_T_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_lnT_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_T_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_TT_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_lnT_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_T_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_lnT_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year.sex_T_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_TT_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year.sex_lnT_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year.sex_T_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year.sex_lnT_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_TT_hatch_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year_sex_TT_p.sat, "Phi"),
  cjs_parm_informative(step5.2_phi$phi.year.sex_TT_hatch_p.sat, "Phi"))

phi_informative <- phi_informative %>% 
  mutate(mod.name = gsub("step5.2_phi\\$", "", mod.name))


phi_informative_wide <- phi_informative %>% 
  pivot_wider(id_cols = mod.name, names_from = parm, values_from = informative85) %>% 
  mutate(uninformative = "")
# can use edit() to quickly, manually select which models are uninformative
phi_informative_wide <- edit(phi_informative_wide)
# sex and year:sex always uninformative
# T2 always uninformative, T1 uninformative when included with T2

# big difference here from step 1.3 is that lnT is informative


phi_informative_wide %>%
  filter(uninformative == TRUE) %>%  
  select(mod.name, uninformative) %>% 
  saveRDS(here("fitted_models/survival/step5.2_uninformative"))

# 5.3.3 then collecting the models into step 5.3 ----
# don't need to refit, just get the appropriate models from step5.1_creched_p and step5.2_creched_phi

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

step5.3 = collect.models(lx=c("phi.sat_p.year_t", "phi.year_T_hatch_p.sat", "phi.year_lnT_hatch_p.sat", "phi.year_T_p.sat", "phi.year_lnT_p.sat"))


model.table(step5.3)
saveRDS(step5.3, here("fitted_models/survival/step5.3_creched_base"))

# step 6 add creching age and size ----
# helper code to build up candidata set code ----

step5.3mods <- readRDS(here("fitted_models/survival/step5.3_creched_base"))

step6_cand_set_base <- step5.3mods %>% 
  model.table() %>% 
  data.frame() %>%
  rownames_to_column("mod.num") %>% 
  filter(DeltaQAICc <= 5) %>% 
  left_join(., names(step5.3mods) %>%
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
# best step 5.3 ----
phi.year_T_hatch_p.sat <- step5.3mods$phi.year_T_hatch_p.sat
phi.year_lnT_hatch_p.sat <- step5.3mods$phi.year_lnT_hatch_p.sat
phi.year_T_p.sat <- step5.3mods$phi.year_T_p.sat
phi.year_lnT_p.sat <- step5.3mods$phi.year_lnT_p.sat
# best step 5.3 + creche age ----
phi.year_T_hatch_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.age")
phi.year_lnT_hatch_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.age")
phi.year_T_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.age")
phi.year_lnT_crage_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.age")
# best step 5.3 + creching mass ----
phi.year_T_hatch_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.mass")
phi.year_lnT_hatch_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.mass")
phi.year_T_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.mass")
phi.year_lnT_crmass_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.mass")
# best step 5.3 + creching tibiotarsus ----
phi.year_T_hatch_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.tib")
phi.year_lnT_hatch_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.tib")
phi.year_T_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.tib")
phi.year_lnT_crtib_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.tib")
# best step 5.3 + creching flipper ----
phi.year_T_hatch_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + res.htch + cr.flip")
phi.year_lnT_hatch_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + res.htch + cr.flip")
phi.year_T_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + Time + cr.flip")
phi.year_lnT_crflip_p.sat <- run_penguin_models(phi.stru = "SEASON + log(Time + 1) + cr.flip")

step6_cr_timing_surv = collect.models(lx=c("phi.year_T_hatch_p.sat", "phi.year_lnT_hatch_p.sat", "phi.year_T_p.sat", "phi.year_lnT_p.sat", "#crage", "phi.year_T_hatch_crage_p.sat", "phi.year_lnT_hatch_crage_p.sat", "phi.year_T_crage_p.sat", "phi.year_lnT_crage_p.sat", "#crmass", "phi.year_T_hatch_crmass_p.sat", "phi.year_lnT_hatch_crmass_p.sat", "phi.year_T_crmass_p.sat", "phi.year_lnT_crmass_p.sat", "#crtib", "phi.year_T_hatch_crtib_p.sat", "phi.year_lnT_hatch_crtib_p.sat", "phi.year_T_crtib_p.sat", "phi.year_lnT_crtib_p.sat", "#crflip", "phi.year_T_hatch_crflip_p.sat", "phi.year_lnT_hatch_crflip_p.sat", "phi.year_T_crflip_p.sat", "phi.year_lnT_crflip_p.sat"))



step6_cr_timing_surv$model.table %>% 
  select(-model)

saveRDS(step6_cr_timing_surv, here("fitted_models/survival/step6_cr_timing_surv"))

