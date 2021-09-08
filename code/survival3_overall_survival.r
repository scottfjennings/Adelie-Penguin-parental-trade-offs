


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
##	making models without mark.wrapper
##  I think this method of building and running models will be easier because of the npar adjustments I have to make
##  using mark.wrapper it is hard to keep track of which model is which, 
##  and I don't really have that many models to build.


# Determine p structure ----
run.p.models<-function(p.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = list(formula = ~SEASON * sex + time), p = p.stru), output = FALSE, chat = 1.25)
}
# resighting varies by interaction between season, sex and time -- most general model
peng.phi.gen.p.gen = run.p.models(list(formula = ~SEASON * sex + time))
# resighting varies by season
peng.phi.gen.p.SEASON = run.p.models(list(formula = ~SEASON))
# resighting varies by sex
peng.phi.gen.p.sex = run.p.models(list(formula = ~sex))
# resighting varies by time
peng.phi.gen.p.time = run.p.models(list(formula = ~time))
# resighting varies by season + sex
peng.phi.gen.p.SEASON.pl.sex = run.p.models(list(formula = ~SEASON + sex))
# resighting varies by season + time
peng.phi.gen.p.SEASON.pl.time = run.p.models(list(formula = ~SEASON + time))
# resighting varies by sex + time
peng.phi.gen.p.sex.pl.time = run.p.models(list(formula = ~sex + time))
# resighting does not vary -- intercept only
peng.phi.gen.p.dot = run.p.models(list(formula = ~1))
# resighting varies by linear time
peng.phi.gen.p.T = run.p.models(list(formula = ~Time))
# resighting varies by quadratic time
peng.phi.gen.p.TT = run.p.models(list(formula = ~Time + I(Time^2)))
# resighting varies by natural log time
peng.phi.gen.p.lnT = run.p.models(list(formula = ~log(Time+1)))
# resighting varies by whether or not chick was in creche stage on that day
peng.phi.gen.p.cr = run.p.models(list(formula = ~in.cr))


p_mod_names <- c("peng.phi.gen.p.gen", "peng.phi.gen.p.SEASON", "peng.phi.gen.p.sex", "peng.phi.gen.p.time", "peng.phi.gen.p.SEASON.pl.sex", "peng.phi.gen.p.SEASON.pl.time", "peng.phi.gen.p.sex.pl.time", "peng.phi.gen.p.dot", "peng.phi.gen.p.T", "peng.phi.gen.p.TT", "peng.phi.gen.p.lnT", "peng.phi.gen.p.cr")

##  now collect models into one object to make a model comparison table
peng.p=collect.models(p_mod_names)

saveRDS(peng.p, here("fitted_models/p_mod_list"))

# peng.p2 <- readRDS(here("fitted_models/p_mod_list"))


	
# determine phi structure ----
###	now moving on to modeling Phi, with the best p structure from above, which is SEASON+time

# for my thesis and first round of publication draft, Phi strucutre began with 2 steps: first just SEASON and Sex structure (table 4.4 in thesis; table 3b in pub ms), then relative hatch date (tables 4.5 and 3c).
# in editing ms version we decided to combine these 2 steps for a cleaner and more readable analysis; this required fitting more models than the sum of the 2 original candidate sets.


run.phi.model<-function(phi.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru, p = list(formula = ~SEASON + time)), output = FALSE, chat = 1.25)
}

## NO RUN	Phi OLD step 1: determine effect of SEASON and sex ----
# Phi.SEASON.sex = run.phi.model(list(formula = ~SEASON * sex))
# Phi.SEASON_sex = run.phi.model(list(formula = ~SEASON + sex))
# Phi.SEASON = run.phi.model(list(formula = ~SEASON))
# Phi.sex = run.phi.model(list(formula = ~sex))
# Phi.dot = run.phi.model(list(formula = ~1))


# phi_step1 = collect.models(lx=c("Phi.SEASON.sex", "Phi.SEASON_sex", "Phi.SEASON", "Phi.sex", "Phi.dot"))
# saveRDS(phi_step1, here("fitted_models/phi_step1"))
# phi_aictable1 <- model.table(phi_step1, use.lnl=TRUE) 
# saveRDS(phi_aictable1, here("fitted_models/phi_aictable1"))
# phi_aictable1 <- readRDS(here("fitted_models/phi_aictable1"))
# these models fitted and saved 7/8/21


## NO RUN phi OLD Step 2:	is relative hatch date important? ----
# Phi.y.s=run.phi.model(list(formula=~SEASON*sex))
# Phi.res.htch_y.s=run.phi.model(list(formula=~res.htch+SEASON*sex))
# Phi.res.htch.y.s=run.phi.model(list(formula=~res.htch*SEASON*sex))
# Phi.res.htch2_y.s=run.phi.model(list(formula=~res.htch+I(res.htch^2)+SEASON*sex))

# phi_step2=collect.models(lx=c("Phi.y.s", "Phi.res.htch_y.s", "Phi.res.htch.y.s", "Phi.res.htch2_y.s"))
# saveRDS(phi_step2, here("fitted_models/phi_step2"))
# phi_aictable2 <- model.table(phi_step2, use.lnl=TRUE) 
# saveRDS(phi_aictable2, here("fitted_models/phi_aictable2"))
# these models fitted and saved 7/8/21

# Phi NEW step 1; SEASON, Sex and relative hatch date ----
# just SEASON and Sex (this was the original step 1 candidate set)
Phi.SEASON.sex = run.phi.model(list(formula = ~SEASON * sex))
Phi.SEASON_sex = run.phi.model(list(formula = ~SEASON + sex))
Phi.SEASON = run.phi.model(list(formula = ~SEASON))
Phi.sex = run.phi.model(list(formula = ~sex))

# SEASON * Sex with 3 residual hatch date options
Phi.SEASON.sex_hatch = run.phi.model(list(formula = ~ SEASON * sex + res.htch))
Phi.SEASON.sex.hatch = run.phi.model(list(formula = ~ SEASON * sex * res.htch))
Phi.SEASON.sex_hatch2 = run.phi.model(list(formula = ~ SEASON * sex + res.htch + I(res.htch^2)))

# SEASON + Sex with 3 residual hatch date options
Phi.SEASON_sex_hatch = run.phi.model(list(formula = ~ SEASON + sex + res.htch))
Phi.SEASON_sex.hatch = run.phi.model(list(formula = ~ SEASON + sex * res.htch))
Phi.SEASON_sex_hatch2 = run.phi.model(list(formula = ~ SEASON + sex + res.htch + I(res.htch^2)))

# SEASON with 3 residual hatch date options
Phi.SEASON_hatch = run.phi.model(list(formula = ~ SEASON + res.htch))
Phi.SEASON.hatch = run.phi.model(list(formula = ~ SEASON * res.htch))
Phi.SEASON_hatch2 = run.phi.model(list(formula = ~ SEASON + res.htch + I(res.htch^2)))

# Sex with 3 residual hatch date options
Phi.Sex_hatch = run.phi.model(list(formula = ~ sex + res.htch))
Phi.Sex.hatch = run.phi.model(list(formula = ~ sex * res.htch))
Phi.Sex_hatch2 = run.phi.model(list(formula = ~ sex + res.htch + I(res.htch^2)))

# residual hatch date only
Phi.hatch = run.phi.model(list(formula = ~res.htch))
Phi.hatch2 = run.phi.model(list(formula = ~res.htch + I(res.htch^2)))

# intercept only
Phi.dot = run.phi.model(list(formula = ~1))


phi_step1=collect.models(lx=c("Phi.SEASON.sex", "Phi.SEASON_sex", "Phi.SEASON", "Phi.sex", 
                          "Phi.SEASON.sex_hatch", "Phi.SEASON.sex.hatch", "Phi.SEASON.sex_hatch2",
                          "Phi.SEASON_sex_hatch", "Phi.SEASON.sex.hatch", "Phi.SEASON_sex_hatch2",
                          "Phi.SEASON_hatch", "Phi.SEASON.hatch", "Phi.SEASON_hatch2",
                          "Phi.Sex_hatch", "Phi.Sex.hatch", "Phi.Sex_hatch2",
                          "Phi.hatch", "Phi.hatch2",
                          "Phi.dot"))

saveRDS(phi_step1, here("fitted_models/phi_step1"))
phi_step1_2 <- readRDS(here("fitted_models/phi_step1"))


# phi_aictable1 <- model.table(phi_step1, use.lnl=TRUE) 
# saveRDS(phi_aictable1, here("fitted_models/phi_aictable1"))
# 7/8/21 these models run and saved

## phi	Step 2: take best structures from above and investigate hatch date and time effects ----
##	including time constraints based on time-varying covariates
#inc r =time-varying for in cr or not on specific day
#dayold = time-varying for age on specific day;

# considering 3 structures from step 2, those with DAIC <= 2, (would be 3 structures using models with model weight >= 0.1):

# ~ SEASON + res.htch - step 1 DAICc = 0
Phi.SEASON_hatch = run.phi.model(list(formula = ~ SEASON + res.htch))
Phi.SEASON_hatch_t = run.phi.model(list(formula = ~ SEASON + res.htch + time))	
Phi.SEASON_hatch_T = run.phi.model(list(formula = ~ SEASON + res.htch + Time))
Phi.SEASON_hatch_TT = run.phi.model(list(formula = ~ SEASON + res.htch + Time + I(Time^2)))
Phi.SEASON_hatch_lnT = run.phi.model(list(formula = ~ SEASON + res.htch + log(Time+1)))
Phi.SEASON_hatch_t_incr = run.phi.model(list(formula = ~ SEASON + res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.SEASON_hatch_t_age = run.phi.model(list(formula = ~ SEASON + res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.SEASON_hatch_t.incr = run.phi.model(list(formula = ~ SEASON + res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.SEASON_hatch_t.age = run.phi.model(list(formula = ~ SEASON + res.htch + time:dayold))		#8	same slope, diff intercept for dayold

# ~ SEASON * sex + res.htch- step 1 DAIC = 1.126095
Phi.SEASON.sex_hatch = run.phi.model(list(formula = ~ SEASON * sex + res.htch))
Phi.SEASON.sex_hatch_t = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time))	
Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time))
Phi.SEASON.sex_hatch_TT = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + I(Time^2)))
Phi.SEASON.sex_hatch_lnT = run.phi.model(list(formula = ~ SEASON * sex + res.htch + log(Time+1)))
Phi.SEASON.sex_hatch_t_incr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.SEASON.sex_hatch_t_age = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.SEASON.sex_hatch_t.incr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.SEASON.sex_hatch_t.age = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time:dayold))		#8	same slope, diff intercept for dayold
	
# ~SEASON * res.htch - step 1 DAIC = 1.984
Phi.SEASON.hatch = run.phi.model(list(formula = ~SEASON * res.htch))
Phi.SEASON.hatch_t = run.phi.model(list(formula = ~ SEASON * res.htch + time))	
Phi.SEASON.hatch_T = run.phi.model(list(formula = ~ SEASON * res.htch + Time))
Phi.SEASON.hatch_TT = run.phi.model(list(formula = ~ SEASON * res.htch + Time + I(Time^2)))
Phi.SEASON.hatch_lnT = run.phi.model(list(formula = ~ SEASON * res.htch + log(Time+1)))
Phi.SEASON.hatch_t_incr = run.phi.model(list(formula = ~ SEASON * res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.SEASON.hatch_t_age = run.phi.model(list(formula = ~ SEASON * res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.SEASON.hatch_t.incr = run.phi.model(list(formula = ~ SEASON * res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.SEASON.hatch_t.age = run.phi.model(list(formula = ~ SEASON * res.htch + time:dayold))



phi_step2=collect.models(lx=c("Phi.SEASON_hatch", "Phi.SEASON_hatch_t", "Phi.SEASON_hatch_T", "Phi.SEASON_hatch_TT", "Phi.SEASON_hatch_lnT", "Phi.SEASON_hatch_t_incr", "Phi.SEASON_hatch_t_age", "Phi.SEASON_hatch_t.incr", "Phi.SEASON_hatch_t.age",
                              #
                              "Phi.SEASON.sex_hatch", "Phi.SEASON.sex_hatch_t", "Phi.SEASON.sex_hatch_T", "Phi.SEASON.sex_hatch_TT", "Phi.SEASON.sex_hatch_lnT", "Phi.SEASON.sex_hatch_t_incr", "Phi.SEASON.sex_hatch_t_age", "Phi.SEASON.sex_hatch_t.incr", "Phi.SEASON.sex_hatch_t.age",
                              #
                              "Phi.SEASON.hatch", "Phi.SEASON.hatch_t", "Phi.SEASON.hatch_T", "Phi.SEASON.hatch_TT", "Phi.SEASON.hatch_lnT", "Phi.SEASON.hatch_t_incr", "Phi.SEASON.hatch_t_age", "Phi.SEASON.hatch_t.incr", "Phi.SEASON.hatch_t.age"))

saveRDS(phi_step2, here("fitted_models/phi_step2"))

# model averaging ----

step2_mod_av = model.average(phi_step2, "Phi", vcv=TRUE, drop = FALSE)

phi_step2 <- readRDS(here("fitted_models/phi_step2"))

saveRDS(step2_mod_av, here("fitted_models/step2_mod_av"))

