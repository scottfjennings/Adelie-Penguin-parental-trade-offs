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




# goes to root folder for project -- not sure why need to output these so commented out for now
# export.MARK(penguins.process, "processed_penguins_ch")
# export.chdata(penguins.process, filename="penguin", replace=TRUE)

## GOF testing ----
# GOF testing with RELEASE requires no individual covariates and 
# requires dots in the capture history to be replaced with 0. dot -> 0 is supposed to happen automatically, but for some reason isn't, so doing it manually
peng.for.gof=subset(penguins, select=c("ch", "freq", "sex", "SEASON"))

peng.for.gof <- peng.for.gof %>% 
  mutate(ch = gsub("\\.", "0", ch))

peng.for.gof.proc=process.data(peng.for.gof,model="CJS",groups=c("SEASON", "sex"))

release.gof(peng.for.gof.proc)
# note, release.gof creates mxxx.tmp file
# divide total Chi.square by total df to get estimate of c-hat
# > 179.5428/144
# [1] 1.246825
c_hat = round(179.5428/144, 2)

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

# peng.p <- readRDS(here("fitted_models/p_mod_list"))


	
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


# phi_aictable1 <- model.table(phi_step1, use.lnl=TRUE) 
# saveRDS(phi_aictable1, here("fitted_models/phi_aictable1"))
# 7/8/21 these models run and saved

## phi	Step 2: take best structures from above and investigate hatch date and time effects ----
##	including time constraints based on time-varying covariates
#inc r =time-varying for in cr or not on specific day
#dayold = time-varying for age on specific day;

# considering 4 structures from step 2, those with DAIC <= 2, (would be 3 structures using models with model weight >= 0.1):

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

# ~ SEASON * sex + res.htch- step 1 DAIC = 0.458
Phi.SEASON.sex_hatch = run.phi.model(list(formula = ~ SEASON * sex + res.htch))
Phi.SEASON.sex_hatch_t = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time))	
Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time))
Phi.SEASON.sex_hatch_TT = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + I(Time^2)))
Phi.SEASON.sex_hatch_lnT = run.phi.model(list(formula = ~ SEASON * sex + res.htch + log(Time+1)))
Phi.SEASON.sex_hatch_t_incr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.SEASON.sex_hatch_t_age = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.SEASON.sex_hatch_t.incr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.SEASON.sex_hatch_t.age = run.phi.model(list(formula = ~ SEASON * sex + res.htch + time:dayold))		#8	same slope, diff intercept for dayold
	
# ~res.htch - step 1 DAIC = 1.384
Phi.hatch = run.phi.model(list(formula = ~res.htch))
Phi.hatch_t = run.phi.model(list(formula = ~ res.htch + time))	
Phi.hatch_T = run.phi.model(list(formula = ~ res.htch + Time))
Phi.hatch_TT = run.phi.model(list(formula = ~res.htch + Time + I(Time^2)))
Phi.hatch_lnT = run.phi.model(list(formula = ~ res.htch + log(Time+1)))
Phi.hatch_t_incr = run.phi.model(list(formula = ~ res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.hatch_t_age = run.phi.model(list(formula = ~ res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.hatch_t.incr = run.phi.model(list(formula = ~ res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.hatch_t.age = run.phi.model(list(formula = ~ res.htch + time:dayold))

# SEASON + sex + res.htch - step 1 DAIC = 1.987
Phi.SEASON_sex_hatch = run.phi.model(list(formula = ~ SEASON + sex + res.htch))
Phi.SEASON_sex_hatch_t = run.phi.model(list(formula = ~ SEASON + sex + res.htch + time))	
Phi.SEASON_sex_hatch_T = run.phi.model(list(formula = ~ SEASON + sex + res.htch + Time))
Phi.SEASON_sex_hatch_TT = run.phi.model(list(formula = ~ SEASON + sex + res.htch + Time + I(Time^2)))
Phi.SEASON_sex_hatch_lnT = run.phi.model(list(formula = ~ SEASON + sex + res.htch + log(Time+1)))
Phi.SEASON_sex_hatch_t_incr = run.phi.model(list(formula = ~ SEASON + sex + res.htch + time + in.cr))		#5	#same intercept, diff slope for in.cr
Phi.SEASON_sex_hatch_t_age = run.phi.model(list(formula = ~ SEASON + sex + res.htch + time + dayold))		#6	same intercept, diff slope for dayold
Phi.SEASON_sex_hatch_t.incr = run.phi.model(list(formula = ~ SEASON + sex + res.htch + time:in.cr))		#same slope, diff intercept for in.cr
Phi.SEASON_sex_hatch_t.age = run.phi.model(list(formula = ~ SEASON + sex + res.htch + time:dayold))

phi_step2=collect.models(lx=c("Phi.SEASON_hatch", "Phi.SEASON_hatch_t", "Phi.SEASON_hatch_T", "Phi.SEASON_hatch_TT", "Phi.SEASON_hatch_lnT", "Phi.SEASON_hatch_t_incr", "Phi.SEASON_hatch_t_age", "Phi.SEASON_hatch_t.incr", "Phi.SEASON_hatch_t.age",
                              #
                              "Phi.SEASON.sex_hatch", "Phi.SEASON.sex_hatch_t", "Phi.SEASON.sex_hatch_T", "Phi.SEASON.sex_hatch_TT", "Phi.SEASON.sex_hatch_lnT", "Phi.SEASON.sex_hatch_t_incr", "Phi.SEASON.sex_hatch_t_age", "Phi.SEASON.sex_hatch_t.incr", "Phi.SEASON.sex_hatch_t.age",
                              #
                              "Phi.hatch", "Phi.hatch_t", "Phi.hatch_T", "Phi.hatch_TT", "Phi.hatch_lnT", "Phi.hatch_t_incr", "Phi.hatch_t_age", "Phi.hatch_t.incr", "Phi.hatch_t.age",
                              #
                              "Phi.SEASON_sex_hatch", "Phi.SEASON_sex_hatch_t", "Phi.SEASON_sex_hatch_T", "Phi.SEASON_sex_hatch_TT", "Phi.SEASON_sex_hatch_lnT", "Phi.SEASON_sex_hatch_t_incr", "Phi.SEASON_sex_hatch_t_age", "Phi.SEASON_sex_hatch_t.incr", "Phi.SEASON_sex_hatch_t.age"))




saveRDS(phi_step2, here("fitted_models/phi_step2"))

phi_step2 <- readRDS(here("fitted_models/phi_step2"))


#----

step2_mod_av = model.average(phi_step2, "Phi", vcv=TRUE, drop = FALSE)

saveRDS(step2_mod_av, here("fitted_models/step2_mod_av"))

