

# run the same candidate set as in ch4_creche_age_size_survival.R, but only for chicks that actually reached the creche stage. 
# so can't compare these models to the best models from previous steps


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
# .inp file is created by hand in a spreadsheet editor, then saved at a .txt, then saved as a .inp

# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/
penguins=convert.inp(here("data/all_covs_sex_groups.inp"), 
					group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), 
					covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56", 
					 "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64", 
					 "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72",
					 "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80",
					 "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", 
					 "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96",
					 "dayold97", "dayold98", "av.food", "av.trip.length", "cr.age", "fail.age", "weight", 
					"flipper", "tibiotar", "CH_B", "CH_S", "mean.N", "cr.mass", 
					"res.htch", "cr.flip", "cr.tib", "did.cr", "fate.age" ),use.comments=FALSE) 

		 
penguins = filter(penguins, did.cr == 1)


# Next create the processed dataframe and the design data. We’ll use a group 
# variable for colony so it can be used in the set of models for Phi. Factor 
# variables (covariates with a small finite set of values) are best handled by using 
# them to define groups in the data. 
# need to set begin.time to 51 because RMark considers the first dayold to be the first release dayold, but really I want to consider recaptures beginning on dayold 52.
penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 50) 
penguins.ddl = make.design.data(penguins.process)

# function to run models
run.phi.model<-function(phi.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru, p = list(formula = ~SEASON + time)), output = FALSE, chat = 1.25)
}

# how does the size and age at which chicks enter the creche stage relate to survival during creche stage?
# 7/9/21 there are a bunch of models here that I'm not sure should be included. I can't find any documentation why these models were added, and they don't appear in any draft. commenting them out for now, only keeping the ones that are in all the ms drafts

#Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~SEASON * sex + res.htch + Time)) #best from overall survival; not sure it make sense for this to be here; same data?
# Phi.dot	 =	run.phi.model(list(formula = ~1))
					
Phi.cr.mass	 =	run.phi.model(list(formula	 =	~cr.mass)	) # cr mass
Phi.cr.flip	 =	run.phi.model(list(formula	 =	~cr.flip)	) # cr flip
Phi.cr.tib	 =	run.phi.model(list(formula	 =	~cr.tib)	) # cr tib
Phi.fate.age	 =	run.phi.model(list(formula	 =	~fate.age)	) # cr age

Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time)) # best overall survival
Phi.SEASON.sex_hatch_T_tibgr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + tibiotar)) # best growth survival


Phi.SEASON.sex_hatch_T_crmass = run.phi.model(list(formula =	~ SEASON * sex + res.htch + Time + cr.mass)) # cr mass + best overall
Phi.SEASON.sex_hatch_T_crflip =	run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + cr.flip)) # cr flip + best overall
Phi.SEASON.sex_hatch_T_crtib = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + cr.tib)) # cr tib + best overall
Phi.SEASON.sex_hatch_T_fateage = run.phi.model(list(formula =	~ SEASON * sex + res.htch + Time + fate.age)) # cr age + best overall
					

cr_age_size_models = collect.models(lx=c("Phi.cr.mass", "Phi.cr.flip", "Phi.cr.tib", "Phi.fate.age", "Phi.SEASON.sex_hatch_T", "Phi.SEASON.sex_hatch_T_tibgr", "Phi.SEASON.sex_hatch_T_crmass", "Phi.SEASON.sex_hatch_T_crtib", "Phi.SEASON.sex_hatch_T_crflip"))

saveRDS(cr_age_size_models, here("fitted_models/cr_age_size_models_cronly"))

