




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


# filter out chicks that did not creche or did not get crech size estimated
penguins <- penguins %>% 
  filter(!is.na(cr.age), !is.na(cr.mass))
	

penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 50) 
penguins.ddl = make.design.data(penguins.process)

# best model for overall survival is Phi(~SEASON + res.htch + Time + I(Time^2))p(~SEASON + time)




# function to run models
run.phi.model<-function(phi.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru, p = list(formula = ~SEASON + time)), output = FALSE, chat = 1.25)
}

# how does the size and age at which chicks enter the creche stage relate to survival during creche stage?
# 7/9/21 there are a bunch of models here that I'm not sure should be included. I can't find any documentation why these models were added, and they don't appear in any draft. commenting them out for now, only keeping the ones that are in all the ms drafts

#Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~SEASON * sex + res.htch + Time)) #best from overall survival; not sure it make sense for this to be here; same data?
# Phi.dot	 =	run.phi.model(list(formula = ~1))
					
Phi.cr.mass	 =	run.phi.model(list(formula	 =	~ cr.mass)	) # cr mass
Phi.cr.flip	 =	run.phi.model(list(formula	 =	~ cr.flip)	) # cr flip
Phi.cr.tib	 =	run.phi.model(list(formula	 =	~ cr.tib)	) # cr tib
Phi.crage	 =	run.phi.model(list(formula	 =	~ cr.age)	) # cr age

Phi.SEASON_hatch_TT = run.phi.model(list(formula = ~ SEASON + res.htch + Time + I(Time^2))) # best overall survival
Phi.SEASON_hatch_TT_tibgr = run.phi.model(list(formula = ~ SEASON + res.htch + Time + I(Time^2) + tib.gr)) # best growth survival


Phi.SEASON_hatch_TT_crmass = run.phi.model(list(formula =	~ SEASON + res.htch + Time + I(Time^2) + cr.mass)) # cr mass + best overall
Phi.SEASON_hatch_TT_crflip =	run.phi.model(list(formula = ~ SEASON + res.htch + Time + I(Time^2) + cr.flip)) # cr flip + best overall
Phi.SEASON_hatch_TT_crtib = run.phi.model(list(formula = ~ SEASON + res.htch + Time + I(Time^2) + cr.tib)) # cr tib + best overall
Phi.SEASON_hatch_TT_crage = run.phi.model(list(formula =	~ SEASON + res.htch + Time + I(Time^2) + cr.age)) # cr age + best overall
					
Phi.dot = run.phi.model(list(formula = ~1))

cr_age_size_models = collect.models(lx=c("Phi.cr.mass", "Phi.cr.flip", "Phi.cr.tib", "Phi.crage", "Phi.SEASON_hatch_TT", "Phi.SEASON_hatch_TT_tibgr", "Phi.SEASON_hatch_TT_crmass", "Phi.SEASON_hatch_TT_crtib", "Phi.SEASON_hatch_TT_crflip", "Phi.SEASON_hatch_TT_crage", "Phi.dot"))

model.table(cr_age_size_models)


saveRDS(cr_age_size_models, here("fitted_models/survival/cr_age_size_models"))

