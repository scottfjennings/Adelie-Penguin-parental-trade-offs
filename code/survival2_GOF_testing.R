

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
					 "in.cr97", "in.cr98", "cr.age", "fail.age", "cr.mass", "res.htch", "cr.flip", "cr.tib"), use.comments = TRUE) 


		 
##	setting NA values to mean of that field					
# These are only missing for chicks that didn't creche; pretty sure this not needed. commenting out, if no problems running analysis then can delete 
# penguins$cr.mass[is.na(penguins$cr.mass)]=1357.6				
# penguins$cr.tib[is.na(penguins$cr.tib)]=111.5				
# penguins$cr.flip[is.na(penguins$cr.flip)]=114.13				

##	actually need cr.age for chicks that didn't cr set to something other than 99
##	trying mean cr.age now (calculated after taking out 99's)
#penguins$cr.age <- replace(penguins$cr.age, penguins$cr.age==99, penguins$fail.age)
#penguins$cr.age <- as.numeric(penguins$cr.age)
#penguins$cr.age[is.na(penguins$cr.age)]=penguins$fail.age		#20

#write.csv(penguins, "resighting_survival/temp_penguins.csv")

#penguins$never.cr=1-penguins$did.cr
#penguins$did.cr=as.numeric(penguins$did.cr)
#penguins$never.cr=as.numeric(penguins$never.cr)


# Next create the processed dataframe and the design data. We’ll use a group 
# variable for colony so it can be used in the set of models for Phi. Factor 
# variables (covariates with a small finite set of values) are best handled by using 
# them to define groups in the data. 
# need to set begin.time to 51 because RMark considers the first dayold to be the first release dayold, but really I want to consider recaptures beginning on dayold 52.
penguins.process = process.data(penguins, model = "CJS", groups = c("SEASON", "sex"), begin.time = 50) 
penguins.ddl = make.design.data(penguins.process)

penguins.ddl$Phi %>% 
  data_frame() %>% 
  select(par.index, sex, SEASON) %>% 
  saveRDS(here("data/seas_sex_par_index"))

# RUN TO HERE TO PREP DATA FOR SURVIVAL ANALYSIS



sex_seas_phi_pims <- penguins.ddl$Phi
saveRDS(sex_seas_phi_pims, here("fitted_models/sex_seas_phi_pims"))

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





#c hat by hand
peng.phi.gen.p.gen$results$deviance/peng.phi.gen.p.gen$results$deviance.df







