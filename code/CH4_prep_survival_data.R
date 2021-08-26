

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
# .inp file is created by hand in a spreadsheet editor, then saved at a .txt, then saved as a .inp

# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/

## !!!NOTE 8/25/21-  did.cr in all_covs_sex_groups.inp is NOT correct. I am currently dealing with this downstream of reading in all_covs_sex_groups.inp, rather than remaking it.

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

		 
##	setting NA values to mean of that field					

penguins$cr.mass[is.na(penguins$cr.mass)]=1357.6				
penguins$cr.tib[is.na(penguins$cr.tib)]=111.5				
penguins$cr.flip[is.na(penguins$cr.flip)]=114.13				
penguins$mean.N[is.na(penguins$mean.N)]=10.76

##	actually need cr.age for chicks that didn't cr set to something other than 99
##	trying mean cr.age now (calculated after taking out 99's)
#penguins$cr.age <- replace(penguins$cr.age, penguins$cr.age==99, penguins$fail.age)
#penguins$cr.age <- as.numeric(penguins$cr.age)
#penguins$cr.age[is.na(penguins$cr.age)]=penguins$fail.age		#20

#write.csv(penguins, "resighting_survival/temp_penguins.csv")

## add time-varying covariates for creched or not
penguins$in.cr50=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold50<penguins$cr.age,0,1))
penguins$in.cr51=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold51<penguins$cr.age,0,1))
penguins$in.cr52=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold52<penguins$cr.age,0,1))
penguins$in.cr53=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold53<penguins$cr.age,0,1))
penguins$in.cr54=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold54<penguins$cr.age,0,1))
penguins$in.cr55=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold55<penguins$cr.age,0,1))
penguins$in.cr56=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold56<penguins$cr.age,0,1))
penguins$in.cr57=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold57<penguins$cr.age,0,1))
penguins$in.cr58=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold58<penguins$cr.age,0,1))
penguins$in.cr59=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold59<penguins$cr.age,0,1))
penguins$in.cr60=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold60<penguins$cr.age,0,1))
penguins$in.cr61=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold61<penguins$cr.age,0,1))
penguins$in.cr62=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold62<penguins$cr.age,0,1))
penguins$in.cr63=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold63<penguins$cr.age,0,1))
penguins$in.cr64=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold64<penguins$cr.age,0,1))
penguins$in.cr65=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold65<penguins$cr.age,0,1))
penguins$in.cr66=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold66<penguins$cr.age,0,1))
penguins$in.cr67=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold67<penguins$cr.age,0,1))
penguins$in.cr68=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold68<penguins$cr.age,0,1))
penguins$in.cr69=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold69<penguins$cr.age,0,1))
penguins$in.cr70=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold70<penguins$cr.age,0,1))
penguins$in.cr71=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold71<penguins$cr.age,0,1))
penguins$in.cr72=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold72<penguins$cr.age,0,1))
penguins$in.cr73=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold73<penguins$cr.age,0,1))
penguins$in.cr74=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold74<penguins$cr.age,0,1))
penguins$in.cr75=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold75<penguins$cr.age,0,1))
penguins$in.cr76=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold76<penguins$cr.age,0,1))
penguins$in.cr77=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold77<penguins$cr.age,0,1))
penguins$in.cr78=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold78<penguins$cr.age,0,1))
penguins$in.cr79=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold79<penguins$cr.age,0,1))
penguins$in.cr80=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold80<penguins$cr.age,0,1))
penguins$in.cr81=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold81<penguins$cr.age,0,1))
penguins$in.cr82=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold82<penguins$cr.age,0,1))
penguins$in.cr83=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold83<penguins$cr.age,0,1))
penguins$in.cr84=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold84<penguins$cr.age,0,1))
penguins$in.cr85=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold85<penguins$cr.age,0,1))
penguins$in.cr86=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold86<penguins$cr.age,0,1))
penguins$in.cr87=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold87<penguins$cr.age,0,1))
penguins$in.cr88=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold88<penguins$cr.age,0,1))
penguins$in.cr89=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold89<penguins$cr.age,0,1))
penguins$in.cr90=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold90<penguins$cr.age,0,1))
penguins$in.cr91=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold91<penguins$cr.age,0,1))
penguins$in.cr92=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold92<penguins$cr.age,0,1))
penguins$in.cr93=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold93<penguins$cr.age,0,1))
penguins$in.cr94=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold94<penguins$cr.age,0,1))
penguins$in.cr95=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold95<penguins$cr.age,0,1))
penguins$in.cr96=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold96<penguins$cr.age,0,1))
penguins$in.cr97=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold97<penguins$cr.age,0,1))
penguins$in.cr98=ifelse(penguins$did.cr==0, 0, ifelse(penguins$dayold98<penguins$cr.age,0,1))
 

penguins$never.cr=1-penguins$did.cr
penguins$did.cr=as.numeric(penguins$did.cr)
penguins$never.cr=as.numeric(penguins$never.cr)


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







