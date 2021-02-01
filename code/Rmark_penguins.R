#need to tell R that Mark is not stored in C:\Program Files

#personal computer
#MarkPath="C:/Users/Scott/Mark/"


library(RMark)
library(tidyverse)

setwd("C:/Users/Scott/Documents/THESIS/Data")


# the following script is adapted for my data from appendix C of the Mark book

# CJS analysis of penguin chick survival
# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/
penguins=convert.inp("resighting_survival/RMark_analysis/all_covs_cr_groups.inp", 
					group.df=data.frame(SEASON=c(rep("1213",2),rep("1314",2)),creched=rep(c("Y","N"),2)), 
					covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56", 
					 "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64", 
					 "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72",
					 "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80",
					 "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", 
					 "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96",
					 "dayold97", "dayold98", "sex", "av.food", "av.trip.length", "cr.age", "fail.age", "weight.slope40", 
					"flipper.slope40", "tibiotar.slope35", "CH_B", "CH_S", "mean.N", "cr.mass", 
					"res.htch", "cr.flip", "cr.tib", "did.cr", "fate.age" ),use.comments=FALSE) %>% 
  rename(res.htch = res.htch)

		 
##	setting NA values to mean of that field					

penguins$cr.mass[is.na(penguins$cr.mass)]=0		#1357.6				
penguins$cr.tib[is.na(penguins$cr.tib)]=0		#111.5				
penguins$cr.flip[is.na(penguins$cr.flip)]=0		#114.13				
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
 
#penguins$did.cr=as.numeric(penguins$did.cr)

penguins$never.cr=1-penguins$did.cr
penguins$did.cr=as.numeric(penguins$did.cr)
penguins$never.cr=as.numeric(penguins$never.cr)

# Next create the processed dataframe and the design data. We’ll use a group 
# variable for colony so it can be used in the set of models for Phi. Factor 
# variables (covariates with a small finite set of values) are best handled by using 
# them to define groups in the data. 
# need to set begin.time to 51 because RMark considers the first dayold to be the first release dayold, but really I want to consider recaptures beginning on dayold 52.
penguins.process=process.data(penguins,model="CJS",groups=c("SEASON", "creched"),begin.time=50) 

penguins.ddl=make.design.data(penguins.process)



###########3 RUN TO HERE TO START



# to view the design matrix for a given model 
write.csv(penguins.p.results[[1]]$design.matrix, "modtimes_dm.csv")

C:/Users/Scott/Documents/THESIS/Data/resighting_survival

PIMS(penguins.p.results[[10]], "p", simplified=FALSE)



penguins.p.results[[4]]$results$beta

#########################################
##	making models without mark.wrapper
##  I think this method of building and running models will be easier because of the npar adjustments I have to make
##  using mark.wrapper it is hard to keep track of which model is which, 
##  and I don't really have that many models to build.



run.model<-function(phi.stru, p.stru) {
mark(penguins.process,penguins.ddl,model.parameters=list(Phi=phi.stru,p=p.stru), output=FALSE)
}
peng.phi.gen.p.gen=run.model(list(formula=~SEASON*sex*time),list(formula=~SEASON*sex*time))
peng.phi.gen.p.SEASON=run.model(list(formula=~SEASON*sex*time),list(formula=~SEASON))
peng.phi.gen.p.sex=run.model(list(formula=~SEASON*sex*time),list(formula=~sex))
peng.phi.gen.p.time=run.model(list(formula=~SEASON*sex*time),list(formula=~time))
peng.phi.gen.p.SEASON.pl.sex=run.model(list(formula=~SEASON*sex*time),list(formula=~SEASON+sex) 	)
peng.phi.gen.p.SEASON.pl.time=run.model(list(formula=~SEASON*sex*time),list(formula=~SEASON+time) 	)
peng.phi.gen.p.sex.pl.time=run.model(list(formula=~SEASON*sex*time),list(formula=~sex+time))
peng.phi.gen.p.dot=run.model(list(formula=~SEASON*sex*time),list(formula=~1))
peng.phi.gen.p.T=run.model(list(formula=~SEASON*sex*time),list(formula=~(Time+1)))
peng.phi.gen.p.TT=run.model(list(formula=~SEASON*sex*time),list(formula=~(Time+1)+I((Time+1)^2)))
peng.phi.gen.p.lnT=run.model(list(formula=~SEASON*sex*time),list(formula=~log(Time+1)))
peng.phi.gen.p.cr=run.model(list(formula=~SEASON*sex*time),list(formula=~cr))

##  now collect models into one object to make a model comparison table
peng.p=collect.models()
model.table(peng.p, use.lnl=TRUE) 

## real and beta estimates
##	!! need to make sure xls file is created first
library(xlsx)
write.p.est<-function (mod, mod.name) {
mod.beta=as.data.frame(mod$results$beta)
write.csv(mod.beta, "temp/mod.beta.csv")
mod.beta.b=read.csv("temp/mod.beta.csv")
file.remove("temp/mod.beta.csv")
mod.beta.b$parm.type="beta"
colnames(mod.beta.b)[1] <- 'parm'

mod.real=as.data.frame(mod$results$real)
write.csv(mod.real, "temp/mod.real.csv")
mod.real.b=read.csv("temp/mod.real.csv")
file.remove("temp/mod.real.csv")
mod.real.b$parm.type="real"
colnames(mod.real.b)[1] <- 'parm'
mod.real.a=subset(mod.real.b, select=c("parm", "estimate", "se", "lcl", "ucl", "parm.type"))

mod.est=rbind(mod.beta.b, mod.real.a)

write.xlsx(mod.est, "figs_summaries/CH3/bydate/bydate_p_est.xlsx", sheetName=mod.name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
}						

write.p.est(peng.phi.gen.p.gen, "peng.phi.gen.p.gen")
write.p.est(peng.phi.gen.p.SEASON, "peng.phi.gen.p.SEASON")
write.p.est(peng.phi.gen.p.sex, "peng.phi.gen.p.sex")
write.p.est(peng.phi.gen.p.time, "peng.phi.gen.p.time")
write.p.est(peng.phi.gen.p.SEASON.pl.sex, "peng.phi.gen.p.SEASON.pl.sex")
write.p.est(peng.phi.gen.p.SEASON.pl.time, "peng.phi.gen.p.SEASON.pl.time")
write.p.est(peng.phi.gen.p.sex.pl.time, "peng.phi.gen.p.sex.pl.time")
write.p.est(peng.phi.gen.p.dot, "peng.phi.gen.p.dot")
write.p.est(peng.phi.gen.p.T, "peng.phi.gen.p.T")
write.p.est(peng.phi.gen.p.TT, "peng.phi.gen.p.TT")
write.p.est(peng.phi.gen.p.lnT, "peng.phi.gen.p.lnT")
write.p.est(peng.phi.gen.p.cr, "peng.phi.gen.p.cr")

##--
write.csv(peng.phi.gen.p.gen$design.matrix, "figs_summaries/CH3/bydate/temp/p.gen_dm_p.csv")
write.csv(peng.phi.gen.p.SEASON$design.matrix, "figs_summaries/CH3/bydate/temp/p.season_dm.csv")
write.csv(peng.phi.gen.p.sex$design.matrix, "figs_summaries/CH3/bydate/temp/p.sex_dm.csv")
write.csv(peng.phi.gen.p.time$design.matrix, "figs_summaries/CH3/bydate/temp/p.time_dm.csv")
write.csv(peng.phi.gen.p.SEASON.pl.sex$design.matrix, "figs_summaries/CH3/bydate/temp/p.SEASsex_dm.csv")
write.csv(peng.phi.gen.p.SEASON.pl.time$design.matrix, "figs_summaries/CH3/bydate/temp/p.SEAS.time.csv")
write.csv(peng.phi.gen.p.sex.pl.time$design.matrix, "figs_summaries/CH3/bydate/temp/p.sex.time.csv")
write.csv(peng.phi.gen.p.dot$design.matrix, "figs_summaries/CH3/bydate/temp/p.dot.csv")
write.csv(peng.phi.gen.p.T$design.matrix, "figs_summaries/CH3/bydate/temp/p.T.csv")
write.csv(peng.phi.gen.p.TT$design.matrix, "figs_summaries/CH3/bydate/temp/p.TT.csv")
write.csv(peng.phi.gen.p.lnT$design.matrix, "figs_summaries/CH3/bydate/temp/p.lnT.csv")
write.csv(peng.phi.gen.p.cr$design.matrix, "figs_summaries/CH3/bydate/temp/p.cr.csv")


	
################################
###	now moving on to modeling Phi, with the best p structure from above, which is SEASON+time

run.model<-function(phi.stru, p.stru) {
mark(penguins.process,penguins.ddl,model.parameters=list(Phi=phi.stru,p=list(formula=~SEASON+time)), output=FALSE)
}

##	first step will be to determine effect of SEASON and sex
Phi.SEASON.x.sex=run.model(list(formula=~SEASON*sex))
Phi.SEASON_sex=run.model(list(formula=~SEASON+sex))
Phi.SEASON=run.model(list(formula=~SEASON))
Phi.sex=run.model(list(formula=~sex))
Phi.dot=run.model(list(formula=~1))


peng.phi=collect.models()
model.table(peng.phi, use.lnl=TRUE) 

Phi.SEASON.x.sex$design.matrix

write.csv(Phi.cr.x.age_lnT$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.cr.x.age_lnT.csv")
write.csv(Phi.SEASON_sex$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.SEASON_sex.csv")
write.csv(Phi.SEASON$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.SEASON.csv")
write.csv(Phi.sex$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.sex.csv")

##	dot was best

##--------
##	then take best structure from step 1 and investigate time effects

Phi.t_y.s=run.model(list(formula = ~time+SEASON*sex))
Phi.T_y.s =run.model(list(formula=~Time+SEASON*sex))
Phi.TT_y.s=run.model(list(formula=~Time+I(Time^2)+SEASON*sex))
Phi.lnT_y.s=run.model(list(formula=~log(Time+1)+SEASON*sex))
peng.phi=collect.models(lx=c("Phi.t_y.s", "Phi.T_y.s", "Phi.TT_y.s", "Phi.lnT_y.s"))
model.table(peng.phi, use.lnl=TRUE) 



##	now run the non-creche age/size covariates of actual interest with lnT
Phi.lnT_y.s=run.model(list(formula=~log(Time+1)+SEASON*sex))		#keep for comparison
Phi.res.htch_lnT_y.s=run.model(list(formula=~res.htch+log(Time+1)+SEASON*sex))		#relative hatch day- test for importance of synchrony in season
Phi.in.cr_lnT_y.s=run.model(list(formula=~in.cr+log(Time+1)+SEASON*sex))		#in.cr=time-varying for in cr or not on specific day
Phi.in.cr_age_lnT_y.s=run.model(list(formula=~in.cr+dayold++SEASON*sex))		#dayold=time-varying for age on specific day; I think dayold and lnT can't be together
#Phi.in.cr.x.age_lnT_y.s=run.model(list(formula=~in.cr*dayold+log(Time+1)+SEASON*sex))
peng.phi=collect.models(lx=c("Phi.lnT_y.s", "Phi.res.htch_lnT_y.s", "Phi.in.cr_lnT_y.s", "Phi.in.cr_age_lnT_y.s"))
model.table(peng.phi, use.lnl=TRUE) 


# growth rates				 	
Phi.massgr=run.model(list(formula=~weight.slope40), list(formula=~SEASON+time))
Phi.flipgr=run.model(list(formula=~flipper.slope40), list(formula=~SEASON+time))
Phi.tibgr=run.model(list(formula=~tibiotar.slope35), list(formula=~SEASON+time))
Phi.massgr_in.cr_lnT_y.s=run.model(list(formula=~weight.slope40+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.flipgr_in.cr_lnT_y.s=run.model(list(formula=~flipper.slope40+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.tibgr_in.cr_lnT_y.s=run.model(list(formula=~tibiotar.slope35+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.in.cr_lnT_y.s=run.model(list(formula=~in.cr+log(Time+1)+SEASON*sex))
peng.phi=collect.models(lx=c("Phi.massgr", "Phi.flipgr", "Phi.tibgr",
							"Phi.massgr_in.cr_lnT_y.s", "Phi.flipgr_in.cr_lnT_y.s", "Phi.tibgr_in.cr_lnT_y.s", 
							"Phi.in.cr_lnT_y.s"))
model.table(peng.phi, use.lnl=TRUE) 



#creche age/sizes	
Phi.in.cr_lnT_y.s=run.model(list(formula=~in.cr+log(Time+1)+SEASON*sex))
Phi.tibgr_in.cr_lnT_y.s=run.model(list(formula=~tibiotar.slope35+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.dot=run.model(list(formula=~1))
				
Phi.cr.mass=run.model(list(formula=~creched+creched:cr.mass), list(formula=~SEASON+time))
Phi.cr.flip=run.model(list(formula=~creched+did.cr:cr.flip), list(formula=~SEASON+time))
Phi.cr.tib=run.model(list(formula=~creched+did.cr:cr.tib), list(formula=~SEASON+time))
Phi.fate.age=run.model(list(formula=~creched+did.cr:fate.age), list(formula=~SEASON+time))
					
Phi.fate.age_cr.mass=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.mass), list(formula=~SEASON+time))
Phi.fate.age_cr.flip=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.flip), list(formula=~SEASON+time))
Phi.fate.age_cr.tib=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.tib), list(formula=~SEASON+time))

Phi.cr.mass_lnT_y.s=run.model(list(formula=~creched+did.cr:cr.mass+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.cr.flip_lnT_y.s=run.model(list(formula=~creched+did.cr:cr.flip+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.cr.tib_lnT_y.s=run.model(list(formula=~creched+did.cr:cr.tib+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.fate.age_lnT_y.s=run.model(list(formula=~creched+did.cr:fate.age+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
					
Phi.fate.age_cr.mass_lnT_y.s=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.mass+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.fate.age_cr.flip_lnT_y.s=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.flip+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))
Phi.fate.age_cr.tib_lnT_y.s=run.model(list(formula=~creched+did.cr:fate.age+did.cr:cr.tib+in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))

peng.phi=collect.models(lx=c("Phi.in.cr_lnT_y.s", "Phi.tibgr_in.cr_lnT_y.s", "Phi.dot",
							"Phi.cr.mass", "Phi.cr.flip", "Phi.cr.tib", "Phi.fate.age", 					
							"Phi.fate.age_cr.mass", "Phi.fate.age_cr.flip", "Phi.fate.age_cr.tib", 
							"Phi.cr.mass_lnT_y.s", "Phi.cr.flip_lnT_y.s", "Phi.cr.tib_lnT_y.s", "Phi.fate.age_lnT_y.s",
							"Phi.fate.age_cr.mass_lnT_y.s", "Phi.fate.age_cr.flip_lnT_y.s", "Phi.fate.age_cr.tib_lnT_y.s"))
model.table(peng.phi, use.lnl=TRUE) 



library(xlsx)
write.phi.est<-function (mod, mod.name) {
mod.beta=as.data.frame(mod$results$beta)
write.csv(mod.beta, "temp/mod.beta.csv")
mod.beta.b=read.csv("temp/mod.beta.csv")
file.remove("temp/mod.beta.csv")
mod.beta.b$parm.type="beta"
colnames(mod.beta.b)[1] <- 'parm'

mod.real=as.data.frame(mod$results$real)
write.csv(mod.real, "temp/mod.real.csv")
mod.real.b=read.csv("temp/mod.real.csv")
file.remove("temp/mod.real.csv")
mod.real.b$parm.type="real"
colnames(mod.real.b)[1] <- 'parm'
mod.real.a=subset(mod.real.b, select=c("parm", "estimate", "se", "lcl", "ucl", "parm.type"))

mod.est=rbind(mod.beta.b, mod.real.a)

write.xlsx(mod.est, "figs_summaries/CH3/cr_groups/Phi_est.xlsx", sheetName=mod.name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
}	
##	write beta and real estimates for the models that make up >=90% of model weight, plus dot model
##	for models with cr age and sizes=0 for chicks that didn't creche
write.phi.est(Phi.cr.mass_lnT_y.s, "Phi.cr.mass_lnT_y.s")
write.phi.est(Phi.cr.flip_lnT_y.s, "Phi.cr.flip_lnT_y.s")
write.phi.est(Phi.cr.tib_lnT_y.s, "Phi.cr.tib_lnT_y.s")
write.phi.est(Phi.fate.age_cr.flip_lnT_y.s, "Phi.fate.age_cr.flip_lnT_y.s")
write.phi.est(Phi.fate.age_cr.mass_lnT_y.s, "Phi.fate.age_cr.mass_lnT_y.s")
write.phi.est(Phi.fate.age_cr.tib_lnT_y.s, "Phi.fate.age_cr.tib_lnT_y.s")
write.phi.est(Phi.dot, "Phi.dot")


##-----
write.csv(Phi.cr.mass_lnT_y.s$design.matrix, "figs_summaries/CH3/cr_groups/Phi.cr.mass_lnT_y.s_DM.csv")

write.csv(Phi.sex_fate.age_cr.mass$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.sex_fate.age_cr.mass.csv")
write.csv(Phi.SEASON_fate.age_cr.mass$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.SEASON_fate.age_cr.mass.csv")
write.csv(Phi.dot$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.dot.csv")



## then build and run models, and make parameter adjustment							
peng.phi.gen.p.SEASON.pl.time=adjust.parameter.count(peng.phi.gen.p.SEASON.pl.time,231)
peng.phi.SEASON.p.SEASON.pl.time=adjust.parameter.count(peng.phi.SEASON.p.SEASON.pl.time,49)
peng.phi.sex.p.SEASON.pl.time=adjust.parameter.count(peng.phi.sex.p.SEASON.pl.time,49)
peng.phi.time.p.SEASON.pl.time=adjust.parameter.count(peng.phi.time.p.SEASON.pl.time,93)
peng.phi.SEASON.pl.sex.p.SEASON.pl.time=adjust.parameter.count(peng.phi.SEASON.pl.sex.p.SEASON.pl.time,50)
peng.phi.SEASON.pl.time.p.SEASON.pl.time=adjust.parameter.count(peng.phi.SEASON.pl.time.p.SEASON.pl.time,94)
peng.phi.sex.pl.time.p.SEASON.pl.time=adjust.parameter.count(peng.phi.sex.pl.time.p.SEASON.pl.time,94)
peng.phi.dot.p.SEASON.pl.time=adjust.parameter.count(peng.phi.dot.p.SEASON.pl.time,48)
peng.phi.T.p.SEASON.pl.time=adjust.parameter.count(peng.phi.T.p.SEASON.pl.time,49)
peng.phi.TT.p.SEASON.pl.time=adjust.parameter.count(peng.phi.TT.p.SEASON.pl.time,50)
peng.phi.lnT.p.SEASON.pl.time=adjust.parameter.count(peng.phi.lnT.p.SEASON.pl.time,50)
peng.phi.times.p.SEASON.pl.time=adjust.parameter.count(peng.phi.times.p.SEASON.pl.time,97)



peng.phi=collect.models()
model.table(peng.phi, use.lnl=TRUE) 


write.csv(peng.phi.T.p.SEASON.pl.time$design.matrix, "phi.T.p.SEASON.pl.time.csv")



## plotting estimates from some top models
library(ggplot2)
# 
estimates=Phi.cr.mass_lnT_y.s$results$real
parm=factor(c (rep("Phi", 192), rep("p", 96)))
group=factor(c	(rep("Phi1213.notcr", 48),rep("Phi1314.notcr", 48),
				rep("Phi1213.cr", 48),rep("Phi1314.cr", 48),
				rep("p1213", 48),rep("p1314", 48)))
dayold=factor(c(rep(paste("dayold", 50:97, sep=""),6)))
estimates.a=as.data.frame(cbind(estimates, dayold, parm, group))
estimates=subset(estimates.a, parm=="Phi")
estimates=droplevels(estimates)
ucl= estimates$ucl 
lcl= estimates$lcl 
limits <- aes(ymax = ucl , ymin=lcl)
ggplot(estimates, aes(dayold, estimate)) + geom_point(aes(color=group))+
										#geom_errorbar(limits, width=0.25) + 
										labs(title = "Estimates")
										
ggsave("figs_summaries/CH3/cr_groups/Phi.cr.mass_lnT_y.s_plot.jpg")
	
# Phi(lnT)p(SEASON+time)
estimates=peng.phi.lnT.p.SEASON.pl.time$results$real
parm=factor(c(rep("Phi", 46),rep("p1213", 46),rep("p1314", 46)))
dayold=factor(c(rep(paste("dayold", 52:97, sep=""),3)))
estimates.a=as.data.frame(cbind(estimates, dayold, parm))
ucl= estimates.a$ucl 
lcl= estimates.a$lcl 
limits <- aes(ymax = ucl , ymin=lcl)
ggplot(estimates.a, aes(dayold, estimate)) + geom_point()+
										geom_errorbar(limits, width=0.25) + 
										facet_grid(parm ~ .)+ 
										labs(title = "Estimates for Phi(lnT)p(SEASON+time)")
ggsave("C:/Users/Scott/Documents/THESIS/Data/resighting_survival/Phi(lnT)p(SEASON+time).pdf")


# Phi(TT)p(SEASON+time)
estimates=peng.phi.TT.p.SEASON.pl.time$results$real
parm=factor(c(rep("Phi", 46),rep("p1213", 46),rep("p1314", 46)))
dayold=factor(c(rep(paste("dayold", 52:97, sep=""),3)))
estimates.a=as.data.frame(cbind(estimates, dayold, parm))
ucl= estimates.a$ucl 
lcl= estimates.a$lcl 
limits <- aes(ymax = ucl , ymin=lcl)
ggplot(estimates.a, aes(dayold, estimate)) + geom_point()+
										geom_errorbar(limits, width=0.25) + 
										facet_grid(parm ~ .)+ 
										labs(title = "Estimates for Phi(TT)p(SEASON+time)")
ggsave("C:/Users/Scott/Documents/THESIS/Data/resighting_survival/Phi(TT)p(SEASON+time).pdf")



# Phi(times)p(SEASON+time)
estimates=peng.phi.times.p.SEASON.pl.time$results$real
parm=factor(c(rep("Phi", 46),rep("p1213", 46),rep("p1314", 46)))
dayold=factor(c(rep(paste("dayold", 52:97, sep=""),3)))
estimates.a=as.data.frame(cbind(estimates, dayold, parm))
ucl= estimates.a$ucl 
lcl= estimates.a$lcl 
limits <- aes(ymax = ucl , ymin=lcl)
ggplot(estimates.a, aes(dayold, estimate)) + geom_point()+
										geom_errorbar(limits, width=0.25) + 
										facet_grid(parm ~ .)+ 
										labs(title = "Phi(times)p(SEASON+time)")
ggsave("C:/Users/Scott/Documents/THESIS/Data/resighting_survival/Phi(TT)p(SEASON+time).pdf")