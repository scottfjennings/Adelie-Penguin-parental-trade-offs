#need to tell R that Mark is not stored in C:\Program Files

#personal computer
#MarkPath="C:/Users/Scott/Mark/"
library(RMark)

setwd("C:/Users/Scott/Documents/THESIS/Data")
options(scipen = 999)

##		add required packages
library(lme4) 
library(AICcmodavg)
library(grid)
library(ggplot2)
library(plyr)
library(tidyverse)


# the following script is adapted for my data from appendix C of the Mark book

# CJS analysis of penguin chick survival
# Import data (all_covs_groups.inp) and convert it from the MARK inp file format to the \textbf{RMark}
# format using the function convert.inp 
# It is defined with 4 groups: Males in 1213, females in 1213, males in 1314, and females in 1314
# This structure is defined with the group.df argument of convert.inp. 
##C:/Users/jenninsc/Documents/THESIS/Data/resighting_survival/RMark_analysis/
penguins=convert.inp("resighting_survival/RMark_analysis/all_covs_sex_groups.inp", 
					group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), 
					covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56", 
					 "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64", 
					 "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72",
					 "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80",
					 "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", 
					 "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96",
					 "dayold97", "dayold98", "av.food", "av.trip.length", "cr.age", "fail.age", "weight.slope40", 
					"flipper.slope40", "tibiotar.slope35", "CH_B", "CH_S", "mean.N", "cr.mass", 
					"res.htch", "cr.flip", "cr.tib", "did.cr", "fate.age" ),use.comments=FALSE) %>% 
  rename(res.htch = res.htch)

		 
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
penguins.process=process.data(penguins,model="CJS",groups=c("SEASON", "sex"),begin.time=50) 
penguins.ddl=make.design.data(penguins.process)



export.chdata(penguins.process, filename="penguin", replace=TRUE)

## GOF testing
peng.for.gof=subset(penguins, select=c("ch", "freq", "sex", "SEASON"))

count(peng.for.gof, "ch")


peng.for.gof.proc=process.data(peng.for.gof,model="CJS",groups=c("SEASON", "sex"))

penguins.GOF=release.gof(penguins.process)

#########################################
##	making models without mark.wrapper
##  I think this method of building and running models will be easier because of the npar adjustments I have to make
##  using mark.wrapper it is hard to keep track of which model is which, 
##  and I don't really have that many models to build.


#Determine p structure
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
peng.p.adj=adjust.chat(1.24, peng.p)
p_table <- model.table(peng.p.adj, use.lnl=TRUE) 
write.csv(p_table, "C:/Users/Scott/Documents/THESIS/Thesis/CH4/MSdocs/p_table.csv", row.names = F)
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

write.xlsx(mod.est, "figs_summaries/CH3/sex_groups/p_est.xlsx", sheetName=mod.name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
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
mark(penguins.process,penguins.ddl,model.parameters=list(Phi=phi.stru,p=list(formula=~SEASON+time)), output=FALSE, chat=1.24)
}

##	first step will be to determine effect of SEASON and sex
Phi.SEASON.x.sex=run.model(list(formula=~SEASON*sex))
Phi.SEASON_sex=run.model(list(formula=~SEASON+sex))
Phi.SEASON=run.model(list(formula=~SEASON))
Phi.sex=run.model(list(formula=~sex))
Phi.dot=run.model(list(formula=~1))


step1=collect.models(lx=c("Phi.SEASON.x.sex", "Phi.SEASON_sex", "Phi.SEASON", "Phi.sex", "Phi.dot"))
phi_table1 <- model.table(step1, use.lnl=TRUE) 

write.csv(phi_table1, "C:/Users/Scott/Documents/THESIS/Thesis/CH4/MSdocs/phi_table1.csv", row.names = F)
Phi.SEASON.x.sex$design.matrix
##---

export.chdata(peng.for.gof.proc, filename="penguin.gof", replace=TRUE)

export.model(step2)


write.csv(Phi.cr.x.age_lnT$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.cr.x.age_lnT.csv")
write.csv(Phi.SEASON_sex$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.SEASON_sex.csv")
write.csv(Phi.SEASON$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.SEASON.csv")
write.csv(Phi.sex$design.matrix, "figs_summaries/CH3/bydate/temp/Phi.sex.csv")

## Step 2	is relative hatch date important?
Phi.y.s=run.model(list(formula=~SEASON*sex))
Phi.res.htch_y.s=run.model(list(formula=~res.htch+SEASON*sex))
Phi.res.htch.y.s=run.model(list(formula=~res.htch*SEASON*sex))
Phi.res.htch2_y.s=run.model(list(formula=~res.htch+I(res.htch^2)+SEASON*sex))

step2=collect.models(lx=c("Phi.y.s", "Phi.res.htch_y.s", "Phi.res.htch.y.s", "Phi.res.htch2_y.s"))
model.table(step2, use.lnl=TRUE) 


Phi.SEASON.x.sex = run.model(list(formula = ~SEASON * sex))
Phi.SEASON_sex = run.model(list(formula = ~SEASON + sex))
Phi.SEASON = run.model(list(formula = ~SEASON))
Phi.sex = run.model(list(formula = ~sex))
Phi.dot = run.model(list(formula = ~1))
Phi.y.s = run.model(list(formula = ~SEASON * sex))

Phi.res.htch_y.s = run.model(list(formula = ~res.htch + SEASON * sex))
Phi.res.htch.y.s = run.model(list(formula = ~res.htch * SEASON * sex))
Phi.res.htch2_y.s = run.model(list(formula = ~res.htch + I(res.htch^2) + SEASON * sex))

Phi.res.htch_y_s = run.model(list(formula = ~res.htch + SEASON + sex))
Phi.res.htch.y_s = run.model(list(formula = ~res.htch * SEASON + sex))
Phi.res.htch2_y_s = run.model(list(formula = ~res.htch + I(res.htch^2) + SEASON + sex))

Phi.res.htch_y = run.model(list(formula = ~res.htch + SEASON))
Phi.res.htch.y = run.model(list(formula = ~res.htch * SEASON))
Phi.res.htch2_y = run.model(list(formula = ~res.htch + I(res.htch^2) + SEASON))

Phi.res.htch_s = run.model(list(formula = ~res.htch + sex))
Phi.res.htch.s = run.model(list(formula = ~res.htch * sex))
Phi.res.htch2_s = run.model(list(formula = ~res.htch + I(res.htch^2) + sex))

Phi.res.htch = run.model(list(formula = ~res.htch))
Phi.res.htch2 = run.model(list(formula = ~res.htch + I(res.htch^2)))


step2.1=collect.models(lx=c("Phi.SEASON.x.sex", "Phi.SEASON_sex", "Phi.SEASON", "Phi.sex", "Phi.dot", "Phi.y.s", 
                          "Phi.res.htch_y.s", "Phi.res.htch.y.s", "Phi.res.htch2_y.s",
                          "Phi.res.htch_y_s", "Phi.res.htch.y_s", "Phi.res.htch2_y_s",
                          "Phi.res.htch_y", "Phi.res.htch.y", "Phi.res.htch2_y",
                          "Phi.res.htch_s", "Phi.res.htch.s", "Phi.res.htch2_s",
                          "Phi.res.htch", "Phi.res.htch2"))
phi_combined_base_stru <- model.table(step2.1, use.lnl=TRUE) 
write.csv(phi_combined_base_stru, "C:/Users/Scott/Documents/THESIS/Thesis/CH4/MSdocs/phi_combined_base_stru.csv", row.names = F)


##--------
##	Step 3: take best structure from above and investigate time effects
Phi.res.htch_y.s=run.model(list(formula=~res.htch+SEASON*sex))
Phi.t_hatch_y.s=run.model(list(formula=~time+res.htch+SEASON*sex))		#1
Phi.T_hatch_y.s =run.model(list(formula=~Time+res.htch+SEASON*sex))
Phi.TT_hatch_y.s=run.model(list(formula=~Time+I(Time^2)+res.htch+SEASON*sex))
Phi.lnT_hatch_y.s=run.model(list(formula=~log(Time+1)+res.htch+SEASON*sex))
##	including time constraints based on time-varying covariates
#in.cr=time-varying for in cr or not on specific day
#dayold=time-varying for age on specific day;
Phi.t_in.cr_hatch_y.s=run.model(list(formula=~time+in.cr+res.htch+SEASON*sex))		#5	#same intercept, diff slope for in.cr
Phi.t_age_hatch_y.s=run.model(list(formula=~time+dayold+res.htch+SEASON*sex))		#6	same intercept, diff slope for dayold
Phi.t.in.cr_hatch_y.s=run.model(list(formula=~time:in.cr+res.htch+SEASON*sex))		#same slope, diff intercept for in.cr
Phi.t.age_hatch_y.s=run.model(list(formula=~time:dayold+res.htch+SEASON*sex))		#8	same slope, diff intercept for dayold
	
step3=collect.models(lx=c("Phi.res.htch_y.s", "Phi.t_hatch_y.s", "Phi.T_hatch_y.s", "Phi.TT_hatch_y.s", "Phi.lnT_hatch_y.s", 
							"Phi.t_in.cr_hatch_y.s", "Phi.t_age_hatch_y.s", "Phi.t.in.cr_hatch_y.s", "Phi.t.age_hatch_y.s"))
model.table(step3, use.lnl=TRUE) 

step3.adj=adjust.chat(1.24, step3)
model.table(step3.adj, use.lnl=TRUE) 

adjust.chat(1.5, step3)
adjust.chat(2, step3)



write.csv(Phi.t_in.cr_hatch_y.s$design.matrix, "C:/Users/Scott/Documents/Classes/FW661/Phi.t_age_hatch_y.s.csv")
write.csv(Phi.t.in.cr_hatch_y.s$design.matrix, "C:/Users/Scott/Documents/Classes/FW661/Phi.t.in.cr_hatch_y.s.csv")
write.csv(Phi.cr.mass_T_hatch_y.s$design.matrix, "C:/Users/Scott/Documents/Classes/FW661/Phi.cr.mass_T_hatch_y.s.csv")


##--model averaging for step 3 models- this does real estimates
step3.mod.av=model.average(step3, "Phi", vcv=TRUE)

f1213.step3.est=subset(step3.mod.av$estimates[1:48,], select=c("time", "estimate", "se", "lcl", "ucl", "group", "SEASON", "sex"))
f1314.step3.est=subset(step3.mod.av$estimates[1177:1224,], select=c("time", "estimate", "se", "lcl", "ucl", "group", "SEASON", "sex"))
m1213.step3.est=subset(step3.mod.av$estimates[2353:2400,], select=c("time", "estimate", "se", "lcl", "ucl", "group", "SEASON", "sex"))
m1314.step3.est=subset(step3.mod.av$estimates[3529:3576,], select=c("time", "estimate", "se", "lcl", "ucl", "group", "SEASON", "sex"))
step3.est=rbind(f1213.step3.est, f1314.step3.est, m1213.step3.est, m1314.step3.est)

write.csv(step3.est, "figs_summaries/CH3/sex_groups/step3.est.chat1.24.csv")

head(f1213.step3.est)
################################################################
# Step 4 growth rates	
Phi.dot	 =	run.model(list(formula	 =	~1)	)			 	
Phi.massgr=run.model(list(formula=~weight.slope40), list(formula=~SEASON+time))
Phi.flipgr=run.model(list(formula=~flipper.slope40), list(formula=~SEASON+time))
Phi.tibgr=run.model(list(formula=~tibiotar.slope35), list(formula=~SEASON+time))
Phi.massgr_T_hatch_y.s=run.model(list(formula=~weight.slope40+Time+res.htch+SEASON*sex), list(formula=~SEASON+time))
Phi.flipgr_T_hatch_y.s=run.model(list(formula=~flipper.slope40+Time+res.htch+SEASON*sex), list(formula=~SEASON+time))
Phi.tibgr_T_hatch_y.s=run.model(list(formula=~tibiotar.slope35+Time+res.htch+SEASON*sex), list(formula=~SEASON+time))
Phi.T_y.s=run.model(list(formula=~Time+res.htch+SEASON*sex))
step4=collect.models(lx=c("Phi.massgr", "Phi.flipgr", "Phi.tibgr",
							"Phi.massgr_T_hatch_y.s", "Phi.flipgr_T_hatch_y.s", "Phi.tibgr_T_hatch_y.s", 
							"Phi.T_y.s", "Phi.dot"))
model.table(step4, use.lnl=TRUE) 

step4.chat1.24=adjust.chat(1.24, step4)
step4.chat1.5=adjust.chat(1.5, step4)
step4.chat2=adjust.chat(2, step4)

write.phi.est(Phi.massgr, "Phi.massgr")
write.phi.est(Phi.flipgr, "Phi.flipgr") 
write.phi.est(Phi.tibgr, "Phi.tibgr")
write.phi.est(Phi.massgr_T_hatch_y.s, "Phi.massgr_T_hatch_y.s") 
write.phi.est(Phi.flipgr_T_hatch_y.s, "Phi.flipgr_T_hatch_y.s") 
write.phi.est(Phi.tibgr_T_hatch_y.s, "Phi.tibgr_T_hatch_y.s") 
write.phi.est(Phi.T_y.s, "Phi.T_y.s")

write.xlsx(summary.mark(Phi.tibgr_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step4_Phi_est.chat1.5.xlsx", sheetName="Phi.tibgr_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
write.xlsx(summary.mark(Phi.massgr_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step4_Phi_est.chat1.5.xlsx", sheetName="Phi.massgr_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			


#creche age/sizes	
Phi.in.cr_y.s	 =	run.model(list(formula	 =	~Time+res.htch+SEASON*sex)	)
Phi.tibgr_in.cr.t_y.s	 =	run.model(list(formula	 =	~tibiotar.slope35+Time+res.htch+SEASON*sex)	)
Phi.dot	 =	run.model(list(formula	 =	~1)	)
					
Phi.cr.mass	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.mass)	)
Phi.cr.flip	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.flip)	)
Phi.cr.tib	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.tib)	)
Phi.fate.age	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age)	)
					
Phi.fate.age_cr.mass	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.mass)	)
Phi.fate.age_cr.flip	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.flip)	)
Phi.fate.age_cr.tib	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.tib)	)
					
Phi.cr.mass_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.mass+Time+res.htch+SEASON*sex)	)
Phi.cr.flip_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.flip+Time+res.htch+SEASON*sex)	)
Phi.cr.tib_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.tib+Time+res.htch+SEASON*sex)	)
Phi.fate.age_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+Time+res.htch+SEASON*sex)	)
					
Phi.fate.age_cr.mass_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.mass+Time+res.htch+SEASON*sex)	)
Phi.fate.age_cr.flip_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.flip+Time+res.htch+SEASON*sex)	)
Phi.fate.age_cr.tib_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.tib+Time+res.htch+SEASON*sex)	)
					
Phi.cr.mass_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.mass+Time+res.htch+SEASON*sex)	)
Phi.cr.flip_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.flip+Time+res.htch+SEASON*sex)	)
Phi.cr.tib_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.tib+Time+res.htch+SEASON*sex)	)
Phi.fate.age_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:fate.age+Time+res.htch+SEASON*sex)	)

##-- take out T from best models
Phi.cr.mass_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.mass+res.htch+SEASON*sex)	)
Phi.cr.flip_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.flip+res.htch+SEASON*sex)	)
Phi.cr.tib_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.tib+res.htch+SEASON*sex)	)


step5=collect.models(lx=c("Phi.in.cr_y.s", "Phi.tibgr_in.cr.t_y.s", "Phi.dot",
							"Phi.cr.mass", "Phi.cr.flip", "Phi.cr.tib", "Phi.fate.age",
							"Phi.fate.age_cr.mass", "Phi.fate.age_cr.flip", "Phi.fate.age_cr.tib",
							"Phi.cr.mass_T_hatch_y.s", "Phi.cr.flip_T_hatch_y.s", "Phi.cr.tib_T_hatch_y.s", "Phi.fate.age_T_hatch_y.s",
							"Phi.fate.age_cr.mass_T_hatch_y.s", "Phi.fate.age_cr.flip_T_hatch_y.s", "Phi.fate.age_cr.tib_T_hatch_y.s",
							"Phi.cr.mass_hatch_y.s", "Phi.cr.flip_hatch_y.s", "Phi.cr.tib_hatch_y.s"))
model.table(step5, use.lnl=TRUE) 

adjust.chat(1.24, step5)
adjust.chat(1.5, step5)
adjust.chat(2, step5)

test=run.model(list(formula=~did.cr+did.cr:cr.mass+time:in.cr+log(Time+1)+SEASON*sex), list(formula=~SEASON+time))


best_real <- compute.real(Phi.cr.flip_T_hatch_y.s, data = penguins) %>% 
  mutate(estimate = round(estimate, 5),
         se = round(se, 5),
         lcl = round(lcl, 5),
         ucl = round(ucl, 5))
best_real2 <- best_real[1:192,]
best_real2 <- best_real2 %>% 
  mutate(group = c(rep("1213_F", 48), rep("1314_F", 48), rep("1213_M", 48), rep("1314_M", 48)))

write.xlsx(summary.mark(Phi.cr.mass_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step5_Phi_est.chat1.5.xlsx", sheetName="Phi.cr.mass_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
write.xlsx(summary.mark(Phi.cr.flip_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step5_Phi_est.chat1.5.xlsx", sheetName="Phi.cr.flip_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
write.xlsx(summary.mark(Phi.cr.tib_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step5_Phi_est.chat1.24.xlsx", sheetName="Phi.cr.tib_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
write.xlsx(summary.mark(Phi.fate.age_cr.flip_T_hatch_y.s)$beta, "figs_summaries/CH3/sex_groups/step5_Phi_est.chat1.24.xlsx", sheetName="Phi.fate.age_cr.flip_T_hatch_y.s", col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
##--
Phi.cr.mass_T_hatch_y.s	 =	run.model(list(formula	 =	~did.cr+did.cr:cr.mass+Time+res.htch+SEASON*sex)	)

pred1000=covariate.predictions(Phi.cr.mass_T_hatch_y.s, data=data.frame(cr.mass=1000), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
#pred1500=covariate.predictions(Phi.cr.mass_T_hatch_y.s, data=data.frame(cr.mass=1500), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
pred2000=covariate.predictions(Phi.cr.mass_T_hatch_y.s, data=data.frame(cr.mass=2000), indices=c(1:48, 1177:1224, 2353:2400, 3529:3576), alpha=0.025)
pred10=as.data.frame(pred1000$estimates)
#pred15=as.data.frame(pred1500$estimates)
pred20=as.data.frame(pred2000$estimates)
pred10$pred=1000
#pred15$pred=1500
pred20$pred=2000
pred.all=rbind(pred10, pred20)
pred.all = within(pred.all, { 
sex=ifelse(model.index<97,"Female","Male")
SEASON=factor(c	(rep("2012", 48),rep("2013", 48),
				rep("2012", 48),rep("2013", 48),
				rep("2012", 48),rep("2013", 48),
				rep("2012", 48),rep("2013", 48)))
day=as.numeric(c(rep(50:97,8)))	
pred=as.factor(pred)			
})
write.csv(pred.all, "figs_summaries/CH3/sex_groups/step5.est.csv")
pred.all= read.csv("figs_summaries/CH3/sex_groups/step5.est.csv")
pred.all$day=as.numeric(pred.all$day)

ggplot(data=pred.all, aes(x=day, y=estimate))+
	geom_point(size=0)+
	geom_line(size=1, data=subset(pred.all, pred=="1000"), color="red") +
	geom_line(size=1, data=subset(pred.all, pred=="2000"), color="blue") +
	geom_line(aes(y=ucl), linetype="dotted",size=1, data=subset(pred.all, pred=="1000"), color="red") +
	geom_line(aes(y=lcl), linetype="dotted",size=1, data=subset(pred.all, pred=="1000"), color="red") +
	geom_line(aes(y=ucl), linetype="dotted",size=1, data=subset(pred.all, pred=="2000"), color="blue") +
	geom_line(aes(y=lcl), linetype="dotted",size=1, data=subset(pred.all, pred=="2000"), color="blue") +
	scale_colour_discrete(name="Creching\nmass (g)")+
	geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(sex~SEASON)+
	xlab("Day of Season")+
	ylim(0.65,1)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					#,strip.text = element_blank() 
					,strip.background = element_blank()
					#,plot.margin = unit( c(-1,1,0,0) , units = "lines" )
					,axis.ticks = element_blank()
					#,axis.text.x = element_blank()
					,axis.text = element_text(size=15)
					,axis.title = element_text(size=20)
					,axis.title.y = element_blank()
					) 
					



library(xlsx)
write.phi.est<-function (mod, mod.name) {
mod.beta=as.data.frame(summary.mark(mod)$beta)
write.csv(mod.beta, "temp/mod.beta.csv")
mod.beta.b=read.csv("temp/mod.beta.csv")
file.remove("temp/mod.beta.csv")
mod.beta.b$parm.type="beta"
colnames(mod.beta.b)[1] <- 'parm'

mod.real=as.data.frame(summary.mark(mod)$reals)
write.csv(mod.real, "temp/mod.real.csv")
mod.real.b=read.csv("temp/mod.real.csv")
file.remove("temp/mod.real.csv")
mod.real.b$parm.type="real"
colnames(mod.real.b)[1] <- 'parm'
mod.real.a=subset(mod.real.b, select=c("parm", "estimate", "se", "lcl", "ucl", "parm.type"))

mod.est=rbind(mod.beta.b, mod.real.a)

write.xlsx(mod.est, "figs_summaries/CH3/sex_groups/step4_Phi_est.chat1.24.xlsx", sheetName=mod.name, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)			
}	
##	write beta and real estimates for the models that make up >=90% of model weight, plus dot model
##	for models with cr age and sizes=0 for chicks that didn't creche






write.phi.est(Phi.cr.flip_T_hatch_y.s, "Phi.cr.flip_T_hatch_y.s")
write.phi.est(Phi.cr.mass_T_hatch_y.s, "Phi.cr.mass_T_hatch_y.s")
write.phi.est(Phi.cr.tib_T_hatch_y.s, "Phi.cr.tib_T_hatch_y.s")
write.phi.est(Phi.fate.age_cr.flip_T_hatch_y.s, "Phi.fate.age_cr.flip_T_hatch_y.s")
write.phi.est(Phi.dot, "Phi.dot")

write.phi.est(Phi.fate.age_cr.tib_lnT_y.s, "Phi.fate.age_cr.tib_lnT_y.s")


##-----

write.csv(Phi.cr.flip_T_hatch_y.s$design.matrix, "figs_summaries/CH3/sex_groups/Phi.cr.flip_T_hatch_y.s.DM.csv")
write.csv(Phi.fate.age_cr.flip_T_hatch_y.s$design.matrix, "figs_summaries/CH3/sex_groups/Phi.fate.age_cr.flip_T_hatch_y.s.DM.csv")



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

mod.av=read.csv("figs_summaries/CH3/sex_groups/step3.mod.av.csv")


## plotting estimates from some top models
library(ggplot2)
estimates=Phi.cr.flip_T_hatch_y.s$results$real
parm=factor(c (rep("Phi", 192), rep("p", 96)))
group=factor(c	(rep("F1213", 48),rep("F1314", 48),
				rep("M1213", 48),rep("M1314", 48),
				rep("p1213", 48),rep("p1314", 48)))
day=factor(c(rep(50:97,6)))
estimates.a=as.data.frame(cbind(estimates, day, parm, group))
estimates=subset(estimates.a, parm=="Phi")
estimates=droplevels(estimates)
ucl= step3.est$ucl 
lcl= step3.est$lcl 
limits <- aes(ymax = ucl , ymin=lcl)
qplot( time, estimate, data=step3.est, group=group) + 
										geom_point(aes(color=group))+
										geom_errorbar(limits, width=0.25) + 
										guides(shape=FALSE)+
										facet_wrap(~group, ncol=2)+
										labs(title = "Model-averaged Survival Estimates-grouped by year and sex; c-hat=1.24")
										
ggsave("figs_summaries/CH3/sex_groups/mod.av_plot.jpg")
	
step3.est=read.csv("figs_summaries/CH3/sex_groups/step3.mod.av.csv")	
step3.est

## now make plot for MS
library(gridExtra)
ucl= step3.est$ucl 
lcl= step3.est$lcl 
limits <- aes(ymax = ucl , ymin=lcl)			
a=qplot( time, estimate, data=subset(step3.est, SEASON=="1213"), shape=sex, group=sex,  ylab="Daily Survival Probability", xlab="")+
	geom_errorbar(limits, width=0.4)+
	geom_vline(xintercept = 67, linetype = "longdash")+
	geom_errorbar(limits, width=0.75, data=subset(step3.est, SEASON=="1213"&sex=="Male"))+
	geom_line(size=0.5)+
	ylim(0.7,1)+
	#scale_y_continuous(breaks=seq(0.7, 1, 0.025))+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,legend.position="none"
					) 

	
b=qplot( time, estimate, data=subset(step3.est, SEASON=="1314"), shape=sex, group=sex,  ylab="", xlab="")+
	geom_errorbar(limits, width=0.5)+
	geom_vline(xintercept = 67, linetype = "longdash")+
	geom_line(size=0.5)+
	geom_errorbar(limits, width=0.75, data=subset(step3.est, SEASON=="1314"&sex=="Male"))+
	ylim(0.7,1)+
	#scale_y_continuous(breaks=seq(0.7, 1, 0.025))+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,legend.position="none"
					) 



e=qplot(time, estimate, data=subset(step3.est, SEASON=="1314"), shape=sex, group=sex,  ylab="daily survival rate", xlab="day of season")+
	geom_errorbar(limits, width=0.25)+
	scale_shape(guide = guide_legend(title = NULL), labels=c("Female", "Male"))+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,plot.margin= unit(c(1, 1, -1, 1), "lines")
					)
					
g_legend<-function(p){
tmp <- ggplotGrob(p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)}

mylegend<-g_legend(e)


grid.arrange(a, b, mylegend,
                        main=textGrob("Model-averaged Survival Estimates", just="top", gp=gpar(fontsize=20)),
						sub=textGrob("Day of season", vjust = -1),
						ncol=3,widths = unit(c(4, 4, .5), "null"))
						
##---

step3.est=read.csv("figs_summaries/CH3/sex_groups/step3.mod.av.csv")	
step3.est$SEASON=as.factor(step3.est$SEASON)

step3.est$SEASON=ifelse(step3.est$SEASON=="1213", 2012, 2013)				


	
ggplot(data=step3.est, aes(x=time, y=estimate))+
	geom_point(size=0)+
	geom_line(size=1.5, data=subset(step3.est, sex=="Male"), color="red") +
	geom_line(size=1.5, data=subset(step3.est, sex=="Female"), color="blue") +
	geom_line(aes(y=ucl), linetype="dotted",size=1.5, data=subset(step3.est, sex=="Male"), color="red") +
	geom_line(aes(y=lcl), linetype="dotted",size=1.5, data=subset(step3.est, sex=="Male"), color="red") +
	geom_line(aes(y=ucl), linetype="dotted",size=1.5, data=subset(step3.est, sex=="Female"), color="blue") +
	geom_line(aes(y=lcl), linetype="dotted",size=1.5, data=subset(step3.est, sex=="Female"), color="blue") +
	#annotate("text", x = 75, y = 0.85, label = "Male", color="red", size=10)+
	#annotate("text", x = 77.5, y = 0.83, label = "Female", color="blue", size=10)+
	geom_vline(xintercept = 67, linetype = "longdash", alpha=0.3, size=1)+
	facet_grid(~SEASON)+
	xlab("Day of Season")+
	ylim(0.82,1)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					#,plot.margin = unit( c(-1,1,0,0) , units = "lines" )
					,axis.ticks = element_blank()
					,axis.text.x = element_blank()
					,axis.text = element_text(size=15)
					,axis.title = element_text(size=20)
					,axis.title.y = element_blank()
					) 
						