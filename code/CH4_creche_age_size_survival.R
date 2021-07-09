
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
						