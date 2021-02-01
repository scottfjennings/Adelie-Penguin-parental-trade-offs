


library(lme4)
library(AICcmodavg)
library(grid)
library(ggplot2)
library(plyr)
library(scatterplot3d)

setwd("C:/Users/Scott/Documents/THESIS/Data")


ana.table=read.csv("chick_growth/ana.table.full.csv")
#but need to turn SEASON back into a factor
ana.table$SEASON<-as.factor(ana.table$SEASON)
ana.table$nest.seas<-as.factor(ana.table$nest.seas)
#sex and growth, with CHICK recoded
ana.table = within(ana.table, { 
		CH_B = as.factor(ifelse(CHICK=="B", "1", "0"))} )
ana.table = within(ana.table, {		
		CH_S = as.factor(ifelse(CHICK=="S", "1", "0"))} )
ana.table$nest.seas=as.character(ana.table$nest.seas)
#str(ana.table)


t.test(weight.slope40~did.cr, data=ana.table, var.equal = TRUE) 
t.test(flipper.slope40~did.cr, data=ana.table, var.equal = TRUE) 
t.test(tibiotar.slope35~did.cr, data=ana.table, var.equal = TRUE) 
t.test(foot.slope35~did.cr, data=ana.table, var.equal = TRUE) 
t.test(bill.slopeall~did.cr, data=ana.table, var.equal = TRUE) 

data=ana.table[complete.cases(ana.table[,17]),]

data=subset(ana.table, did.cr=="1")
data=data[complete.cases(data[,17]),]
data=data[complete.cases(data[,21]),]

tapply(data$cr.mass, data$CH_S, mean)

cor.test(data$weight.slope40, data$resid.hatch)
cor.test(data$av.trip.length, data$cr.age)
boxplot(SEASON~tibiotar.slope35, data=data)


ggplot(data, aes(x = sex, y = cr.mass)) + 
					geom_boxplot() +
					geom_point(aes(colour = av.food))+
					facet_grid(~SEASON)

tapply(data$cr.mass, data$SEASON, mean)
tapply(data$cr.mass, data$SEASON, sd)					

ggplot(data, aes(x = SEASON, y = av.food)) + 
					geom_boxplot() +
					geom_point()


#################################################################
scatterplot3d(data$cr.age, data$cr.mass, data$weight.slope40)

s3d <-scatterplot3d(data$cr.age, data$weight.slope40,data$cr.mass, pch=16, highlight.3d=TRUE,
		xlab="Creche age", ylab="Growth rate (mass)", zlab="Creche mass", type="h")
  
fit 	<- lm(cr.mass~	weight.slope40		 +resid.hatch+	sex+SEASON		 , data=data)

s3d$plane3d(fit)

xlab, ylab, zlab
#################################################################

###creching age step 1
sex	<- lmer(cr.age~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
year	<- lmer(cr.age~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
sex_year	<- lmer(cr.age~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
sex.year	<- lmer(cr.age~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)

m.models <- list( )
m.models[[1]] <- 	sex
m.models[[2]] <- 	year
m.models[[3]] <-	sex_year
m.models[[4]] <- 	sex.year

Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)


#################################################################
##creching age step 2

mass_sex_year	<- lmer(cr.age~	weight.slope40		 +	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
flip_sex_year	<- lmer(cr.age~	flipper.slope40		 +	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)
tib_sex_year	<- lmer(cr.age~	tibiotar.slope35		 +	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)
hatch_sex_year	<- lmer(cr.age~	resid.hatch+	sex+SEASON		 +	(1 | nest.seas), data=data, REML = TRUE)
sex_year	<- lmer(cr.age~	sex+SEASON	 +	(1 | nest.seas) , data=data, REML = FALSE)
mass_hatch_sex_year	<- lmer(cr.age~	weight.slope40		 +resid.hatch+	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
flip_hatch_sex_year	<- lmer(cr.age~	flipper.slope40		 +resid.hatch+	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)
tib_hatch_sex_year	<- lmer(cr.age~	tibiotar.slope35		 +resid.hatch+	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)

int		<- lmer (cr.age~ 1 +	(1 | nest.seas)  , data=data, REML = FALSE)

m.models <- list( )	
m.models[[1]] <- 	mass_sex_year
m.models[[2]] <- 	flip_sex_year
m.models[[3]] <-	tib_sex_year
m.models[[4]] <- 	sex_year
m.models[[5]] <- 	hatch_sex_year
m.models[[6]] <- 	mass_hatch_sex_year
m.models[[7]] <- 	flip_hatch_sex_year
m.models[[8]] <-	tib_hatch_sex_year
m.models[[9]] <-	int

Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

summary(mass_sex_year)
confint(mass_sex_year, method="boot")

summary(hatch_sex_year)
confint(hatch_sex_year, method="boot")

#--
summary(mass_hatch_sex_year)
confint(mass_hatch_sex_year, method="boot")

summary(flip_hatch_sex_year)
confint(flip_hatch_sex_year, method="boot")

summary(tib_hatch_sex_year)
confint(tib_hatch_sex_year, method="boot")
##-- VCA
hatch_sex_year.reml <- lmer(cr.age~	resid.hatch+	sex+SEASON		 +	(1 | nest.seas), data=data, REML = TRUE)	
int.reml <- lmer(cr.age~ 1 + (1 | nest.seas) , data=data, REML = TRUE)		
VarCorr(hatch_sex_year.reml)
VarCorr(int.reml)


((1.0837-0.86213)/1.0837)*100



#################################################################

fm1 <- lmer(cr.age~	resid.hatch+	sex+SEASON		 +	(1 | nest.seas), data=data, REML = TRUE)

newdat <- expand.grid(
		SEASON=c("1213","1314")
		,sex=c("F","M")
		,resid.hatch=seq(-10,10)
		,cr.age= 0
)

mm <- model.matrix(terms(fm1),newdat)
newdat$pred <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
#tvar1 <- pvar1+VarCorr(fm1)$Subject[1]  ## must be adapted for more complex models
newdat <- data.frame(
    newdat
    , plo = newdat$pred-2*sqrt(pvar1)
    , phi = newdat$pred+2*sqrt(pvar1)
)					
cr.age.pred=newdat			

#################################################################

ggplot(data=subset(cr.age.pred, SEASON=="1213"&sex=="F"), aes(x=resid.hatch, y=pred))+
	geom_point(size=1)+
	geom_line(size=1) +
	ylab("Creching age")+
	xlab("Hatch date")+
	geom_line(aes(y=plo), linetype="dotted",size=1) +
	geom_line(aes(y=phi), linetype="dotted",size=1) +
	#geom_ribbon( aes(ymin = lci, ymax = uci), alpha = .15)+
	guides(shape=FALSE)+
	geom_vline(xintercept = 0, linetype = "longdash", alpha=0.3, size=1)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					#,plot.margin = unit( c(-1,1,0,0) , units = "lines" )
					,axis.ticks.y = element_blank()
					,axis.text.x = element_blank()
					,axis.text.y = element_text(size=15)
					,axis.title = element_text(size=20)
					)	

#################################################################
##creching mass step 1
mass.sex	<- lmer(cr.mass~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
mass.year	<- lmer(cr.mass~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
mass.sex_year	<- lmer(cr.mass~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
mass.sex.year	<- lmer(cr.mass~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
m.models <- list( )
m.models[[1]] <- 	mass.sex
m.models[[2]] <- 	mass.year
m.models[[3]] <-	mass.sex_year
m.models[[4]] <- 	mass.sex.year
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

##--##creching mass step 2
mass.growth_sex_year	<- lm(cr.mass~	weight.slope40		 +	sex+SEASON		 , data=data)
mass.sex_year	<- lm(cr.mass~	sex+SEASON		 , data=data)
mass.hatch_sex_year	<- lm(cr.mass~	resid.hatch+	sex+SEASON		 , data=data)
mass.growth_hatch_sex_year	<- lm(cr.mass~	weight.slope40		 +resid.hatch+	sex+SEASON		 , data=data)
mass.age_hatch_sex_year	<- lm(cr.mass~cr.age+resid.hatch+sex+SEASON		 , data=data)
mass.age_growth_hatch_sex_year	<- lm(cr.mass~	cr.age + weight.slope40		 +resid.hatch+	sex+SEASON		 , data=data)
mass.age_sex_year	<- lm(cr.mass~	cr.age +	sex+SEASON		 , data=data)
mass.age_growth_sex_year	<- lm(cr.mass~	cr.age + weight.slope40		 +	sex+SEASON		 , data=data)

m.models <- list( )	
m.models[[1]] <- 	mass.growth_sex_year
m.models[[2]] <- 	mass.sex_year
m.models[[3]] <-	mass.hatch_sex_year	
m.models[[4]] <-	mass.growth_hatch_sex_year	
m.models[[5]] <-	mass.age_hatch_sex_year	
m.models[[6]] <-	mass.age_growth_hatch_sex_year	
m.models[[7]] <-	mass.age_sex_year	
m.models[[8]] <-	mass.age_growth_sex_year	

Modnames <- paste("mod", 1:length(m.models), sep = " ")

Modnames <- c("MassGrowth + Sex + Year", "Sex + Year", "Hatch + Sex + Year", "MassGrowth + Hatch + Sex + Year", "Age + Hatch + Sex + Year", "Age + MassGrowth + Hatch + Sex + Year", "Age + Sex + Year", "Age + MassGrowth + Sex + Year")


print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

summary(mass.growth_hatch_sex_year)
confint(mass.growth_hatch_sex_year)

##--
#################################################################
## no support for RE
##		+	(1 | nest.seas)  , data=data, REML = FALSE) 
flip.sex	<- lm(cr.flip~	sex	, data=data)
flip.year	<- lm(cr.flip~	SEASON	, data=data)
flip.sex_year	<- lm(cr.flip~	sex+SEASON	, data=data)
flip.sex.year	<- lm(cr.flip~	sex*SEASON	 , data=data)
m.models <- list( )
m.models[[1]] <- 	flip.sex
m.models[[2]] <- 	flip.year
m.models[[3]] <-	flip.sex_year
m.models[[4]] <- 	flip.sex.year
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)


##--
flip.growth_sex_year	<- lm(cr.flip~	flipper.slope40		 +	sex+SEASON		 , data=data)
flip.sex_year	<- lm(cr.flip~	sex+SEASON		  , data=data)
flip.hatch_sex_year	<- lm(cr.flip~	resid.hatch+	sex+SEASON			 , data=data)
flip.growth_hatch_sex_year	<- lm(cr.flip~	flipper.slope40		 +resid.hatch+	sex+SEASON			 , data=data)
flip.age_hatch_sex_year	<- lm(cr.flip~	cr.age + resid.hatch+	sex+SEASON		 , data=data)
flip.age_growth_hatch_sex_year	<- lm(cr.flip~	cr.age + flipper.slope40		 +resid.hatch+	sex+SEASON		 , data=data)
flip.age_sex_year	<- lm(cr.flip~	cr.age +	sex+SEASON		 , data=data)
flip.age_growth_sex_year	<- lm(cr.flip~	cr.age + flipper.slope40		 +	sex+SEASON		 , data=data)

m.models <- list( )	
m.models[[1]] <- 	flip.growth_sex_year
m.models[[2]] <- 	flip.sex_year
m.models[[3]] <-	flip.hatch_sex_year
m.models[[4]] <-	flip.growth_hatch_sex_year
m.models[[5]] <-	flip.age_hatch_sex_year	
m.models[[6]] <-	flip.age_growth_hatch_sex_year	
m.models[[7]] <-	flip.age_sex_year	
m.models[[8]] <-	flip.age_growth_sex_year	

#Modnames <- paste("mod", 1:length(m.models), sep = " ")

Modnames <- c("FlipGrowth + Sex + Year", "Sex + Year", "Hatch + Sex + Year", "FlipGrowth + Hatch + Sex + Year", "Age + Hatch + Sex + Year", "Age + FlipGrowth + Hatch + Sex + Year", "Age + Sex + Year", "Age + FlipGrowth + Sex + Year")

print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

summary(flip.growth_hatch_sex_year)
confint(flip.growth_hatch_sex_year)

#################################################################

tib.sex	<- lmer(cr.tib~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
tib.year	<- lmer(cr.tib~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
tib.sex_year	<- lmer(cr.tib~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
tib.sex.year	<- lmer(cr.tib~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
m.models <- list( )
m.models[[1]] <- 	tib.sex
m.models[[2]] <- 	tib.year
m.models[[3]] <-	tib.sex_year
m.models[[4]] <- 	tib.sex.year
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)
##--
tib.growth_sex_year	<- lm(cr.tib~	tibiotar.slope35		 +	sex+SEASON, data=data)
tib.sex_year	<- lm(cr.tib~	sex+SEASON		 , data=data)
tib.hatch_sex_year	<- lm(cr.tib~	resid.hatch+	sex+SEASON  , data=data)
tib.growth_hatch_sex_year	<- lm(cr.tib~	tibiotar.slope35		 +resid.hatch+	sex+SEASON  , data=data)
tib.age_hatch_sex_year	<- lm(cr.tib~cr.age+resid.hatch+sex+SEASON		 , data=data)
tib.age_growth_hatch_sex_year	<- lm(cr.tib~cr.age+tibiotar.slope35+resid.hatch+sex+SEASON		 , data=data)
tib.age_sex_year	<- lm(cr.tib~cr.age+sex+SEASON		 , data=data)
tib.age_growth_sex_year	<- lm(cr.tib~cr.age+tibiotar.slope35+sex+SEASON		 , data=data)

m.models <- list( )	
m.models[[1]] <- 	tib.growth_sex_year
m.models[[2]] <- 	tib.sex_year
m.models[[3]] <-	tib.hatch_sex_year
m.models[[4]] <-	tib.growth_hatch_sex_year
m.models[[5]] <-	tib.age_hatch_sex_year	
m.models[[6]] <-	tib.age_growth_hatch_sex_year	
m.models[[7]] <-	tib.age_sex_year	
m.models[[8]] <-	tib.age_growth_sex_year	

#Modnames <- paste("mod", 1:length(m.models), sep = " ")
Modnames <- c("TibGrowth + Sex + Year", "Sex + Year", "Hatch + Sex + Year", "TibGrowth + Hatch + Sex + Year", "Age + Hatch + Sex + Year", "Age + TibGrowth + Hatch + Sex + Year", "Age + Sex + Year", "Age + TibGrowth + Sex + Year")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

summary(tib.hatch_sex_year)
confint(tib.hatch_sex_year)

summary(tib.growth_hatch_sex_year)
confint(tib.growth_hatch_sex_year)

#################################################################

foot.sex	<- lmer(cr.foot~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
foot.year	<- lmer(cr.foot~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
foot.sex_year	<- lmer(cr.foot~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
foot.sex.year	<- lmer(cr.foot~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
m.models <- list( )
m.models[[1]] <- 	foot.sex
m.models[[2]] <- 	foot.year
m.models[[3]] <-	foot.sex_year
m.models[[4]] <- 	foot.sex.year
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

##--
foot.growth_sex_year	<- lmer(cr.foot~	foot.slope35		 +	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)
foot.sex_year	<- lmer(cr.foot~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
foot.int	<- lmer(cr.foot~	1 +	(1 | nest.seas) , data=data, REML = FALSE)
m.models <- list( )	
m.models[[1]] <- 	foot.growth_sex_year
m.models[[2]] <- 	foot.sex_year
m.models[[3]] <-	foot.int
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)


#################################################################

bill.sex	<- lmer(cr.bill~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
bill.year	<- lmer(cr.bill~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
bill.sex_year	<- lmer(cr.bill~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
bill.sex.year	<- lmer(cr.bill~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
m.models <- list( )
m.models[[1]] <- 	bill.sex
m.models[[2]] <- 	bill.year
m.models[[3]] <-	bill.sex_year
m.models[[4]] <- 	bill.sex.year
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)

##--
bill.growth_sex_year	<- lmer(cr.bill~	bill.slopeall		 +	sex+SEASON		 +	(1 | nest.seas), data=data, REML = FALSE)
bill.sex_year	<- lmer(cr.bill~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
bill.int	<- lmer(cr.bill~	1 +	(1 | nest.seas) , data=data, REML = FALSE)
m.models <- list( )	
m.models[[1]] <- 	bill.growth_sex_year
m.models[[2]] <- 	bill.sex_year
m.models[[3]] <-	bill.int
Modnames <- paste("mod", 1:length(m.models), sep = " ")
print(aictab(cand.set = m.models, modnames = Modnames, sort = TRUE), digits = 4, LL = TRUE)
##--





ggplot(data, aes(x = flipper.slope40, y = cr.age)) +
					geom_point()+
					stat_smooth(method = "lm")


					
#################################################################

mod.a <- lm(cr.mass~	cr.age + weight.slope40		 +	sex+SEASON		 , data=data)
newdat.a <- expand.grid(
SEASON="1213"
		,sex="F"
		,cr.age= 20
		,cr.mass=0
		,weight.slope40=data$weight.slope40
)
pred.a=predict.lm(mod.a, newdat.a, interval="confidence") 
cr.mass.pred.mass.gr=cbind(newdat.a, pred.a)
			
					
					
							
mod.b <- lm(cr.mass~	cr.age + weight.slope40		 +	sex+SEASON		 , data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= data$cr.age
		,cr.mass=0
		,weight.slope40=64
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.mass.pred.cr.age=cbind(newdat, pred)

##		
mod <- lm(cr.flip~	cr.age + flipper.slope40 +	sex+SEASON, data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= 20
		,cr.flip=0
		,flipper.slope40=data$flipper.slope40
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.flip.pred.flip.gr=cbind(newdat, pred)



mod <- lm(cr.flip~	cr.age + flipper.slope40 +	sex+SEASON, data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= data$cr.age
		,cr.flip=0
		,flipper.slope40=3.8
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.flip.pred.cr.age=cbind(newdat, pred)


library(gridExtra)
			
a=ggplot(data=cr.mass.pred.mass.gr, aes(x=weight.slope40, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("MASS")+
	ylab("g")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	ylim(500,2200)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					,plot.margin = unit( c(0,1,0,-0.5) , units = "lines" )
					,axis.ticks = element_blank()
					,axis.text.x = element_blank()
					,axis.text.y = element_text(size=15)
					,axis.title.y = element_text(size=20)
					)

b=ggplot(data=cr.flip.pred.flip.gr, aes(x=flipper.slope40, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("FLIPPER")+
	ylab("mm")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	ylim(90,145)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					,plot.margin = unit( c(-1,1,0,0) , units = "lines" )
					,axis.ticks = element_blank()
					,axis.text.x = element_blank()
					,axis.text.y = element_text(size=15)
					,axis.title.y = element_text(size=20)
					)



c=ggplot(data=cr.mass.pred.cr.age, aes(x=cr.age, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("MASS")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	ylim(500,2200)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					,plot.margin = unit( c(0,1,0,0) , units = "lines" )
					,axis.ticks = element_blank()
					,axis.text.x = element_blank()
					,axis.text.y = element_blank()
					,axis.title.y = element_blank()
					)

d=ggplot(data=cr.flip.pred.cr.age, aes(x=cr.age, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("FLIPPER")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	ylim(90,145)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,panel.grid.major = element_blank()
					,strip.text = element_blank() 
					,strip.background = element_blank()
					,plot.margin = unit( c(-1,1,0,-0.5) , units = "lines" )
					,axis.ticks = element_blank()
					,axis.text.x = element_blank()
					,axis.text.y = element_blank()
					)
					
grid.arrange(arrangeGrob(a, b, ncol=1, sub=textGrob("Growth rate", vjust=0, gp=gpar(fontsize=15))), 
			 arrangeGrob(c, d, ncol=1, sub=textGrob("Creche age", vjust=0, gp=gpar(fontsize=15))),
			left=textGrob(expression("Size at creching"), rot=90, gp=gpar(fontsize=20)),
			ncol=2)

, widths=c(4.5, 4)
##		
mod <- lm(cr.tib~	cr.age + tibiotar.slope35		 +	sex+SEASON		 , data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= 20
		,cr.tib=0
		,tibiotar.slope35=data$tibiotar.slope35
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.tib.pred=cbind(newdat, pred)



library(gridExtra)
			
a=ggplot(data=cr.mass.pred, aes(x=weight.slope40, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("MASS")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

b=ggplot(data=cr.flip.pred, aes(x=flipper.slope40, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("FLIPPER")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

					
c=ggplot(data=cr.tib.pred, aes(x=tibiotar.slope35, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("TIBIOTARSUS")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

					
grid.arrange(a, b, c,
			left=textGrob(expression("Size at creching"), rot=90, gp=gpar(fontsize=15)),
			sub=textGrob("Growth rate", vjust=0, gp=gpar(fontsize=15)),
			ncol=1)

#################################################

					
mod <- lm(cr.mass~	cr.age + weight.slope40		 +	sex+SEASON		 , data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= data$cr.age
		,cr.mass=0
		,weight.slope40=64
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.mass.pred=cbind(newdat, pred)
##		
mod <- lm(cr.flip~	cr.age + flipper.slope40 +	sex+SEASON, data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= data$cr.age
		,cr.flip=0
		,flipper.slope40=3.8
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.flip.pred=cbind(newdat, pred)
##		
mod <- lm(cr.tib~	cr.age + tibiotar.slope35		 +	sex+SEASON		 , data=data)
newdat <- expand.grid(
		SEASON="1213"
		,sex="F"
		,cr.age= data$cr.age
		,cr.tib=0
		,tibiotar.slope35=2.6
)
pred=predict.lm(mod, newdat, interval="confidence") 
cr.tib.pred=cbind(newdat, pred)



library(gridExtra)
			
a=ggplot(data=cr.mass.pred, aes(x=cr.age, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("MASS")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

b=ggplot(data=cr.flip.pred, aes(x=cr.age, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("FLIPPER")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

					
c=ggplot(data=cr.tib.pred, aes(x=cr.age, y=fit))+
	geom_point(size=0)+
	geom_line(size=1) +
	ggtitle("TIBIOTARSUS")+
	ylab("")+
	xlab("")+
	geom_line(aes(y=lwr), linetype="dotted",size=1) +
	geom_line(aes(y=upr), linetype="dotted",size=1) +
	guides(shape=FALSE)+
	theme_bw()+
				theme(
					plot.background = element_blank()
					,panel.grid.minor = element_blank()
					,strip.background = element_blank()
					,plot.margin = unit( c(1,1,0,0) , units = "lines" )
					)

					
grid.arrange(a, b, c,
			left=textGrob(expression("Size at creching"), rot=90, gp=gpar(fontsize=15)),
			sub=textGrob("Creching age", vjust=0, gp=gpar(fontsize=15)),
			ncol=1)

