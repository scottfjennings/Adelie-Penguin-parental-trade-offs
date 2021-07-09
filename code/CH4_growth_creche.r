

library(tidyverse)
library(lme4)
library(AICcmodavg)
library(grid)
options(scipen = 999)

# read data ----
#  ana.table.full.csv is created by the sequence of numbered code files in THESIS\thesis_data_work\code_files\analysis_data_prep
#ana.table=read.csv("data/ana.table.full.csv")
ana.table <- readRDS("data/ana_table_full")
#ana.table.old=read.csv("data/ana.table.full.old.csv")%>%  rename(CHICK_ID = id)
# data management ----
ana.table <- ana.table %>% 
  mutate(SEASON = as.factor(SEASON),
         CH_B = as.factor(ifelse(CHICK=="B", "1", "0")),
         CH_S = as.factor(ifelse(CHICK=="S", "1", "0")))


data = ana.table %>% 
  filter(!is.na(cr.seasday), weight.slope40 >= 0, !is.na(sex))



#data.old <- ana.table.old %>% filter(did.cr == "1", !is.na(sex)) 



# analysis ----
# first creching age as response variable ----
# creching age model selection step 1 ----


sex	<- lmer(cr.age~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
year	<- lmer(cr.age~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
sex_year	<- lmer(cr.age~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
sex.year	<- lmer(cr.age~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)

(cr_age_aic1 <- aictab(list(sex, year, sex_year, sex.year), c("sex", "year", "sex_year", "sex.year")))

saveRDS(cr_age_aic1, "rds/cr_age_aic1")

# year model is best supported, but all are competitive
# sex.year is "singular"

# creching age model selection step 2 ----

mass_year	<- lmer(cr.age ~	weight.slope40 + SEASON +	(1 | nest.seas) , data=data, REML = FALSE)
flip_year	<- lmer(cr.age ~	flipper.slope40	+ SEASON + (1 | nest.seas), data=data, REML = FALSE)
tib_year	<- lmer(cr.age ~	tibiotar.slope35 +	SEASON +	(1 | nest.seas), data=data, REML = FALSE)
hatch_year	<- lmer(cr.age ~	resid.hatch + SEASON +	(1 | nest.seas), data=data, REML = FALSE)
year	<- lmer(cr.age ~	SEASON +	(1 | nest.seas) , data=data, REML = FALSE)
mass_hatch_year	<- lmer(cr.age ~	weight.slope40 + resid.hatch +	SEASON + (1 | nest.seas) , data=data, REML = FALSE)
flip_hatch_year	<- lmer(cr.age ~	flipper.slope40 + resid.hatch +	SEASON + (1 | nest.seas), data=data, REML = FALSE)
tib_hatch_year	<- lmer(cr.age ~	tibiotar.slope35 + resid.hatch +	SEASON + (1 | nest.seas), data=data, REML = FALSE)
int		<- lmer (cr.age~ 1 +	(1 | nest.seas)  , data=data, REML = FALSE)

# all models with resid.hatch have singular fit

(cr_age_aic <- aictab(list(mass_year, flip_year, tib_year, hatch_year, year, mass_hatch_year, flip_hatch_year, tib_hatch_year, int), 
       c("mass_year", "flip_year", "tib_year", "hatch_year", "year", "mass_hatch_year", "flip_hatch_year", "tib_hatch_year", "int")))

saveRDS(cr_age_aic, "rds/cr_age_aic")

# same candidate set but without random effect
mass_year.lm	<- lm(cr.age ~	weight.slope40 + SEASON, data=data)
flip_year.lm	<- lm(cr.age ~	flipper.slope40	+ SEASON, data=data)
tib_year.lm	<- lm(cr.age ~	tibiotar.slope35 +	SEASON, data=data)
hatch_year.lm	<- lm(cr.age ~	resid.hatch + SEASON, data=data)
year.lm	<- lm(cr.age ~	SEASON, data=data)
mass_hatch_year.lm	<- lm(cr.age ~	weight.slope40 + resid.hatch +	SEASON, data=data)
flip_hatch_year.lm	<- lm(cr.age ~	flipper.slope40 + resid.hatch +	SEASON, data=data)
tib_hatch_year.lm	<- lm(cr.age ~	tibiotar.slope35 + resid.hatch +	SEASON, data=data)
int.lm		<- lm(cr.age~ 1, data=data)

(cr_age_aic.lm <- aictab(list(mass_year.lm, flip_year.lm, tib_year.lm, hatch_year.lm, year.lm, mass_hatch_year.lm, flip_hatch_year.lm, tib_hatch_year.lm, int.lm), 
       c("mass_year.lm", "flip_year.lm", "tib_year.lm", "hatch_year.lm", "year.lm", "mass_hatch_year.lm", "flip_hatch_year.lm", "tib_hatch_year.lm", "int.lm")))

saveRDS(cr_age_aic.lm, "rds/cr_age_aic_lm")

# get same model selection results
# and similar beta estimates

# creching age output from best models ----
# coefficients
get_cr_mass_coefCI <- function() {
mass_hatch_year_coefCI <- cbind(fixef(mass_hatch_year), confint(mass_hatch_year, method="boot")[3:6,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "mass_hatch_year")


hatch_year_coefCI <- cbind(fixef(hatch_year), confint(hatch_year, method="boot")[3:5,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "hatch_year")

tib_hatch_year_coefCI <- cbind(fixef(tib_hatch_year), confint(tib_hatch_year, method="boot")[3:6,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "tib_hatch_year")

flip_hatch_year_coefCI <- cbind(fixef(flip_hatch_year), confint(flip_hatch_year, method="boot")[3:6,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "flip_hatch_year")
# combine
coefCI <- rbind(mass_hatch_year_coefCI,
                hatch_year_coefCI,
                tib_hatch_year_coefCI,
                flip_hatch_year_coefCI) %>% 
  mutate(varb = ifelse(varb %in% c("weight.slope40", "tibiotar.slope35", "flipper.slope40"), "growth", varb))
}

# lm
get_cr_mass_coefCI_lm <- function() {
hatch_year.lm_coefCI <- cbind(coef(hatch_year.lm), confint(hatch_year.lm)) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "hatch_year.lm")
  
mass_hatch_year.lm_coefCI <- cbind(coef(mass_hatch_year.lm), confint(mass_hatch_year.lm)) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "mass_hatch_year.lm")


tib_hatch_year.lm_coefCI <- cbind(coef(tib_hatch_year.lm), confint(tib_hatch_year.lm)) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "tib_hatch_year.lm")

flip_hatch_year.lm_coefCI <- cbind(coef(flip_hatch_year.lm), confint(flip_hatch_year.lm)) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "flip_hatch_year.lm")
# combine 
coefCI <- rbind(mass_hatch_year.lm_coefCI,
                hatch_year.lm_coefCI,
                tib_hatch_year.lm_coefCI,
                flip_hatch_year.lm_coefCI) %>% 
  mutate(varb = ifelse(varb %in% c("weight.slope40", "tibiotar.slope35", "flipper.slope40"), "growth", varb))
}

coefCI.lm <- get_cr_mass_coefCI_lm()
# paste coefs and CI and pivot wider
coefCI_wider <- coefCI %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 

coefCI.lm_wider <- coefCI.lm %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 



# creching age lmer VCA ----
# residual variance of intercept only model
int_var <- VarCorr(int) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "int")

# residual variance of best models

mass_hatch_year_var <- VarCorr(mass_hatch_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "mass_hatch_year")
#((5.3157  -4.5217  )/5.3157  )*100

hatch_year_var <- VarCorr(hatch_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>%
  mutate(mod = "hatch_year")

#((4.5997 - 4.5217)/4.5997)*100

tib_hatch_year_var <- VarCorr(tib_hatch_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "tib_hatch_year")
#((4.5787 - 4.5217)/4.5787)*100

flip_hatch_year_var <- VarCorr(flip_hatch_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "flip_hatch_year")
#((4.5831 - 4.5217)/4.5831)*100

# combine variances from all competitive models, calculate variance explained
vars = rbind(mass_hatch_year_var, hatch_year_var, tib_hatch_year_var, flip_hatch_year_var) %>% 
  mutate(int.var = int_var$sdcor,
         var.expl = ((int.var - sdcor)/int.var) * 100)

# creching age rsquared for lms ----
lm.adj.r2 <- rbind(data.frame(mod = "hatch_year.lm", adj.r2 = summary(hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "mass_hatch_year.lm", adj.r2 = summary(mass_hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "tib_hatch_year.lm", adj.r2 = summary(tib_hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "flip_hatch_year.lm", adj.r2 = summary(flip_hatch_year.lm)$adj.r.squared))

# creching age combine coefficients, aic variables, variance explained ----
# mixed mod
coef_aic <- coefCI_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(cr_age_aic) %>% 
  left_join(select(vars, Modnames = mod, var.expl)) %>% 
  select(Modnames, growth, resid.hatch, SEASON1314, K, Delta_AICc, AICcWt, LL, var.expl) %>% 
  mutate(var.expl = round(var.expl, 1),
         across(c(Delta_AICc, AICcWt), ~round(., 2))) %>% 
  arrange(Delta_AICc)

saveRDS(coef_aic, "rds/cr_age_best")

# lm
coef_aic.lm <- coefCI.lm_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(cr_age_aic.lm) %>% 
  left_join(select(lm.adj.r2, Modnames = mod, adj.r2)) %>% 
  select(Modnames, growth, resid.hatch, SEASON1314, K, Delta_AICc, AICcWt, LL, adj.r2) %>% 
  mutate(adj.r2 = round(adj.r2, 3),
         across(c(Delta_AICc, AICcWt), ~round(., 2))) %>% 
  arrange(Delta_AICc)

saveRDS(coef_aic.lm, "rds/cr_age_best_lm")
# creching age plot ----

fm1 <- lmer(cr.age ~	weight.slope40 + resid.hatch +	SEASON + (1 | nest.seas) , data=data, REML = TRUE)

newdat <- expand.grid(
		SEASON=c("1213","1314")
		,sex=c("F")
		,weight.slope40 = seq(0, max(data$weight.slope40), length.out = 10)
		,resid.hatch=seq(-6,10)
		,cr.age= mean(data$cr.age)
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


ggplot(cr.age.pred, group = factor(weight.slope40))+
	geom_line(aes(x = resid.hatch, y = pred, color = factor(weight.slope40)), size=1) +
	ylab("Creching age")+
	xlab("Hatch date")+
	#geom_line(aes(y=plo, color = weight.slope40), linetype="dotted",size=1) +
	#geom_line(aes(y=phi, color = weight.slope40), linetype="dotted",size=1) +
	#geom_ribbon( aes(ymin = lci, ymax = uci), alpha = .15)+
	guides(shape=FALSE)+
	geom_vline(xintercept = 0, linetype = "longdash", alpha=0.3, size=1)+
	theme_bw() +
  facet_wrap(~SEASON)




  
# creching mass model selection step 1 ----


mass.sex	<- lmer(cr.mass~	sex		 +	(1 | nest.seas)  , data=data, REML = FALSE)
mass.year	<- lmer(cr.mass~	SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
mass.sex_year	<- lmer(cr.mass~	sex+SEASON		 +	(1 | nest.seas) , data=data, REML = FALSE)
mass.sex.year	<- lmer(cr.mass~	sex*SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)

(cr_mass_aic1 <- aictab(list(mass.sex, mass.year, mass.sex_year, mass.sex.year), c("mass.sex", "mass.year", "mass.sex_year", "mass.sex.year")))
# sex_year best
saveRDS(cr_mass_aic1, "rds/cr_mass_aic1")


# mass.sex_year best supported

# creching flipper length model selection step 1 ----

flip.sex	<- lmer(cr.flip~	sex	+	(1 | nest.seas)  , data=data, REML = FALSE)
flip.year	<- lmer(cr.flip~	SEASON	+	(1 | nest.seas)  , data=data, REML = FALSE)
flip.sex_year	<- lmer(cr.flip~	sex+SEASON	+	(1 | nest.seas)  , data=data, REML = FALSE)
flip.sex.year	<- lmer(cr.flip~	sex*SEASON	 +	(1 | nest.seas)  , data=data, REML = FALSE)

(cr_flip_aic1 <- aictab(list(flip.sex, flip.year, flip.sex_year, flip.sex.year), c("flip.sex", "flip.year", "flip.sex_year", "flip.sex.year")))

# year best
saveRDS(cr_flip_aic1, "rds/cr_flip_aic1")

#creching tib length model selection step 1 ----

tib.sex	<- lmer(cr.tib~	sex	+	(1 | nest.seas)  , data=data, REML = FALSE)
tib.year	<- lmer(cr.tib~	SEASON	+	(1 | nest.seas)  , data=data, REML = FALSE)
tib.sex_year	<- lmer(cr.tib~	sex+SEASON	+	(1 | nest.seas)  , data=data, REML = FALSE)
tib.sex.year	<- lmer(cr.tib~	sex*SEASON	 +	(1 | nest.seas)  , data=data, REML = FALSE)

(cr_tib_aic1 <- aictab(list(tib.sex, tib.year, tib.sex_year, tib.sex.year), c("tib.sex", "tib.year", "tib.sex_year", "tib.sex.year")))
# sex_year best
saveRDS(cr_tib_aic1, "rds/cr_tib_aic1")


# creche size model selection step 2  ----
# same best structure for all morph in step 1, so can do step 2 with same function for all morphs
cr_size_step2_cand_set_sex_year <- function(data) {

growth_sex_year	<- lmer(cr.size~	growth.rate		 +	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
sex_year	<- lmer(cr.size~	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
hatch_sex_year	<- lmer(cr.size~	resid.hatch+	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
growth_hatch_sex_year	<- lmer(cr.size~	growth.rate		 +resid.hatch+	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
age_sex_year	<- lmer(cr.size~	cr.age +	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
age_growth_sex_year	<- lmer(cr.size~	cr.age + growth.rate		 +	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
int	<- lmer(cr.size~	1		 +	(1 | nest.seas)  , data=data, REML = FALSE)

cr_size_aic <- aictab(list(growth_sex_year, sex_year, hatch_sex_year, growth_hatch_sex_year, age_sex_year, age_growth_sex_year, int), c("growth_sex_year", "sex_year", "hatch_sex_year", "growth_hatch_sex_year", "age_sex_year", "age_growth_sex_year", "int")) %>% 
  data.frame() 

return(cr_size_aic)
}

cr_mass_aic <- data %>%  
      select(cr.size = cr.mass, growth.rate = weight.slope40, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_sex_year() 
cr_mass_aic
# age_growth_sex_year best by 20.41 aicc

saveRDS(cr_mass_aic, "rds/cr_mass_aic")

cr_tib_aic <- data %>% 
      select(cr.size = cr.tib, growth.rate = tibiotar.slope35, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_sex_year() 
# age_sex_year best
# age_growth_sex_year daicc = 0.62 but uninformative. next best daicc = 13.18

saveRDS(cr_tib_aic, "rds/cr_tib_aic")


cr_size_step2_cand_set_year <- function(data) {

growth_year	<- lmer(cr.size ~	growth.rate + SEASON + (1 | nest.seas), data=data, REML = FALSE)
year	<- lmer(cr.size~	SEASON + (1 | nest.seas), data=data, REML = FALSE)
hatch_year	<- lmer(cr.size ~ resid.hatch +	SEASON + (1 | nest.seas), data=data, REML = FALSE)
growth_hatch_year	<- lmer(cr.size ~	growth.rate + resid.hatch +	SEASON + (1 | nest.seas), data=data, REML = FALSE)
age_year	<- lmer(cr.size ~ cr.age + SEASON +	(1 | nest.seas), data=data, REML = FALSE)
age_growth_year	<- lmer(cr.size ~	cr.age + growth.rate + SEASON + (1 | nest.seas), data=data, REML = FALSE)
int	<- lmer(cr.size ~ 1 + (1 | nest.seas), data=data, REML = FALSE)

cr_size_aic <- aictab(list(growth_year, year, hatch_year, growth_hatch_year, age_year, age_growth_year, int), c("growth_year", "year", "hatch_year", "growth_hatch_year", "age_year", "age_growth_year", "int")) %>% 
  data.frame()

return(cr_size_aic)
}

cr_flip_aic <- data %>% 
      select(cr.size = cr.flip, growth.rate = flipper.slope40, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_year()
# age_growth_year best by 11.15 aicc

saveRDS(cr_flip_aic, "rds/cr_flip_aic")

# best models output ----
# creching mass size best model output ----

mass.age_growth_sex_year	<- lmer(cr.mass~	cr.age + weight.slope40		 +	sex+SEASON		 +	(1 | nest.seas)  , data=data, REML = FALSE)
mass.int	<- lmer(cr.mass~	1		 +	(1 | nest.seas)  , data=data, REML = FALSE)
# coefficients
mass.age_growth_sex_year_coefCI <- cbind(fixef(mass.age_growth_sex_year), confint(mass.age_growth_sex_year, method="boot")[3:7,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "age_growth_sex_year")

# VCA
# residual variance of intercept only
mass_best_var_expl <- function() {
mass.int_var <- VarCorr(mass.int) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "int")
# residual variance of best model
mass.age_growth_sex_year_var <- VarCorr(mass.age_growth_sex_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "age_growth_sex_year")

# calculate percent variance explained
mass.vars = mass.age_growth_sex_year_var %>% 
  mutate(int.var = mass.int_var$sdcor,
         var.expl = ((int.var - sdcor)/int.var) * 100)
}

mass_vars <- mass_best_var_expl()

mass.coefCI_wider <- mass.age_growth_sex_year_coefCI %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 


mass.coef_aic <- mass.coefCI_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(., cr_mass_aic) %>% 
  left_join(select(mass_vars, Modnames = mod, var.expl)) %>% 
  select(Modnames, cr.age, weight.slope40, sexM, SEASON1314, K, Delta_AICc, AICcWt, LL, var.expl) %>% 
  mutate(var.expl = round(var.expl, 1),
         across(c(Delta_AICc, AICcWt), ~round(., 2)))

saveRDS(mass.coef_aic, "rds/cr_mass_best")


##--
# creching flipper length best model output ----

flip.age_growth_year	<- lmer(cr.flip ~	cr.age + flipper.slope40 + SEASON +	(1 | nest.seas), data=data, REML = FALSE)
flip.int	<- lmer(cr.flip ~	1 +	(1 | nest.seas)  , data=data, REML = FALSE)
# coefficients
flip.age_growth_year_coefCI <- cbind(fixef(flip.age_growth_year), confint(flip.age_growth_year, method="boot")[3:6,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "age_growth_year")

# VCA
# residual variance of intercept only
flip.int_var <- VarCorr(flip.int) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "int")
# residual variance of best model
flip.age_growth_year_var <- VarCorr(flip.age_growth_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "age_growth_year")

# calculate percent variance explained
flip.vars = flip.age_growth_year_var %>% 
  mutate(int.var = flip.int_var$sdcor,
         var.expl = ((int.var - sdcor)/int.var) * 100)


flip.coefCI_wider <- flip.age_growth_year_coefCI %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 


flip.coef_aic <- flip.coefCI_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(., cr_flip_aic) %>% 
  left_join(select(flip.vars, Modnames = mod, var.expl)) %>% 
  select(Modnames, cr.age, flipper.slope40, SEASON1314, K, Delta_AICc, AICcWt, LL, var.expl) %>% 
  mutate(var.expl = round(var.expl, 1),
         across(c(Delta_AICc, AICcWt), ~round(., 2)))


saveRDS(flip.coef_aic, "rds/cr_flip_best")




# creching tib length best model output ----

tib.age_sex_year	<- lmer(cr.tib ~ cr.age + sex + SEASON +	(1 | nest.seas)  , data=data, REML = FALSE)
tib.int	<- lmer(cr.tib ~	1		 +	(1 | nest.seas)  , data=data, REML = FALSE)
# coefficients
tib.age_sex_year_coefCI <- cbind(fixef(tib.age_sex_year), confint(tib.age_sex_year, method="boot")[3:6,]) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "age_sex_year")

# VCA
# residual variance of intercept only
tib.int_var <- VarCorr(tib.int) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "int")
# residual variance of best model
tib.age_sex_year_var <- VarCorr(tib.age_sex_year) %>% 
  data.frame() %>% 
  filter(grp == "Residual") %>% 
  select(sdcor) %>% 
  mutate(mod = "age_sex_year")

# calculate percent variance explained
tib.vars = tib.age_sex_year_var %>% 
  mutate(int.var = tib.int_var$sdcor,
         var.expl = ((int.var - sdcor)/int.var) * 100)


tib.coefCI_wider <- tib.age_sex_year_coefCI %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 


tib.coef_aic <- tib.coefCI_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(., cr_tib_aic) %>% 
  left_join(select(tib.vars, Modnames = mod, var.expl)) %>% 
  select(Modnames, cr.age, sexM, SEASON1314, K, Delta_AICc, AICcWt, LL, var.expl) %>% 
  mutate(var.expl = round(var.expl, 1),
         across(c(Delta_AICc, AICcWt), ~round(., 2)))


saveRDS(tib.coef_aic, "rds/cr_tib_best")





