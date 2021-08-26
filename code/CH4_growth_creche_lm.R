

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


# analysis ----
# first creching age as response variable ----
# creching age model selection step 1 ----


sex.lm	<- lm(cr.age~	sex, data = data)
year.lm	<- lm(cr.age~	SEASON, data = data)
sex_year.lm	<- lm(cr.age~	sex+SEASON, data = data)
sex.year.lm	<- lm(cr.age~	sex*SEASON, data = data)

(cr_age_aic.lm1 <- aictab(list(sex.lm, year.lm, sex_year.lm, sex.year.lm), c("sex.lm", "year.lm", "sex_year.lm", "sex.year.lm")))

saveRDS(cr_age_aic.lm1, "rds/cr_age_aic_lm1")


# creching age model selection step 2 ----
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



# creching age rsquared for lms ----
lm.adj.r2 <- rbind(data.frame(mod = "hatch_year.lm", adj.r2 = summary(hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "mass_hatch_year.lm", adj.r2 = summary(mass_hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "tib_hatch_year.lm", adj.r2 = summary(tib_hatch_year.lm)$adj.r.squared),
                   data.frame(mod = "flip_hatch_year.lm", adj.r2 = summary(flip_hatch_year.lm)$adj.r.squared))

# creching age combine coefficients, aic variables, variance explained ----
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




  
# creching mass model selection step 1 ----


mass.sex.lm	<- lm(cr.mass~	sex, data = data)
mass.year.lm	<- lm(cr.mass~	SEASON, data = data)
mass.sex_year.lm	<- lm(cr.mass~	sex+SEASON, data = data)
mass.sex.year.lm	<- lm(cr.mass~	sex*SEASON, data = data)

(cr_mass_aic.lm1 <- aictab(list(mass.sex.lm, mass.year.lm, mass.sex_year.lm, mass.sex.year.lm), c("mass.sex.lm", "mass.year.lm", "mass.sex_year.lm", "mass.sex.year.lm")))

# sex_year best
saveRDS(cr_mass_aic.lm1, "rds/cr_mass_aic_lm1")


# mass.sex_year best supported

# creching flipper length model selection step 1 ----

flip.sex.lm	<- lm(cr.flip~	sex, data = data)
flip.year.lm	<- lm(cr.flip~	SEASON, data = data)
flip.sex_year.lm	<- lm(cr.flip~	sex+SEASON, data = data)
flip.sex.year.lm	<- lm(cr.flip~	sex*SEASON, data = data)

(cr_flip_aic.lm1 <- aictab(list(flip.sex.lm, flip.year.lm, flip.sex_year.lm, flip.sex.year.lm), c("flip.sex.lm", "flip.year.lm", "flip.sex_year.lm", "flip.sex.year.lm")))

# sex_year best
saveRDS(cr_flip_aic.lm1, "rds/cr_flip_aic_lm1")

#creching tib length model selection step 1 ----

tib.sex.lm	<- lm(cr.tib~	sex, data = data)
tib.year.lm	<- lm(cr.tib~	SEASON, data = data)
tib.sex_year.lm	<- lm(cr.tib~	sex+SEASON, data = data)
tib.sex.year.lm	<- lm(cr.tib~	sex*SEASON, data = data)

(cr_tib_aic.lm1 <- aictab(list(tib.sex.lm, tib.year.lm, tib.sex_year.lm, tib.sex.year.lm), c("tib.sex.lm", "tib.year.lm", "tib.sex_year.lm", "tib.sex.year.lm")))

# sex_year best
saveRDS(cr_tib_aic.lm1, "rds/cr_tib_aic_lm1")


# creche size model selection step 2  ----
# same best structure for all morph in step 1, so can do step 2 with same function for all morphs
cr_size_step2_cand_set_sex_year <- function(data) {

growth_sex_year	<- lm(cr.size~	growth.rate		 +	sex+SEASON		 , data = data)
sex_year	<- lm(cr.size~	sex+SEASON		 , data = data)
hatch_sex_year	<- lm(cr.size~	resid.hatch+	sex+SEASON		 , data = data)
growth_hatch_sex_year	<- lm(cr.size~	growth.rate		 +resid.hatch+	sex+SEASON		 , data = data)
age_sex_year	<- lm(cr.size~	cr.age +	sex+SEASON		 , data = data)
age_growth_sex_year	<- lm(cr.size~	cr.age + growth.rate		 +	sex+SEASON		 , data = data)
int	<- lm(cr.size~	1		 , data = data)

cr_size_aic <- aictab(list(growth_sex_year, sex_year, hatch_sex_year, growth_hatch_sex_year, age_sex_year, age_growth_sex_year, int), c("growth_sex_year", "sex_year", "hatch_sex_year", "growth_hatch_sex_year", "age_sex_year", "age_growth_sex_year", "int")) %>% 
  data.frame() 

return(cr_size_aic)
}

cr_mass_aic <- data %>%  
      select(cr.size = cr.mass, growth.rate = weight.slope40, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_sex_year() 
cr_mass_aic
# age_growth_sex_year best by 21.17 aicc

saveRDS(cr_mass_aic, "rds/cr_mass_aic_lm")

cr_flip_aic <- data %>% 
      select(cr.size = cr.flip, growth.rate = flipper.slope40, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_sex_year() 
# age_growth_sex_year best by 11.8

saveRDS(cr_flip_aic, "rds/cr_flip_aic_lm")



cr_tib_aic <- data %>% 
      select(cr.size = cr.tib, growth.rate = tibiotar.slope35, sex, SEASON, resid.hatch, nest.seas, cr.age) %>% 
      cr_size_step2_cand_set_sex_year() 
# age_sex_year best
# age_growth_sex_year daicc = 0.64 but uninformative. next best daicc = 13.16

saveRDS(cr_tib_aic, "rds/cr_tib_aic")


# best models output ----
# creching mass best model output ----

mass.age_growth_sex_year	<- lm(cr.mass~	cr.age + weight.slope40 +	sex + SEASON, data = data)
mass.int	<- lm(cr.mass~	1, data = data)

# creching age rsquared for lms ----
mass.adj.r2 <- rbind(data.frame(Modnames = "mass.age_growth_sex_year", adj.r2 = summary(mass.age_growth_sex_year)$adj.r.squared))

# coefficients
mass.age_growth_sex_year_coefCI <- cbind(coef(mass.age_growth_sex_year), confint(mass.age_growth_sex_year)) %>% 
  data.frame() %>% 
  rename(est = 1, lci = 2, uci = 3) %>% 
  rownames_to_column("varb") %>% 
  mutate(mod = "age_growth_sex_year.lm")

mass.coefCI_wider <- mass.age_growth_sex_year_coefCI %>% 
  filter(varb != "(Intercept)") %>% 
  mutate(across(c(est, lci, uci), ~round(., 2))) %>% 
  mutate(coef.ci = paste(est, " (", lci, " to ", uci, ")", sep = "")) %>% 
  select(mod, varb, coef.ci) %>% 
  pivot_wider(id_cols = mod, names_from = varb, values_from = coef.ci) 


mass.coef_aic <- mass.coefCI_wider %>% 
  rename(Modnames = mod) %>% 
  left_join(., cr_mass_aic) %>% 
  left_join(mass.adj.r2) %>% 
  select(Modnames, cr.age, weight.slope40, sexM, SEASON1314, K, Delta_AICc, AICcWt, LL, var.expl) %>% 
  mutate(var.expl = round(var.expl, 1),
         across(c(Delta_AICc, AICcWt), ~round(., 2)))

saveRDS(mass.coef_aic, "rds/cr_mass_best")


##--
# creching flipper length best model output ----

flip.age_growth_year	<- lm(cr.flip ~	cr.age + flipper.slope40 + SEASON +	(1 | nest.seas), data=data, REML = FALSE)
flip.int	<- lm(cr.flip ~	1 , data = data)
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

tib.age_sex_year	<- lm(cr.tib ~ cr.age + sex + SEASON , data = data)
tib.int	<- lm(cr.tib ~	1		 , data = data)
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





