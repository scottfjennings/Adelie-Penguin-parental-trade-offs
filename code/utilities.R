






# define some helper functions ----
# generic function for runing CJS models
run.models<-function(phi.stru, p.stru = "SEASON * sex + time") {
  phi.stru.list <- list(formula = formula(paste("~", phi.stru)))
  p.stru.list <- list(formula = formula(paste("~", p.stru)))
zmod <- mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru.list, p = p.stru.list), output = FALSE, chat = 1.25)
return(zmod)
}


# Arnold 2010 recommends removing uninformative parameters if 85% CI overlap 0. 85% is more consistent with AIC selection (95% CI is too broad and we end up excluding variables based on different standards). RMark doesn't provide 85% CI, so here I calculate them based on multiplying se by 1.44 (from Z table for 85% CI). 
# To help ensure that multiplying se by some Z value provides appropriate CI, I also include a test for whether 95% CI calculated by this method are equal to 95% provided by R/RMark.

# linear model version
lm_parm_informative <- function(zmodel) {
  zmod.name = deparse(substitute(zmodel))
  
  zmodel_ci <- confint(zmodel)%>% 
                data.frame() %>% 
                rownames_to_column("parm") %>% 
    rename(lcl = X2.5..,
           ucl = X97.5..)
  
  coef(summary(zmodel)) %>% 
  data.frame() %>% 
  rownames_to_column("parm") %>% 
  filter(!grepl("ntercept", parm)) %>% 
    select(parm, Estimate, Std..Error) %>% 
    left_join(., zmodel_ci) %>% 
    mutate(lci85 = Estimate - (1.44 * Std..Error),
           uci85 = Estimate + (1.44 * Std..Error),
           ci95check = lcl == Estimate - (1.96 * Std..Error) & ucl == Estimate + (1.96 * Std..Error),
           informative95 = ifelse((lcl < 0 & ucl < 0) | (lcl > 0 & ucl > 0), TRUE, FALSE),
           informative85 = ifelse((lci85 < 0 & uci85 < 0) | (lci85 > 0 & uci85 > 0), TRUE, FALSE),
           mod.name = zmod.name)
}




cjs_parm_informative <- function(zmodel, phi.p) {
  zmod.name = deparse(substitute(zmodel))
  
  coef(zmodel) %>% 
  data.frame() %>% 
  rownames_to_column("parm") %>% 
  filter(grepl(paste(phi.p, ":", sep = ""), parm), !grepl("ntercept", parm), !grepl("time", parm)) %>% 
    mutate(lci85 = estimate - (1.44 * se),
           uci85 = estimate + (1.44 * se),
           ci95check = lcl == estimate - (1.96 * se) & ucl == estimate + (1.96 * se),
           informative95 = ifelse((lcl < 0 & ucl < 0) | (lcl > 0 & ucl > 0), TRUE, FALSE),
           informative85 = ifelse((lci85 < 0 & uci85 < 0) | (lci85 > 0 & uci85 > 0), TRUE, FALSE),
           mod.name = zmod.name)
}
#