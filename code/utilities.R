






# define some helper functions ----
# generic function for runing CJS models
run.models<-function(phi.stru, p.stru = "SEASON * sex + time") {
  phi.stru.list <- list(formula = formula(paste("~", phi.stru)))
  p.stru.list <- list(formula = formula(paste("~", p.stru)))
zmod <- mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru.list, p = p.stru.list), output = FALSE, chat = zchat)
return(zmod)
}


# Arnold 2010 recommends removing uninformative parameters if 85% CI overlap 0. 85% is more consistent with AIC selection (95% CI is too broad and we end up excluding variables based on different standards). RMark doesn't provide 85% CI, so here I calculate them based on multiplying se by 1.44 (from Z table for 85% CI). 
# To help ensure that multiplying se by some Z value provides appropriate CI, I also include a test for whether 95% CI calculated by this method are equal to 95% provided by R/RMark.

# linear model version
lm_parm_informative <- function(zmodel, zmod.name) {
 # zmod.name = deparse(substitute(zmodel)) 
 #zmod.name = names(zmodel)  
 #zmodel <- zmodel[[1]]
  zmodel_ci <- confint(zmodel)%>% 
                data.frame() %>% 
                rownames_to_column("parm") %>% 
    rename(lcl = X2.5..,
           ucl = X97.5..)
  
mod_call <- as.character(zmodel$call[2])
mod_call <- sub(".*\\~ ", "", mod_call)
  
parms_informative <- coef(summary(zmodel)) %>% 
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
           mod.name = zmod.name,
           #mod.name = sub(".*\\$ ", "", mod.name),
           mod.call = mod_call)
return(parms_informative)
}



cjs_parm_informative <- function(zmodel, phi.p) {
#  zmod.name = deparse(substitute(zmodel))
  
  model_name <- as.character(zmodel$model.name[1])
model_name <- sub(".*\\~ ", "", model_name)
  
  
  coef(zmodel) %>% 
  data.frame() %>% 
  rownames_to_column("parm") %>% 
  filter(grepl(paste(phi.p, ":", sep = ""), parm), !grepl("ntercept", parm), !grepl("time", parm)) %>% 
    mutate(lci85 = estimate - (1.44 * se),
           uci85 = estimate + (1.44 * se),
           ci95check = lcl == estimate - (1.96 * se) & ucl == estimate + (1.96 * se),
           informative95 = ifelse((lcl < 0 & ucl < 0) | (lcl > 0 & ucl > 0), TRUE, FALSE),
           informative85 = ifelse((lci85 < 0 & uci85 < 0) | (lci85 > 0 & uci85 > 0), TRUE, FALSE),
          # mod.name = zmod.name,
           model.name = model_name)
}
#

get_mod_call <- function(zmodel) {

  zmod.name = deparse(substitute(zmodel)) 
zmod.name <- sub(".*\\$ ", "", zmod.name)
  
  
mod_call <- as.character(zmodel$call[2])
mod_call <- sub(".*\\~ ", "", mod_call)
  
mod_call_df <- data.frame(mod.call = mod_call,
                          mod.name = zmod.name)
  return(mod_call_df)
}



#' Change variable names in model call to nicer versions for output in model tables
#' 
#'  Adds a nice 
#'
#' @param zmod.call column in data frame with the results of mod
#'
#' @return
#' @export
#'
#' @examples
mod_call_to_structure <- function(zmod.call) {
  zmod.call = gsub("SEASON", "Yr", zmod.call)
  zmod.call = gsub("sex", "Sex", zmod.call)
  zmod.call = gsub("hatch", "Hatch", zmod.call)
  zmod.call = gsub("weight.slope40", "Mass", zmod.call)
  zmod.call = gsub("tibiotar.slope35", "Tibio", zmod.call)
  zmod.call = gsub("flipper.slope40", "Flipper", zmod.call)
  zmod.call = gsub("cr.age", paste("Cr", "\u00E8", "ching age", sep = ""), zmod.call)
}
