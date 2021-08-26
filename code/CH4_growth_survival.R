

#need to tell R that Mark is not stored in C:\Program Files


# packages loaded in CH4_prep_survival_data

run.phi.model<-function(phi.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru, p = list(formula = ~SEASON + time)), output = FALSE, chat = 1.25)
}


# evaluating effect of growth rates on survival
# comparing models with growth rates only to those with growth rates and the best structure form the final step of overall survival modeling:
# ~SEASON * sex + res.htch + Time			 	
Phi.massgr = run.phi.model(list(formula = ~weight))
Phi.flipgr = run.phi.model(list(formula = ~flipper))
Phi.tibgr = run.phi.model(list(formula = ~tibiotar))
Phi.SEASON.sex_hatch_T_massgr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + weight))
Phi.SEASON.sex_hatch_T_flipgr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + flipper))
Phi.SEASON.sex_hatch_T_tibgr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + tibiotar))
Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~SEASON * sex + res.htch + Time)) #best from overall survival
Phi.dot = run.phi.model(list(formula = ~1))

growth_surv_models = collect.models(lx=c("Phi.massgr", "Phi.flipgr", "Phi.tibgr", "Phi.SEASON.sex_hatch_T_massgr", "Phi.SEASON.sex_hatch_T_flipgr", "Phi.SEASON.sex_hatch_T_tibgr", "Phi.SEASON.sex_hatch_T", "Phi.dot"))


saveRDS(growth_surv_models, here("fitted_models/growth_surv_models"))


growth_surv_models <- readRDS(here("fitted_models/growth_surv_models"))
