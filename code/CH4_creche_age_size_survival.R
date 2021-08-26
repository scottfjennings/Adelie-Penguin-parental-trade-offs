


#need to tell R that Mark is not stored in C:\Program Files


# packages loaded in CH4_prep_survival_data

run.phi.model<-function(phi.stru) {
mark(penguins.process, penguins.ddl, model.parameters = list(Phi = phi.stru, p = list(formula = ~SEASON + time)), output = FALSE, chat = 1.25)
}

# how does the size and age at which chicks enter the creche stage relate to survival during creche stage?
# 7/9/21 there are a bunch of models here that I'm not sure should be included. I can't find any documentation why these models were added, and they don't appear in any draft. commenting them out for now, only keeping the ones that are in all the ms drafts

#Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~SEASON * sex + res.htch + Time)) #best from overall survival; not sure it make sense for this to be here; same data?
# Phi.dot	 =	run.phi.model(list(formula = ~1))
					
Phi.cr.mass	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.mass)	) # cr mass
Phi.cr.flip	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.flip)	) # cr flip
Phi.cr.tib	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.tib)	) # cr tib
Phi.fate.age	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age)	) # cr age

Phi.SEASON.sex_hatch_T = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time)) # best overall survival
Phi.SEASON.sex_hatch_T_tibgr = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + tibiotar)) # best growth survival


					
# Phi.fate.age_cr.mass	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.mass)	) # cr age + cr mass
# Phi.fate.age_cr.flip	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.flip)	) # cr age + cr flip
# Phi.fate.age_cr.tib	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.tib)	) # cr age + cr tib
					
Phi.SEASON.sex_hatch_T_crmass = run.phi.model(list(formula =	~ SEASON * sex + res.htch + Time + did.cr + did.cr:cr.mass)) # cr mass + best overall
Phi.SEASON.sex_hatch_T_crflip =	run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + did.cr + did.cr:cr.flip)) # cr flip + best overall
Phi.SEASON.sex_hatch_T_crtib = run.phi.model(list(formula = ~ SEASON * sex + res.htch + Time + did.cr + did.cr:cr.tib)) # cr tib + best overall
Phi.SEASON.sex_hatch_T_fateage = run.phi.model(list(formula =	~ SEASON * sex + res.htch + Time + did.cr + did.cr:fate.age)) # cr age + best overall
					
# Phi.fate.age_cr.mass_T_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.mass+Time+res.htch+SEASON*sex)	) # cr mass + cr age + best overall
# Phi.fate.age_cr.flip_T_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.flip+Time+res.htch+SEASON*sex)	) # cr flip + cr age + best overall
# Phi.fate.age_cr.tib_T_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:fate.age+did.cr:cr.tib+Time+res.htch+SEASON*sex)	) # cr tib + cr age + best overall
					

##-- take out T from best models
# Phi.cr.mass_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.mass+res.htch+SEASON*sex)	)
# Phi.cr.flip_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.flip+res.htch+SEASON*sex)	)
# Phi.cr.tib_hatch_y.s	 =	run.phi.model(list(formula	 =	~did.cr+did.cr:cr.tib+res.htch+SEASON*sex)	)


cr_age_size_models = collect.models(lx=c("Phi.cr.mass", "Phi.cr.flip", "Phi.cr.tib", "Phi.fate.age", "Phi.SEASON.sex_hatch_T", "Phi.SEASON.sex_hatch_T_tibgr", "Phi.SEASON.sex_hatch_T_crmass", "Phi.SEASON.sex_hatch_T_crtib", "Phi.SEASON.sex_hatch_T_crflip"))

saveRDS(cr_age_size_models, here("fitted_models/cr_age_size_models"))

