

library(tidyverse)
library(here)
library(lme4)
library(AICcmodavg)
library(grid)
options(scipen = 999)

# read data ----
#  ana.table.full.csv is created by the sequence of numbered code files in THESIS\thesis_data_work\code_files\analysis_data_prep
#ana.table=read.csv("data/ana.table.full.csv")
ana.table <- readRDS(here("data/ana_table_full"))
#ana.table.old=read.csv("data/ana.table.full.old.csv")%>%  rename(CHICK_ID = id)
# data management ----

data = ana.table %>% 
  filter(!is.na(cr.seasday), !is.na(weight.slope40), !is.na(sex)) %>% 
  mutate(SEASON = as.factor(SEASON))


# analysis ----
# first creching age as response variable ----
# creching age model selection step 1 ----

cr_age_step1_mods <- list(
sex = lm(cr.age~	sex, data = data),
year = lm(cr.age~	SEASON, data = data),
sex_year = lm(cr.age~	sex+SEASON, data = data),
sex.year = lm(cr.age~	sex*SEASON, data = data)
)

saveRDS(cr_age_step1_mods, here("fitted_models/cr_age_step1_mods"))

aictab(list(cr_age_step1_mods$sex, cr_age_step1_mods$year, cr_age_step1_mods$sex_year, cr_age_step1_mods$sex.year), c("sex", "year", "sex_year", "sex.year"))

# year model best supported

# creching age model selection step 2 ----

cr_age_step2_mods <- list(
mass_year	= lm(cr.age ~	weight.slope40 + SEASON, data = data),
flip_year	= lm(cr.age ~	flipper.slope40	+ SEASON, data = data),
tib_year	= lm(cr.age ~	tibiotar.slope35 +	SEASON, data = data),
hatch_year	= lm(cr.age ~	resid.hatch + SEASON, data = data),
year	= lm(cr.age ~	SEASON, data = data),
mass_hatch_year	= lm(cr.age ~	weight.slope40 + resid.hatch +	SEASON, data = data),
flip_hatch_year	= lm(cr.age ~	flipper.slope40 + resid.hatch +	SEASON, data = data),
tib_hatch_year = lm(cr.age ~	tibiotar.slope35 + resid.hatch +	SEASON, data = data),
int	= lm(cr.age~ 1, data = data)
)

aictab(cr_age_step2_mods, names(cr_age_step2_mods))

saveRDS(cr_age_step2_mods, here("fitted_models/cr_age_step2_mods"))


# all models with resid.hatch are competitive (largest dAICc = 2.14), then a break in dAICc to 43.18 for 5th best.

# creching mass model selection step 1 ----

cr_mass_step1_mods <- list(
mass.sex	= lm(cr.mass~	sex, data = data),
mass.year	= lm(cr.mass~	SEASON, data = data),
mass.sex_year	= lm(cr.mass~	sex+SEASON, data = data),
mass.sex.year	= lm(cr.mass~	sex*SEASON, data = data)
)


(aictab(list(cr_mass_step1_mods$mass.sex, cr_mass_step1_mods$mass.year, cr_mass_step1_mods$mass.sex_year, cr_mass_step1_mods$mass.sex.year), c("mass.sex", "mass.year", "mass.sex_year", "mass.sex.year")))

# sex_year best
saveRDS(cr_mass_step1_mods, "fitted_models/cr_mass_step1_mods")


# mass.sex_year best supported
# creching mass model selection step 2 ----
cr_mass_step2_mods = list(
mass.growth_sex_year = lm(cr.mass ~	weight.slope40 +	sex + SEASON, data = data),
mass.sex_year = lm(cr.mass ~	sex + SEASON, data = data),
mass.hatch_sex_year = lm(cr.mass ~	resid.hatch +	sex + SEASON, data = data),
mass.growth_hatch_sex_year	= lm(cr.mass ~	weight.slope40 + resid.hatch + sex + SEASON, data = data),
mass.age_sex_year = lm(cr.mass ~	cr.age + sex + SEASON, data = data),
mass.age_growth_sex_year	= lm(cr.mass ~	cr.age + weight.slope40 + sex + SEASON, data = data),
mass.int	= lm(cr.mass ~	1, data = data)
)

aictab(list(cr_mass_step2_mods$mass.growth_sex_year, cr_mass_step2_mods$mass.sex_year, cr_mass_step2_mods$mass.hatch_sex_year, cr_mass_step2_mods$mass.growth_hatch_sex_year, cr_mass_step2_mods$mass.age_sex_year, cr_mass_step2_mods$mass.age_growth_sex_year, cr_mass_step2_mods$mass.int), c("mass.growth_sex_year", "mass.sex_year", "mass.hatch_sex_year", "mass.growth_hatch_sex_year", "mass.age_sex_year", "mass.age_growth_sex_year", "mass.int"))

saveRDS(cr_mass_step2_mods, here("fitted_models/cr_mass_step2_mods"))

# creching flipper length model selection step 1 ----

cr_flip_step1_mods <- list(
flip.sex = lm(cr.flip~	sex, data = data),
flip.year	= lm(cr.flip~	SEASON, data = data),
flip.sex_year	= lm(cr.flip~	sex+SEASON, data = data),
flip.sex.year	= lm(cr.flip~	sex*SEASON, data = data)
)

aictab(list(cr_flip_step1_mods$flip.sex, cr_flip_step1_mods$flip.year, cr_flip_step1_mods$flip.sex_year, cr_flip_step1_mods$flip.sex.year), c("flip.sex", "flip.year", "flip.sex_year", "flip.sex.year"))

# flip.year best
saveRDS(cr_flip_step1_mods, "fitted_models/cr_flip_step1_mods")
#

# creching flipper model selection step 2 ----
cr_flip_step2_mods = list(
flip.growth_year = lm(cr.flip ~	flipper.slope40 +	SEASON, data = data),
flip.year = lm(cr.flip ~	sex + SEASON, data = data),
flip.hatch_year = lm(cr.flip ~	resid.hatch + SEASON, data = data),
flip.growth_hatch_year	= lm(cr.flip ~	flipper.slope40 + resid.hatch + SEASON, data = data),
flip.age_year = lm(cr.flip ~	cr.age + SEASON, data = data),
flip.age_growth_year	= lm(cr.flip ~	cr.age + flipper.slope40 + SEASON, data = data),
flip.int	= lm(cr.flip ~	1, data = data)
)

aictab(list(cr_flip_step2_mods$flip.growth_year, cr_flip_step2_mods$flip.year, cr_flip_step2_mods$flip.hatch_year, cr_flip_step2_mods$flip.growth_hatch_year, cr_flip_step2_mods$flip.age_year, cr_flip_step2_mods$flip.age_growth_year, cr_flip_step2_mods$flip.int), c("flip.growth_year", "flip.year", "flip.hatch_year", "flip.growth_hatch_year", "flip.age_year", "flip.age_growth_year", "flip.int"))

saveRDS(cr_flip_step2_mods, "fitted_models/cr_flip_step2_mods")

# creching tib length model selection step 1 ----

cr_tib_step1_mods <- list(
tib.sex	= lm(cr.tib~	sex, data = data),
tib.year = lm(cr.tib~	SEASON, data = data),
tib.sex_year = lm(cr.tib~	sex+SEASON, data = data),
tib.sex.year = lm(cr.tib~	sex*SEASON, data = data)
)

(aictab(list(cr_tib_step1_mods$tib.sex, cr_tib_step1_mods$tib.year, cr_tib_step1_mods$tib.sex_year, cr_tib_step1_mods$tib.sex.year), c("tib.sex", "tib.year", "tib.sex_year", "tib.sex.year")))

# tib.year best
saveRDS(cr_tib_step1_mods, "fitted_models/cr_tib_step1_mods")

# creching tibiotarsus model selection step 2 ----
cr_tib_step2_mods = list(
tib.growth_year = lm(cr.tib ~	tibiotar.slope35 +	SEASON, data = data),
tib.year = lm(cr.tib ~	sex + SEASON, data = data),
tib.hatch_year = lm(cr.tib ~	resid.hatch + SEASON, data = data),
tib.growth_hatch_year	= lm(cr.tib ~	tibiotar.slope35 + resid.hatch + SEASON, data = data),
tib.age_year = lm(cr.tib ~	cr.age + SEASON, data = data),
tib.age_growth_year	= lm(cr.tib ~	cr.age + tibiotar.slope35 + SEASON, data = data),
tib.int	= lm(cr.tib ~	1, data = data)
)

aictab(list(cr_tib_step2_mods$tib.growth_year, cr_tib_step2_mods$tib.year, cr_tib_step2_mods$tib.hatch_year, cr_tib_step2_mods$tib.growth_hatch_year, cr_tib_step2_mods$tib.age_year, cr_tib_step2_mods$tib.age_growth_year, cr_tib_step2_mods$tib.int), c("tib.growth_year", "tib.year", "tib.hatch_year", "tib.growth_hatch_year", "tib.age_year", "tib.age_growth_year", "tib.int"))

saveRDS(cr_tib_step2_mods, "fitted_models/cr_tib_step2_mods")

