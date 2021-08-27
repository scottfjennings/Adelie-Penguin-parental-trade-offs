


# create csv for conversion to inp
# according to notes this was originally created by hand in a csv
# need to paste capture history, sex, year, time-varying age, and individual covariates
# sex year groups are indicated by 4 dummy variables in the first 4 spaces after the capture history 1 = m1213, 2 = f1213, 3 = m1314, 4 = f1314
# capture history is in thesis_data_work/ChickGrowth2015/resight.tab.all.csv
# sex, year and individual covs are in thesis_data_work/rds/ana_table_full
# both of these copied to here("data")
# time-varying age needs to be calculated for season days 50-98

library(tidyverse)
library(here)

 
ch <- read.csv(here("data/resight.tab.all.csv")) %>% 
  select(-seen) %>% 
  unite(ch, -contains("chick_id"), sep = "")



ana_table <- readRDS(here("data/ana_table_full")) %>% 
  rename(chick_id = CHICK_ID) %>% 
  filter(fail.age >= 10 | is.na(fail.age), !is.na(sex))

time_var_age <- expand.grid(chick_id = ana_table$chick_id,
                            seas.day = seq(50, 98)) %>% 
  full_join(., ana_table %>% 
              select(chick_id, hatch.seasday, cr.seasday)) %>% 
  mutate(dayold = seas.day - hatch.seasday + 1, # need age to be 1 on hatch day
         dayold = ifelse(dayold < 0, 0, dayold),
         in.cr = ifelse(seas.day >= cr.seasday, 1, 0),
         in.cr = ifelse(is.na(in.cr), 0, in.cr)) %>% 
  arrange(chick_id, seas.day)

time_var_age_wide <- time_var_age %>%
  data.frame() %>%
  pivot_longer(cols = c("dayold", "in.cr")) %>% 
  mutate(col.name = paste(name, seas.day, sep = "")) %>% 
  pivot_wider(id_cols = chick_id, names_from = col.name, values_from = value) %>% 
  select(chick_id, contains("dayold"), contains("in.cr")) %>% 
  unite(ages, -contains("chick_id"), sep = " ")


indiv_covs <- ana_table %>% 
  left_join(time_var_age_wide) %>% 
  mutate(m1213 = ifelse(SEASON == 1213 & sex == "M", 1, 0),
         f1213 = ifelse(SEASON == 1213 & sex == "F", 1, 0),
         m1314 = ifelse(SEASON == 1314 & sex == "M", 1, 0),
         f1314 = ifelse(SEASON == 1314 & sex == "F", 1, 0)) %>% 
  select(chick_id, m1213, f1213, m1314, f1314, ages, cr.age, fail.age, cr.mass, 
					resid.hatch, cr.flip, cr.tib, fail.age) %>% 
  mutate(resid.hatch = round(resid.hatch, 1),
         cr.mass = round(cr.mass, 2),
         cr.flip = round(cr.flip, 2),
         cr.tib = round(cr.tib, 2)) %>% 
  unite(indiv.covs, -contains("chick_id"), sep = " ")

mark_inp <- inner_join(ch, indiv_covs)  %>% 
  distinct() %>%  
  mutate(chick_id = paste("/*", chick_id, "*/  ")) %>% 
  unite(ch, sep = " ") %>% 
  mutate(ch = paste(ch, ";", sep = ""))




write.csv(mark_inp, here("data/mark_in.csv"), row.names = FALSE)

# open the resulting file in notepad++, delete the first line that has only "ch", then do replace all to get rid of all the ".
# finally, do save as and include ".inp" to the end of the file name to save as .inp

##############################

# for comparison, read in inp file used for thesis analysis

# 8/27/21: there are 13 chicks in the new inp that were excluded from the thesis analysis. I have examined here("data/ana_table_full") and can't find anything that separates these chicks from all the others. Also searched through all the folders in thesis_data_work and couldn't find any notes indicating why these chicks were excluded. I am currently running the analysis with these chicks.

penguins_thesis = convert.inp(here("data/all_covs_sex_groups.inp"), group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56",  "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64",  "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72", "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80", "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96", "dayold97", "dayold98", "av.food", "av.trip.length", "cr.age", "fail.age", "weight", "flipper", "tibiotar", "CH_B", "CH_S", "mean.N", "cr.mass", "res.htch", "cr.flip", "cr.tib", "did.cr", "fate.age" ),use.comments=FALSE) 
 

penguins=convert.inp(here("data/mark_in.inp"), 
					group.df=data.frame(sex=rep(c("Male","Female"),2), SEASON=c(rep("1213",2),rep("1314",2))), 
					covariates=c("dayold50", "dayold51", "dayold52", "dayold53", "dayold54", "dayold55", "dayold56", 
					 "dayold57", "dayold58", "dayold59", "dayold60", "dayold61", "dayold62", "dayold63", "dayold64", 
					 "dayold65", "dayold66", "dayold67", "dayold68", "dayold69", "dayold70", "dayold71", "dayold72",
					 "dayold73", "dayold74", "dayold75", "dayold76", "dayold77", "dayold78", "dayold79", "dayold80",
					 "dayold81", "dayold82", "dayold83", "dayold84", "dayold85", "dayold86", "dayold87", "dayold88", 
					 "dayold89", "dayold90", "dayold91", "dayold92", "dayold93", "dayold94", "dayold95", "dayold96",
					 "dayold97", "dayold98","in.cr50", "in.cr51", "in.cr52", "in.cr53", "in.cr54", "in.cr55", "in.cr56", 
					 "in.cr57", "in.cr58", "in.cr59", "in.cr60", "in.cr61", "in.cr62", "in.cr63", "in.cr64", 
					 "in.cr65", "in.cr66", "in.cr67", "in.cr68", "in.cr69", "in.cr70", "in.cr71", "in.cr72",
					 "in.cr73", "in.cr74", "in.cr75", "in.cr76", "in.cr77", "in.cr78", "in.cr79", "in.cr80",
					 "in.cr81", "in.cr82", "in.cr83", "in.cr84", "in.cr85", "in.cr86", "in.cr87", "in.cr88", 
					 "in.cr89", "in.cr90", "in.cr91", "in.cr92", "in.cr93", "in.cr94", "in.cr95", "in.cr96",
					 "in.cr97", "in.cr98", "cr.age", "fail.age", "cr.mass", "res.htch", "cr.flip", "cr.tib"), use.comments = TRUE) 

peng_comp <- anti_join(penguins %>% rownames_to_column("chick_id") %>% select(chick_id, ch),
                       penguins_thesis %>% select(ch)) %>% 
  mutate(excluded = TRUE,
         chick_id = trimws(chick_id)) %>% 
  select(-ch) %>% 
  full_join(ana_table)


