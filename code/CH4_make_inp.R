


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
