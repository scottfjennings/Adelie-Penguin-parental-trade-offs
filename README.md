# thesis_ch4_analysis
code, data, etc. for analysis of Adelie Penguin survival; chapter 4 of thesis



The analysis has 2 major parts, each with multiple sub-parts

1. Growth rates and creche timing
1) Do mass, flipper, and tibiotarsus growth rates predict age at crèche; 
2) Do mass, flipper or tibiotarsus growth rates predict mass, flipper length, or tibiotarsus length at creching?

2. Survival
1) how do daily survival rates change through the chick-rearing period?
2) can survival across the entire provisioning period be predicted by mass, flipper or tibiotarsus growth rates? and 
2) can survival during the crèche period be predicted by age or size at crèching?



Part 1 is accomplished with a series of linear mixed effects models via the code in the following files:
CH4_growth_creche.R

Part 2 is accomplished with Cormack-Jolly-Seber (CJS) survival models via the code in the following files:
Rmark_penguins.R
Rmark_penguins_sex_groups.R





8/27/21: there are 13 chicks in the new inp that were excluded from the thesis analysis. These chicks were excluded because they died before being measured at least 2 times, so I couldn't calculate growth rates. However, growth rates are not considered in any candidate models for estimating overall survival (2.1 above) 