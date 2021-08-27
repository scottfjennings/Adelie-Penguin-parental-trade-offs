# thesis_ch4_analysis
code, data, etc. for analysis of Adelie Penguin survival; chapter 4 of thesis



The analysis has 2 major parts, each with multiple sub-parts

1. Growth rates and creche timing
1) how do mass, flipper, and tibiotarsus growth rates predict age at crèche; 
2) how does mass growth rate predict mass at creching;
3) how does flipper growth rate predict flipper length at creching; and 
4) how does tibiotarsus growth rate predict tibiotarsus length at creching

2. Survival
1) can survival across the entire provisioning period be predicted by mass, flipper or tibiotarsus growth rates; and 
2) can survival during the crèche period be predicted by age or size at crèching



Part 1 is accomplished with a series of linear mixed effects models via the code in the following files:
CH4_growth_creche.R

Part 2 is accomplished with Cormack-Jolly-Seber (CJS) survival models via the code in the following files:
Rmark_penguins.R
Rmark_penguins_sex_groups.R





8/27/21: there are 13 chicks in the new inp that were excluded from the thesis analysis. I have examined here("data/ana_table_full") and can't find anything that separates these chicks from all the others. Also searched through all the folders in thesis_data_work and couldn't find any notes indicating why these chicks were excluded. I am currently running the analysis with these chicks.