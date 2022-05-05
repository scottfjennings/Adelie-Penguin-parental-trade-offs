# thesis_ch4_analysis
code, data, etc. for analysis of Adelie Penguin survival; chapter 4 of thesis



The analysis has 2 major parts, each with multiple sub-parts

1. Growth rates and creche timing
1) Do mass, flipper, and tibiotarsus growth rates predict age at crèche; 
2) Do mass, flipper or tibiotarsus growth rates predict mass, flipper length, or tibiotarsus length at creching, respectively?

2. Survival
1) how do daily survival rates change through the chick-rearing period?
2) can survival across the entire provisioning period be predicted by mass, flipper or tibiotarsus growth rates? and 
2) can survival during the crèche period be predicted by age or size at crèching?



Part 1 is accomplished with a series of linear models via the code in the following files:
creche_timing1_linear_models.R
  * see 9/1/21 note for the file creche_timing1_mixed_models.R 

Part 2 is accomplished with Cormack-Jolly-Seber (CJS) survival models via the code in the following files:
survival1_make_inp.R
survival2_GOF_testing.R
survival3_overall_survival.R
survival4_growth_effect.R
survival5_creche_size_age_effect

All output for the paper is created by the following code/RMD files
output1_model_selection_tables.RMD
output2_creche_timing_coef_tables.Rmd
output3_survival_figures.R



Some data management and analysis notes (reverse chrono order):
9/8/21: At the first stage of modeling Phi for overall survival, there is a natural break in model support between Phi.SEASON.sex_hatch (dQAICc = 1.14, 2nd model) and Phi.SEASON.hatch (dQAICc = 1.99, 3rd ranked). In contrast there isn't a particularly meaningful break in model support spanning dQAICc = 2: Phi.SEASON.hatch dQAICc = 1.99; Phi.SEASON_sex_hatch dQAICc = 2.16. There is also a natural break in model support between Phi.SEASON_hatch2 (dQAICc = 2.18, 5th ranked) and Phi.SEASON.sex_hatch2 (dQAICc = 3.31, 6th ranked), but there are 9 possible new structures in the next step of Phi modeling, so passing all 5 structures from this step forward will make a very large candidate set at the next step, and all these 5 structures are quite similar to each other. So just using the 2 best supported models here.

9/1/21: doubling back to creche timing part of analysis. for creching age, mixed models at step 2 have convergence issues. using lm for all the creche timing models instead. I think I already figured this out when I was working on this in Feb, but I never recorded the logic of this change, so now I had to re-determine this. 

8/27/21: there are 13 chicks in the new inp that were excluded from the thesis analysis. These chicks were excluded because they died before being measured at least 2 times, so I couldn't calculate growth rates. However, growth rates are not considered in any candidate models for estimating overall survival (2.1 above) 
