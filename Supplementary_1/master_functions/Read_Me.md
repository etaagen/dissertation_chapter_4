### Master function descriptions   

All the master functions are called within [master_simulation_f.R](https://github.com/etaagen/dissertation_chapter_4/blob/main/Supplementary_1/master_simulation_f.R) or [master_simulation_bp.R](https://github.com/etaagen/dissertation_chapter_4/blob/main/Supplementary_1/master_simulation_bp.R). To reproduce the simulations we recommend downloading the master scripts and running on a server's Terminal in parallel (the master functions do not need to be downloaded, accessible via GitHub raw HTML links).  

---- 

**ASR input functions** create ASR wrapper formatted input variables from raw data    

Function arguments: 
`founder_data` (file_S1.1.csv or file_S1.2.csv)   
`QTL_n` (number of QTL per chromosome, i.e., 2 or 200)     

*The BP and CV functions were adapted from F and DV, so lots of arguments are just commented out, i.e., conditions 1, 2 and 5, the DV QTL, or the SnpChip*   

* ASR_input_bp_cv.R, for biparental raw data and CV matrix, `ASR_input_CV_bp()`  

* ASR_input_bp_gw.R, for biparental raw data and GW matrix, `ASR_input_max_SNP_bp()`  

* ASR_input_f_cv.R,  for full founder raw data and CV matrix, `ASR_input_CV(), ASR_QTL_per_chr_CV()`  

* ASR_input_f_gw.R,  for full founder raw data and CV matrix, `ASR_input_max_SNP(),  ASR_QTL_per_chr_max_SNP()`   

---- 

**ASR wrapper functions** for given parameters, simulate genomic selection breeding program with increased recombination   

Function arguments:   
`type` ("R" (QTL) or "DV" (QTL))    
`map_scale` (numeric, i.e., 1 for WT, 20 for 20X recombination chromosome wide)   
`peri_scale` (numeric, i.e., 2 for 2X recombination in pericentromere)   
`genMap_input` (an output of ASR input functions, list of SNP positions by chromosome)  
`haplotypes_input` (an output of ASR input functions, list of haplotypes by chromosome)  
`info_ASR_input` (an output of ASR input functions, relates to _R or _DV, for choosing specific SNP and QTL sites)   
`dv_loci` (an output of ASR input functions, optional as not called within function but useful for tracking loci)  
`QTL_effects` (assigned in parallel script, a list of QTL effects to keep constant among the different simulation parameters for a given rep)  
`QTL_effect_ratio` (formatted as "alternate" or c(X, Y) where alternate has QTL effect go +/-/+/- etc. for condtion 4, X is ratio of + QTL effect (select against minor allele), Y is ratio of - QTL effect (select against major allele), takes alues between 0 and 1, ie. c(0.5, 0.5) is half + and half - condtion 3)   
`heritability` (numeric, assigned in parallel script)  
`legend`, (DF name info column, assigned in parallel script)  

* ASR_wrapper_bp_cv.R, wrapper for simulation parameters, bp pop, CV matrix, `ASR_wrapper_CV_bp()`   

* ASR_wrapper_bp_gw.R, wrapper for simulation parameters, bp pop, GW matrix, `ASR_wrapper_bp()`     

* ASR_wrapper_f_cv.R, wrapper for simulation parameters, full pop, CV matrix, `ASR_wrapper_CV()`   

* ASR_wrapper_f_gw.R, wrapper for simulation parameters, full pop, GW matrix, `ASR_wrapper()`  

---- 

**Additional functions**  

* QTL_AF.R, `my_QTL_AF_change()` measures the QTL AF change during genomic selection, called within ASR wrapper functions. Function arguments: `trait` (created in ASR wrapper), `pop_1` (population, start AF), `pop_2` (population, end AF), `simParam` (created in ASR wrapper)  

* QTL_fixation.R, `QTL_fixation()` measures the % of + and - QTL allele fixed in population, called within ASR wrapper functions. Function arguments: `trait` (created in ASR wrapper), `pop` (population, given breeding cycle), `simParam` (created in ASR wrapper), `pop_name` (character, given breeding cycle)   

* summary_stats.R, `my_ASR_csv_bp()` and `my_ASR_csv()` combine the 100 rep parallel output in the master simulation scripts. Functioin arguments `my_ASR_output` (created in master simulation scripts), suffix (keeps track of GW vs CV parallel output), `SNP` (is there a SnpChip? TRUE for GW, FALSE for CV). `my_ASR_summary()` is called within plots_S2 markdown files, and measures the mean, SD and SE for each response variable across the simulation parameter's 100 reps. Function argumants `measurement_df` (given response variable output of parallel script, or found in /Supplementary_2/results_S2.1/, `measurement` (character, name of response variable), `df_name` (character, name of output dataframe).  
