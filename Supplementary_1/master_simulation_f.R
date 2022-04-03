######## master parallel script ######## 
# increased SnpChip for all conditions (WT, 10X, 10X_peri), mean 2,050 SNP / chr
# 2 parallel loops, first is WT vs 10X vs 10X_peri, second is CV matrix 

# set to run in BioHPC
# in BioHPC navigate to workdir
# mkdir et395/results
# copy in founder data and master_functions, and this script 
# screen -S my_code

######## master variables ######## 
replicates = 100
QTL_per_chr = 200 # 2, 200
H2 = 0.8 # 0.2, 0.8
wt_rec = 1
peri_rec = 2 # 2, 20
chr_rec = 2 # 2, 20
i = rep(QTL_per_chr, replicates)

######## load packages ######## 
library(AlphaSimR)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

######## load functions ######## 

# custom ASR_input functions to select QTL and SnpChip (max SNP or causal varaint relationship matrix approach)  
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/ASR_input_max_SNP_WT.R") 
source("/workdir/et395/master_functions/ASR_input_max_SNP_WT.R") 
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/ASR_input_Custom_DV_causal_variant.R")
source("/workdir/et395/master_functions/ASR_input_Custom_DV_causal_variant.R")

# custom ASR_wrapper function (either genomewide or causal variant relationship matrix)
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/ASR_wrapper_GSrec_genomewide.R")
source("/workdir/et395/master_functions/ASR_wrapper_GSrec_genomewide.R")
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/ASR_wrapper_GSrec_causal_variant.R")
source("/workdir/et395/master_functions/ASR_wrapper_GSrec_causal_variant.R")

# function that scales up genetic map in pericentromere, where deleterious variants are most frequent 
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/DV_scaled_map.R")
source("/workdir/et395/master_functions/DV_scaled_map.R")

# function that measures LD between QTL and SNPs - called upon by ASR_wrapper
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/QTL_SNP_LD.R")
#source("/workdir/et395/master_functions/QTL_SNP_LD.R")

# function that measures QTL fixation ratio - called upon by ASR_wrapper
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/QTL_fixation_ratio.R")
source("/workdir/et395/master_functions/QTL_fixation_ratio.R")

# function that measures QTL AF change (end - beginning) for effect size and +/- variant
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/QTL_AF_change.R")
source("/workdir/et395/master_functions/QTL_AF_change.R")

# summary stats and df prep for anova functions 
#source("/Users/ellietaagen/Desktop/github/cr_simulation/Functions/master_functions/summary_stats_anova.R")
source("/workdir/et395/master_functions/summary_stats_anova.R")

######## load data ######## 
#founder_data <- fread("/Users/ellietaagen/Dropbox/LOA 2021/snpeff impute/founder_data_impute_2.csv")
founder_data <- fread("/workdir/et395/founder_data_impute_2.csv")
founder_data <- as.data.frame(founder_data)
#str(founder_data) # genotype 1 = major allele, 0 = minor allele
#for AlphaSimR have to convert founder_data$Pos.cM. to Morgan (currently in cM)
#founder_data$Pos.cM. <- (founder_data$Pos.cM.)/100

######## set up parallel cores ######## 
n.cores <- parallel::detectCores() - 1 
# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


################################################################################################ 
######## my_ASR_output_1: SnpChip/Map WT, SnpChip/Map Peri 10X, SnpChip/Map Chr 10X, genomewide relationship matrix ######## 
my_ASR_output_1 <- foreach(
  i = rep(QTL_per_chr, replicates), #number of QTL, and number of sampling reps to run
  z = 1:length(i), #keep track of rep
  .packages = c('tidyverse', 'data.table', 'AlphaSimR', 'base')
) %dopar% {
  
  #### generate ASR input 
  
  # define the max number of QTL per chromosome
  QTL_n = i[1]  
  
  # ASR_input() takes the founder_data and selects loci for genotype (snpchip) and loci for QTL (random for traditional approach and DV ratio for DV approach), output is ASR input formatted data
  ASR_input_max_SNP(founder_data = founder_data, QTL_n = QTL_n) 
  
  # ASR_QTL_per_chr() makes sure that both info_ASR_input_DV and  info_ASR_input_R have same QTL count on each chromosome (can vary due DV "high" density on some chromosomes, see output QTL_per_chr)
  ASR_QTL_per_chr_max_SNP(info_ASR_input_DV_WT = info_ASR_input_DV_WT, info_ASR_input_R_WT = info_ASR_input_R_WT)
  
  ######### set QTL effects for DV and R   ######### 
  R = newMapPop(genMap = genMap_input_R_WT, haplotypes = haplotype_input_R_WT, inbred = TRUE, ploidy = 2L)
  DV = newMapPop(genMap = genMap_input_DV_WT, haplotypes = haplotype_input_DV_WT, inbred = TRUE, ploidy = 2L)
  # set the simulation parameters
  SP_R = SimParam$new(R)
  SP_DV = SimParam$new(DV)
  # get QTL effects
  SP_R$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  R_QTL_effects = SP_R$traits[[1]]@addEff
  SP_DV$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  DV_QTL_effects = SP_DV$traits[[1]]@addEff
  
  #### ASR wrapper 
  # R1, WT genetic map, QTL effect c(1, 0)
  rep_name = "R_WT_1"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # R2, WT genetic map, QTL effect c(0.6666, 0.3333)   
  rep_name = "R_WT_2"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # R3, WT genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_WT_3"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # R4, WT genetic map, QTL effect "alternate"  
  rep_name = "R_WT_4"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # R5, WT genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_WT_5"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # R1, 10X scaled pericentromere genetic map, QTL effect c(1, 0)
  rep_name = "R_Peri_1"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # R2, 10X scaled pericentromere genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "R_Peri_2"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # R3, 10X scaled pericentromere genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_Peri_3"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # R4, 10X scaled pericentromere genetic map, QTL effect "alternate" 
  rep_name = "R_Peri_4"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # R5, 10X scaled pericentromere genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_Peri_5"
  ASR_wrapper(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # R1, 10X scaled genetic map, QTL effect c(1, 0)
  rep_name = "R_Chr_1"
  ASR_wrapper(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # R2, 10X scaled genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "R_Chr_2"
  ASR_wrapper(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # R3, 10X scaled genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_Chr_3"
  ASR_wrapper(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # R4, 10X scaled genetic map, QTL effect "alternate" 
  rep_name = "R_Chr_4"
  ASR_wrapper(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # R5, 10X scaled genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_Chr_5"
  ASR_wrapper(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R_WT, haplotypes_input = haplotype_input_R_WT, info_ASR_input = info_ASR_input_R_WT, dv_loci = dv_loci, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  #### DV #### 
  # DV1, WT genetic map, QTL effect c(1, 0)
  rep_name = "DV_WT_1"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # DV2, WT genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_WT_2"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # DV3, WT genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_WT_3"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # DV4, WT genetic map, QTL effect "alternate"
  rep_name = "DV_WT_4"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # DV5, WT genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_WT_5"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # DV1, 10X scaled pericentromere DV genetic map, QTL effect c(1, 0) 
  rep_name = "DV_Peri_1"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # DV2, 10X scaled pericentromere DV genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_Peri_2"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # DV3, 10X scaled pericentromere DV genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_Peri_3"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # DV4, 10X scaled pericentromere DV genetic map, QTL effect "alternate"
  rep_name = "DV_Peri_4"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # DV5, 10X scaled pericentromere DV genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_Peri_5"
  ASR_wrapper(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # DV1, 10X scaled  DV genetic map, QTL effect c(1, 0) 
  rep_name = "DV_Chr_1"
  ASR_wrapper(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # DV2, 10X scaled  DV genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_Chr_2"
  ASR_wrapper(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # DV3, 10X scaled  DV genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_Chr_3"
  ASR_wrapper(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # DV4, 10X scaled  DV genetic map, QTL effect "alternate"
  rep_name = "DV_Chr_4"
  ASR_wrapper(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # DV5, 10X scaled  DV genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_Chr_5"
  ASR_wrapper(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV_WT, haplotypes_input = haplotype_input_DV_WT, info_ASR_input = info_ASR_input_DV_WT, dv_loci = dv_loci, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  
  # combine results
  my_GG_R_data = rbind(GG_R_WT_1, GG_R_WT_2, GG_R_WT_3, GG_R_WT_4, GG_R_WT_5,
                       GG_R_Peri_1, GG_R_Peri_2, GG_R_Peri_3, GG_R_Peri_4, GG_R_Peri_5,
                       GG_R_Chr_1, GG_R_Chr_2, GG_R_Chr_3, GG_R_Chr_4, GG_R_Chr_5) %>% as.data.frame()
  my_GV_R_data = rbind(GV_R_WT_1, GV_R_WT_2, GV_R_WT_3, GV_R_WT_4, GV_R_WT_5,
                       GV_R_Peri_1, GV_R_Peri_2, GV_R_Peri_3, GV_R_Peri_4, GV_R_Peri_5,
                       GV_R_Chr_1, GV_R_Chr_2, GV_R_Chr_3, GV_R_Chr_4, GV_R_Chr_5) %>% as.data.frame()
  my_PA_R_data = rbind(PA_R_WT_1, PA_R_WT_2, PA_R_WT_3, PA_R_WT_4, PA_R_WT_5,
                       PA_R_Peri_1, PA_R_Peri_2, PA_R_Peri_3, PA_R_Peri_4, PA_R_Peri_5,
                       PA_R_Chr_1, PA_R_Chr_2, PA_R_Chr_3, PA_R_Chr_4, PA_R_Chr_5) %>% as.data.frame()
  my_BE_R_data = rbind(BE_R_WT_1, BE_R_WT_2, BE_R_WT_3, BE_R_WT_4, BE_R_WT_5,
                       BE_R_Peri_1, BE_R_Peri_2, BE_R_Peri_3, BE_R_Peri_4, BE_R_Peri_5,
                       BE_R_Chr_1, BE_R_Chr_2, BE_R_Chr_3, BE_R_Chr_4, BE_R_Chr_5) %>% as.data.frame()
  # my_LD_R_data = rbind(LD_R_WT_1, LD_R_WT_2, LD_R_WT_3, LD_R_WT_4, LD_R_WT_5,
  #                      LD_R_Peri_1, LD_R_Peri_2, LD_R_Peri_3, LD_R_Peri_4, LD_R_Peri_5,
  #                      LD_R_Chr_1, LD_R_Chr_2, LD_R_Chr_3, LD_R_Chr_4, LD_R_Chr_5) %>% as.data.frame()
  my_QTL_R_data = rbind(QTL_R_WT_1, QTL_R_WT_2, QTL_R_WT_3, QTL_R_WT_4, QTL_R_WT_5,
                        QTL_R_Peri_1, QTL_R_Peri_2, QTL_R_Peri_3, QTL_R_Peri_4, QTL_R_Peri_5,
                        QTL_R_Chr_1, QTL_R_Chr_2, QTL_R_Chr_3, QTL_R_Chr_4, QTL_R_Chr_5) %>% as.data.frame()
  my_QTL_AF_R_data = rbind(QTL_AF_R_WT_1, QTL_AF_R_WT_2, QTL_AF_R_WT_3, QTL_AF_R_WT_4, QTL_AF_R_WT_5,
                        QTL_AF_R_Peri_1, QTL_AF_R_Peri_2, QTL_AF_R_Peri_3, QTL_AF_R_Peri_4, QTL_AF_R_Peri_5,
                        QTL_AF_R_Chr_1, QTL_AF_R_Chr_2, QTL_AF_R_Chr_3, QTL_AF_R_Chr_4, QTL_AF_R_Chr_5) %>% as.data.frame()
  
  my_GG_DV_data = rbind(GG_DV_WT_1, GG_DV_WT_2, GG_DV_WT_3, GG_DV_WT_4, GG_DV_WT_5,
                        GG_DV_Peri_1, GG_DV_Peri_2, GG_DV_Peri_3, GG_DV_Peri_4, GG_DV_Peri_5,
                        GG_DV_Chr_1, GG_DV_Chr_2, GG_DV_Chr_3, GG_DV_Chr_4, GG_DV_Chr_5) %>% as.data.frame()
  my_GV_DV_data = rbind(GV_DV_WT_1, GV_DV_WT_2, GV_DV_WT_3, GV_DV_WT_4, GV_DV_WT_5,
                        GV_DV_Peri_1, GV_DV_Peri_2, GV_DV_Peri_3, GV_DV_Peri_4, GV_DV_Peri_5,
                        GV_DV_Chr_1, GV_DV_Chr_2, GV_DV_Chr_3, GV_DV_Chr_4, GV_DV_Chr_5) %>% as.data.frame()
  my_PA_DV_data = rbind(PA_DV_WT_1, PA_DV_WT_2, PA_DV_WT_3, PA_DV_WT_4, PA_DV_WT_5,
                        PA_DV_Peri_1, PA_DV_Peri_2, PA_DV_Peri_3, PA_DV_Peri_4, PA_DV_Peri_5,
                        PA_DV_Chr_1, PA_DV_Chr_2, PA_DV_Chr_3, PA_DV_Chr_4, PA_DV_Chr_5) %>% as.data.frame()
  my_BE_DV_data = rbind(BE_DV_WT_1, BE_DV_WT_2, BE_DV_WT_3, BE_DV_WT_4, BE_DV_WT_5,
                        BE_DV_Peri_1, BE_DV_Peri_2, BE_DV_Peri_3, BE_DV_Peri_4, BE_DV_Peri_5,
                        BE_DV_Chr_1, BE_DV_Chr_2, BE_DV_Chr_3, BE_DV_Chr_4, BE_DV_Chr_5) %>% as.data.frame()
  # my_LD_DV_data = rbind(LD_DV_WT_1, LD_DV_WT_2, LD_DV_WT_3, LD_DV_WT_4, LD_DV_WT_5,
  #                       LD_DV_Peri_1, LD_DV_Peri_2, LD_DV_Peri_3, LD_DV_Peri_4, LD_DV_Peri_5,
  #                       LD_DV_Chr_1, LD_DV_Chr_2, LD_DV_Chr_3, LD_DV_Chr_4, LD_DV_Chr_5) %>% as.data.frame()
  my_QTL_DV_data = rbind(QTL_DV_WT_1, QTL_DV_WT_2, QTL_DV_WT_3, QTL_DV_WT_4, QTL_DV_WT_5,
                         QTL_DV_Peri_1, QTL_DV_Peri_2, QTL_DV_Peri_3, QTL_DV_Peri_4, QTL_DV_Peri_5,
                         QTL_DV_Chr_1, QTL_DV_Chr_2, QTL_DV_Chr_3, QTL_DV_Chr_4, QTL_DV_Chr_5) %>% as.data.frame()
  my_QTL_AF_DV_data = rbind(QTL_AF_DV_WT_1, QTL_AF_DV_WT_2, QTL_AF_DV_WT_3, QTL_AF_DV_WT_4, QTL_AF_DV_WT_5,
                           QTL_AF_DV_Peri_1, QTL_AF_DV_Peri_2, QTL_AF_DV_Peri_3, QTL_AF_DV_Peri_4, QTL_AF_DV_Peri_5,
                           QTL_AF_DV_Chr_1, QTL_AF_DV_Chr_2, QTL_AF_DV_Chr_3, QTL_AF_DV_Chr_4, QTL_AF_DV_Chr_5) %>% as.data.frame()
 
  # add simulation rep column
  
  my_QTL_AF_R_data$rep <- z
  my_QTL_AF_DV_data$rep <- z
  my_QTL_AF_data <- rbind(my_QTL_AF_R_data, my_QTL_AF_DV_data)
    
  my_GG_data = rbind(my_GG_R_data, my_GG_DV_data) 
  rep <- rep(z, nrow(my_GG_data)) 
  my_GG_data = cbind(rep, my_GG_data)
  
  my_GV_data = rbind(my_GV_R_data, my_GV_DV_data)
  rep <- rep(z, nrow(my_GV_data)) 
  my_GV_data = cbind(rep, my_GV_data)
  
  my_PA_data = rbind(my_PA_R_data, my_PA_DV_data)
  rep <- rep(z, nrow(my_PA_data)) 
  my_PA_data = cbind(rep, my_PA_data)
  
  my_BE_data = rbind(my_BE_R_data, my_BE_DV_data)
  rep <- rep(z, nrow(my_BE_data)) 
  my_BE_data = cbind(rep, my_BE_data)
  
  # my_LD_data = rbind(my_LD_R_data, my_LD_DV_data)
  # rep <- rep(z, nrow(my_LD_data)) 
  # my_LD_data = cbind(rep, my_LD_data)
  # 
  my_QTL_data = rbind(my_QTL_R_data, my_QTL_DV_data)
  rep <- rep(z, nrow(my_QTL_data)) 
  my_QTL_data = cbind(rep, my_QTL_data)
  
  
  
  # make list format for output
  my_output <- list(my_GG_data, my_GV_data, my_PA_data, my_BE_data, my_QTL_data, my_QTL_AF_data) #my_LD_data
}

######## my_ASR_output_2: causal variant relationship matrix ########
my_ASR_output_2 <- foreach(
  i = rep(QTL_per_chr, replicates), #number of QTL, and number of sampling reps to run
  z = 1:length(i), #keep track of rep
  .packages = c('tidyverse', 'data.table', 'AlphaSimR', 'base')
) %dopar% {

  #### generate ASR input

  # define the max number of QTL per chromosome
  QTL_n = i[1]

  # ASR_input() takes the founder_data and selects loci for genotype (snpchip) and loci for QTL (random for traditional approach and DV ratio for DV approach), output is ASR input formatted data
  ASR_input_CV(founder_data = founder_data, QTL_n = QTL_n)

  # ASR_QTL_per_chr() makes sure that both info_ASR_input_DV and  info_ASR_input_R have same QTL count on each chromosome (can vary due DV "high" density on some chromosomes, see output QTL_per_chr)
  ASR_QTL_per_chr_CV(info_ASR_input_DV = info_ASR_input_DV, info_ASR_input_R = info_ASR_input_R)

  ######### set QTL effects for DV and R   #########
  R = newMapPop(genMap = genMap_input_R, haplotypes = haplotype_input_R, inbred = TRUE, ploidy = 2L)
  DV = newMapPop(genMap = genMap_input_DV, haplotypes = haplotype_input_DV, inbred = TRUE, ploidy = 2L)
  # set the simulation parameters
  SP_R = SimParam$new(R)
  SP_DV = SimParam$new(DV)
  # get QTL effects
  SP_R$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  R_QTL_effects = SP_R$traits[[1]]@addEff
  SP_DV$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
  DV_QTL_effects = SP_DV$traits[[1]]@addEff

  #### ASR wrapper CV
  # R1, WT genetic map, QTL effect c(1, 0)
  rep_name = "R_WT_1"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R,  QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)

  # R2, WT genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "R_WT_2"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R,  QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)

  # R3, WT genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_WT_3"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R,  QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)

  # R4, WT genetic map, QTL effect "alternate"
  rep_name = "R_WT_4"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R,  QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)

  # R5, WT genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_WT_5"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R,  QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)

  # R1, 10X scaled pericentromere genetic map, QTL effect c(1, 0)
  rep_name = "R_Peri_1"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # R2, 10X scaled pericentromere genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "R_Peri_2"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # R3, 10X scaled pericentromere genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_Peri_3"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # R4, 10X scaled pericentromere genetic map, QTL effect "alternate" 
  rep_name = "R_Peri_4"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # R5, 10X scaled pericentromere genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_Peri_5"
  ASR_wrapper_CV(type = "R", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # R1, 10X scaled genetic map, QTL effect c(1, 0)
  rep_name = "R_Chr_1"
  ASR_wrapper_CV(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # R2, 10X scaled genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "R_Chr_2"
  ASR_wrapper_CV(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # R3, 10X scaled genetic map, QTL effect c(0.5, 0.5)
  rep_name = "R_Chr_3"
  ASR_wrapper_CV(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # R4, 10X scaled genetic map, QTL effect "alternate" 
  rep_name = "R_Chr_4"
  ASR_wrapper_CV(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # R5, 10X scaled genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "R_Chr_5"
  ASR_wrapper_CV(type = "R", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_R, haplotypes_input = haplotype_input_R, info_ASR_input = info_ASR_input_R, QTL_effects = R_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  #### DV ####
  # DV1, WT genetic map, QTL effect c(1, 0)
  rep_name = "DV_WT_1"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV,  QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)

  # DV2, WT genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_WT_2"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV,  QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)

  # DV3, WT genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_WT_3"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV,  QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)

  # DV4, WT genetic map, QTL effect "alternate"
  rep_name = "DV_WT_4"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV,  QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)

  # DV5, WT genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_WT_5"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV,  QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)

  # DV1, 10X scaled pericentromere DV genetic map, QTL effect c(1, 0) 
  rep_name = "DV_Peri_1"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # DV2, 10X scaled pericentromere DV genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_Peri_2"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # DV3, 10X scaled pericentromere DV genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_Peri_3"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # DV4, 10X scaled pericentromere DV genetic map, QTL effect "alternate"
  rep_name = "DV_Peri_4"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # DV5, 10X scaled pericentromere DV genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_Peri_5"
  ASR_wrapper_CV(type = "DV", map_scale = wt_rec, peri_scale = peri_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # DV1, 10X scaled  DV genetic map, QTL effect c(1, 0) 
  rep_name = "DV_Chr_1"
  ASR_wrapper_CV(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(1, 0), heritability = H2, legend = rep_name)
  
  # DV2, 10X scaled  DV genetic map, QTL effect c(0.6666, 0.3333)
  rep_name = "DV_Chr_2"
  ASR_wrapper_CV(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.6666, 0.3333), heritability = H2, legend = rep_name)
  
  # DV3, 10X scaled  DV genetic map, QTL effect c(0.5, 0.5)
  rep_name = "DV_Chr_3"
  ASR_wrapper_CV(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.5, 0.5), heritability = H2, legend = rep_name)
  
  # DV4, 10X scaled  DV genetic map, QTL effect "alternate"
  rep_name = "DV_Chr_4"
  ASR_wrapper_CV(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = "alternate", heritability = H2, legend = rep_name)
  
  # DV5, 10X scaled  DV genetic map, QTL effect c(0.3333, 0.6666)
  rep_name = "DV_Chr_5"
  ASR_wrapper_CV(type = "DV", map_scale = chr_rec, peri_scale = wt_rec, genMap_input = genMap_input_DV, haplotypes_input = haplotype_input_DV, info_ASR_input = info_ASR_input_DV, QTL_effects = DV_QTL_effects, QTL_effect_ratio = c(0.3333, 0.6666), heritability = H2, legend = rep_name)
  
  # combine results
  my_GG_R_data = rbind(GG_R_WT_1, GG_R_WT_2, GG_R_WT_3, GG_R_WT_4, GG_R_WT_5,
                       GG_R_Peri_1, GG_R_Peri_2, GG_R_Peri_3, GG_R_Peri_4, GG_R_Peri_5,
                       GG_R_Chr_1, GG_R_Chr_2, GG_R_Chr_3, GG_R_Chr_4, GG_R_Chr_5) %>% as.data.frame()
  my_GV_R_data = rbind(GV_R_WT_1, GV_R_WT_2, GV_R_WT_3, GV_R_WT_4, GV_R_WT_5,
                       GV_R_Peri_1, GV_R_Peri_2, GV_R_Peri_3, GV_R_Peri_4, GV_R_Peri_5,
                       GV_R_Chr_1, GV_R_Chr_2, GV_R_Chr_3, GV_R_Chr_4, GV_R_Chr_5) %>% as.data.frame()
  my_PA_R_data = rbind(PA_R_WT_1, PA_R_WT_2, PA_R_WT_3, PA_R_WT_4, PA_R_WT_5,
                       PA_R_Peri_1, PA_R_Peri_2, PA_R_Peri_3, PA_R_Peri_4, PA_R_Peri_5,
                       PA_R_Chr_1, PA_R_Chr_2, PA_R_Chr_3, PA_R_Chr_4, PA_R_Chr_5) %>% as.data.frame()
  my_BE_R_data = rbind(BE_R_WT_1, BE_R_WT_2, BE_R_WT_3, BE_R_WT_4, BE_R_WT_5,
                       BE_R_Peri_1, BE_R_Peri_2, BE_R_Peri_3, BE_R_Peri_4, BE_R_Peri_5,
                       BE_R_Chr_1, BE_R_Chr_2, BE_R_Chr_3, BE_R_Chr_4, BE_R_Chr_5) %>% as.data.frame()
  # my_LD_R_data = rbind(LD_R_WT_1, LD_R_WT_2, LD_R_WT_3, LD_R_WT_4, LD_R_WT_5,
  #                      LD_R_Peri_1, LD_R_Peri_2, LD_R_Peri_3, LD_R_Peri_4, LD_R_Peri_5,
  #                      LD_R_Chr_1, LD_R_Chr_2, LD_R_Chr_3, LD_R_Chr_4, LD_R_Chr_5) %>% as.data.frame()
  my_QTL_R_data = rbind(QTL_R_WT_1, QTL_R_WT_2, QTL_R_WT_3, QTL_R_WT_4, QTL_R_WT_5,
                        QTL_R_Peri_1, QTL_R_Peri_2, QTL_R_Peri_3, QTL_R_Peri_4, QTL_R_Peri_5,
                        QTL_R_Chr_1, QTL_R_Chr_2, QTL_R_Chr_3, QTL_R_Chr_4, QTL_R_Chr_5) %>% as.data.frame()
  my_QTL_AF_R_data = rbind(QTL_AF_R_WT_1, QTL_AF_R_WT_2, QTL_AF_R_WT_3, QTL_AF_R_WT_4, QTL_AF_R_WT_5,
                           QTL_AF_R_Peri_1, QTL_AF_R_Peri_2, QTL_AF_R_Peri_3, QTL_AF_R_Peri_4, QTL_AF_R_Peri_5,
                           QTL_AF_R_Chr_1, QTL_AF_R_Chr_2, QTL_AF_R_Chr_3, QTL_AF_R_Chr_4, QTL_AF_R_Chr_5) %>% as.data.frame()

  my_GG_DV_data = rbind(GG_DV_WT_1, GG_DV_WT_2, GG_DV_WT_3, GG_DV_WT_4, GG_DV_WT_5,
                        GG_DV_Peri_1, GG_DV_Peri_2, GG_DV_Peri_3, GG_DV_Peri_4, GG_DV_Peri_5,
                        GG_DV_Chr_1, GG_DV_Chr_2, GG_DV_Chr_3, GG_DV_Chr_4, GG_DV_Chr_5) %>% as.data.frame()
  my_GV_DV_data = rbind(GV_DV_WT_1, GV_DV_WT_2, GV_DV_WT_3, GV_DV_WT_4, GV_DV_WT_5,
                        GV_DV_Peri_1, GV_DV_Peri_2, GV_DV_Peri_3, GV_DV_Peri_4, GV_DV_Peri_5,
                        GV_DV_Chr_1, GV_DV_Chr_2, GV_DV_Chr_3, GV_DV_Chr_4, GV_DV_Chr_5) %>% as.data.frame()
  my_PA_DV_data = rbind(PA_DV_WT_1, PA_DV_WT_2, PA_DV_WT_3, PA_DV_WT_4, PA_DV_WT_5,
                        PA_DV_Peri_1, PA_DV_Peri_2, PA_DV_Peri_3, PA_DV_Peri_4, PA_DV_Peri_5,
                        PA_DV_Chr_1, PA_DV_Chr_2, PA_DV_Chr_3, PA_DV_Chr_4, PA_DV_Chr_5) %>% as.data.frame()
  my_BE_DV_data = rbind(BE_DV_WT_1, BE_DV_WT_2, BE_DV_WT_3, BE_DV_WT_4, BE_DV_WT_5,
                        BE_DV_Peri_1, BE_DV_Peri_2, BE_DV_Peri_3, BE_DV_Peri_4, BE_DV_Peri_5,
                        BE_DV_Chr_1, BE_DV_Chr_2, BE_DV_Chr_3, BE_DV_Chr_4, BE_DV_Chr_5) %>% as.data.frame()
  # my_LD_DV_data = rbind(LD_DV_WT_1, LD_DV_WT_2, LD_DV_WT_3, LD_DV_WT_4, LD_DV_WT_5,
  #                       LD_DV_Peri_1, LD_DV_Peri_2, LD_DV_Peri_3, LD_DV_Peri_4, LD_DV_Peri_5,
  #                       LD_DV_Chr_1, LD_DV_Chr_2, LD_DV_Chr_3, LD_DV_Chr_4, LD_DV_Chr_5) %>% as.data.frame()
  my_QTL_DV_data = rbind(QTL_DV_WT_1, QTL_DV_WT_2, QTL_DV_WT_3, QTL_DV_WT_4, QTL_DV_WT_5,
                         QTL_DV_Peri_1, QTL_DV_Peri_2, QTL_DV_Peri_3, QTL_DV_Peri_4, QTL_DV_Peri_5,
                         QTL_DV_Chr_1, QTL_DV_Chr_2, QTL_DV_Chr_3, QTL_DV_Chr_4, QTL_DV_Chr_5) %>% as.data.frame()
  my_QTL_AF_DV_data = rbind(QTL_AF_DV_WT_1, QTL_AF_DV_WT_2, QTL_AF_DV_WT_3, QTL_AF_DV_WT_4, QTL_AF_DV_WT_5,
                            QTL_AF_DV_Peri_1, QTL_AF_DV_Peri_2, QTL_AF_DV_Peri_3, QTL_AF_DV_Peri_4, QTL_AF_DV_Peri_5,
                            QTL_AF_DV_Chr_1, QTL_AF_DV_Chr_2, QTL_AF_DV_Chr_3, QTL_AF_DV_Chr_4, QTL_AF_DV_Chr_5) %>% as.data.frame()
   
  # # add simulation rep column

  my_QTL_AF_R_data$rep <- z
  my_QTL_AF_DV_data$rep <- z
  my_QTL_AF_data <- rbind(my_QTL_AF_R_data, my_QTL_AF_DV_data)
  
  my_GG_data = rbind(my_GG_R_data, my_GG_DV_data)
  rep <- rep(z, nrow(my_GG_data))
  my_GG_data = cbind(rep, my_GG_data)

  my_GV_data = rbind(my_GV_R_data, my_GV_DV_data)
  rep <- rep(z, nrow(my_GV_data))
  my_GV_data = cbind(rep, my_GV_data)

  my_PA_data = rbind(my_PA_R_data, my_PA_DV_data)
  rep <- rep(z, nrow(my_PA_data))
  my_PA_data = cbind(rep, my_PA_data)

  my_BE_data = rbind(my_BE_R_data, my_BE_DV_data)
  rep <- rep(z, nrow(my_BE_data))
  my_BE_data = cbind(rep, my_BE_data)

  # my_LD_data = rbind(my_LD_R_data, my_LD_DV_data)
  # rep <- rep(z, nrow(my_LD_data))
  # my_LD_data = cbind(rep, my_LD_data)

  my_QTL_data = rbind(my_QTL_R_data, my_QTL_DV_data)
  rep <- rep(z, nrow(my_QTL_data))
  my_QTL_data = cbind(rep, my_QTL_data)


  # make list format for output
  my_output <- list(my_GG_data, my_GV_data, my_PA_data, my_BE_data, my_QTL_data, my_QTL_AF_data)  #my_LD_data,, my_DVAF_data
}

######## summarize the my_ASR_output lists ######## 

my_ASR_csv(my_ASR_output = my_ASR_output_1, suffix = 1, SNP = TRUE)
my_ASR_csv(my_ASR_output = my_ASR_output_2, suffix = 2, SNP = FALSE)

# for ANOVAs
gg <- rbind(gg_1, gg_2)
gv <- rbind(gv_1, gv_2)
be <- rbind(be_1, be_2)
pa <- rbind(pa_1, pa_2)
#ld <- ld_1
qtl <- rbind(qtl_1, qtl_2)
qtl_af <- rbind(qtl_af_1, qtl_af_2)

fwrite(gg, "/workdir/et395/results/gg.csv")
fwrite(gv, "/workdir/et395/results/gv.csv")
fwrite(be, "/workdir/et395/results/be.csv")
fwrite(pa, "/workdir/et395/results/pa.csv")
#fwrite(ld, "/workdir/et395/results/ld.csv")
fwrite(qtl, "/workdir/et395/results/qtl.csv")
fwrite(qtl_af, "/workdir/et395/results/qtl_af.csv")

# stop cluster
parallel::stopCluster(cl = my.cluster)
