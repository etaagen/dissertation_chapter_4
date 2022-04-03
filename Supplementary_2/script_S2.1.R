##### how does changing the genetic map change the recombination frequency in our bp data? 
# scroll to bottom for results 


# load packages
library(AlphaSimR)
library(tidyverse)
library(data.table)

######## load founder biparental data and functions ######## 
# measuring the F1 -> DH CO is simplest with biparental pop

# download bp data: https://github.com/etaagen/dissertation_chapter_4/blob/main/Supplementary_1/file_S1.2.csv.zip
# set file path to your download directory
founder_bp <- fread("/your_file_path")

# load custom function
source("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_1/master_functions/ASR_input_bp.R") 

QTL_per_chr = 200
H2 = 0.8

ASR_input_max_SNP_bp(founder_data = founder_bp, QTL_n = QTL_per_chr) 

genMap = genMap_input_R_WT #genMap_input_R_Chr
haplotypes = haplotype_input_R_WT #haplotype_input_R_Chr

# create random QTL effects
R = newMapPop(genMap = genMap, haplotypes = haplotypes, inbred = TRUE, ploidy = 2L)
# set the simulation parameters
SP_R = SimParam$new(R)
# get QTL effects
SP_R$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
R_QTL_effects = SP_R$traits[[1]]@addEff

map_scale = 20 # 2, 20 (1 is WT)
info_ASR_input = info_ASR_input_R_WT
QTL_effects = R_QTL_effects
QTL_effect_ratio = c(0.5, 0.5)
heritability = H2

founderPop = newMapPop(genMap = genMap, haplotypes = haplotypes, inbred = TRUE, ploidy = 2L)
genMap_change = founderPop@genMap
# Multiply genetic length by map_scale
my_genMap_change = lapply(genMap_change, function(x) map_scale*x)
# Replace genetic map
founderPop@genMap = my_genMap_change

# check matches expectations
founderPop@genMap[[1]] %>% tail()

# set the simulation parameters
SP = SimParam$new(founderPop)

########## manually add trait (QTL) and SnpChip (genotype) positions based on ASR_input() results  

# set up additive trait QTL (yield)
# assign QTL_n loci for QTl per chr
SP$addTraitA(nQtlPerChr = QTL_per_chr$Freq, mean = 0, var = 1)
# Extract the trait
trait = SP$traits[[1]]

# loci we know we want to assign to QTL 
QTL_loci <- info_ASR_input %>% filter(Type == "QTL") %>% select(marker_pos)
# change loci
trait@lociLoc = QTL_loci$marker_pos
# change the QTL_effects to be constant among studies
trait@addEff = QTL_effects 
# change the additive effects signs based on QTL_effect_ratio
# "alternate" means 1/-1/1/-1 etc.
if(QTL_effect_ratio[1] == "alternate"){
  # make everything absolute value first
  trait@addEff = abs(trait@addEff)
  change_effect = length(trait@addEff)
  # new_effect_ratio alternates +/- for each QTL
  new_effect_ratio <- c(rep(x = c(1, -1), length.out = change_effect)) 
  # length(new_effect_ratio) may not equal change_effect
  if(length(new_effect_ratio) == change_effect){
    trait@addEff = trait@addEff*new_effect_ratio
  } 
  # or instead of "alternate", could be a ratio that's not fully random, and has some ratio of selecting against minor allele
} else if(QTL_effect_ratio[1] != "alternate"){ # if not fully random
  # QTL_effect_ratio[1] is ratio that's +
  # QTL_effect_ratio[2] is ratio that -
  # make everything absolute value first
  trait@addEff = abs(trait@addEff)
  change_effect = length(trait@addEff)
  # new_effect_ratio takes half of QTL_effect_ratio[1] as + and -, and all of QTL_effect_ratio[2] is +
  new_effect_ratio <- c(rep.int(1, times = round((QTL_effect_ratio[1])*change_effect)), 
                        rep.int(-1, times = round(QTL_effect_ratio[2]*change_effect)))
  new_effect = sample(new_effect_ratio, size = length(new_effect_ratio), replace = FALSE)
  if(length(new_effect) == change_effect){
    trait@addEff = trait@addEff*new_effect
  } else if (length(new_effect) - change_effect < 0) { # if new_effect too small, add remaining QTL_effects
    difference = length(new_effect) - change_effect
    new_effect = c(new_effect, 
                   sample(new_effect_ratio, size = abs(difference), replace = TRUE))
    trait@addEff = trait@addEff*new_effect
  } else if (length(new_effect) - change_effect > 0) { # if new_effect too large, remove QTL_effects
    difference = length(new_effect) - change_effect 
    new_effect <- new_effect[1:(length(new_effect)-abs(difference))]
    trait@addEff = trait@addEff*new_effect
  }
}

# Replace trait loci
SP$switchTrait(traitPos=1, lociMap=trait)
# NOTE WE COME BACK AND RE-ASSIGN
# need the variance of the founder pop to correct for different QTL effect dist. under different conditions

# newPop(founderPop) takes into account SP settings
parents = newPop(founderPop, simParam = SP) 
parents_Var <- varA(pop = parents, simParam = SP)

# corrected QTL effects to set the same variance among conditons
corrected_addEff <- trait@addEff / as.vector(sqrt(parents_Var))
# Replace QTL additive effects 
trait@addEff <- corrected_addEff
# Replace trait loci
SP$resetPed()
SP$switchTrait(traitPos=1, lociMap=trait)

# set up SnpChip (genotypes)
# assign SnpChip_n loci for SNPs per chr 
SnpChip_n <- table(info_ASR_input$Chr, info_ASR_input$Type) %>% as.data.frame() %>% filter(Var2 == "SnpChip") %>% select(Freq)
SP$addSnpChip(SnpChip_n$Freq)
# Extract the SnpChip
genotype = SP$snpChips[[1]]
# loci we know we want to assign to genotype
SnpChip_loci <- info_ASR_input %>% filter(Type == "SnpChip") %>% select(marker_pos)
# change loci
genotype@lociLoc = SnpChip_loci$marker_pos
# Replace genotype loci
SP$snpChips[[1]]@lociLoc = genotype@lociLoc

# set error variance with heritability ~ to yield
SP$setVarE(H2 = heritability)

# keep track of recombination 
# SP$setTrackRec(TRUE) # not sure of utility yet ... 

assign("SP", SP, envir = .GlobalEnv)

# ---------------------------------------------------
# cycles of phenotypic selection burnin
# ---------------------------------------------------
# newPop(founderPop) takes into account SP settings
parents = newPop(founderPop, simParam = SP) 

# specify cross to be made
crossPlan <- matrix(1:2, nrow=1, ncol=2)
F1 = makeCross(pop = parents, crossPlan = crossPlan, simParam = SP) 
#F1_AF = pullSegSiteGeno(F1, simParam = SP) %>% as.data.frame()

# DH, make nDH progeny
nDH <- 100
DH = makeDH(F1, nDH = nDH, simParam=SP)
DH_AF = pullSegSiteGeno(DH, simParam = SP) %>% as.data.frame()

# not very pretty, but it works to turn into list for each chr so we can loop over
chr_length <- NULL
for(i in 1:length(haplotypes)){
  SNP_count <- ncol(haplotypes[[i]])
  chr_length <- c(chr_length, SNP_count)
}
chr_length <- cumsum(chr_length) + c(0, rep(1, 19), 0)
# # pull chr columns for list object
chr1A <- DH_AF[, 1:chr_length[1]]
chr1B <- DH_AF[, (chr_length[1]+1):(chr_length[2]-1)]
chr1D <- DH_AF[, (chr_length[2]):(chr_length[3]-1)]
chr2A <- DH_AF[, (chr_length[3]):(chr_length[4]-1)]
chr2B <- DH_AF[, (chr_length[4]):(chr_length[5]-1)]
chr2D <- DH_AF[, (chr_length[5]):(chr_length[6]-1)]
chr3A <- DH_AF[, (chr_length[6]):(chr_length[7]-1)]
chr3B <- DH_AF[, (chr_length[7]):(chr_length[8]-1)]
chr3D <- DH_AF[, (chr_length[8]):(chr_length[9]-1)]
chr4A <- DH_AF[, (chr_length[9]):(chr_length[10]-1)]
chr4B <- DH_AF[, (chr_length[10]):(chr_length[11]-1)]
chr4D <- DH_AF[, (chr_length[11]):(chr_length[12]-1)]
chr5A <- DH_AF[, (chr_length[12]):(chr_length[13]-1)]
chr5B <- DH_AF[, (chr_length[13]):(chr_length[14]-1)]
chr5D <- DH_AF[, (chr_length[14]):(chr_length[15]-1)]
chr6A <- DH_AF[, (chr_length[15]):(chr_length[16]-1)]
chr6B <- DH_AF[, (chr_length[16]):(chr_length[17]-1)]
chr6D <- DH_AF[, (chr_length[17]):(chr_length[18]-1)]
chr7A <- DH_AF[, (chr_length[18]):(chr_length[19]-1)]
chr7B <- DH_AF[, (chr_length[19]):(chr_length[20]-1)]
chr7D <- DH_AF[, (chr_length[20]):(chr_length[21])]

DH_AF <- list(chr1A, chr1B, chr1D,
              chr2A, chr2B, chr2D,
              chr3A, chr3B, chr3D,
              chr4A, chr4B, chr4D,
              chr5A, chr5B, chr5D,
              chr6A, chr6B, chr6D,
              chr7A, chr7B, chr7D)

total_recombination = NULL
for(h in 1:length(DH_AF)){ # each chr
  my_DH_AF = t(DH_AF[[h]]) %>% as.data.frame()
  DH_recombination = NULL
  for(j in 1:ncol(my_DH_AF)){ # each DH
    for(i in 1:(nrow(my_DH_AF)-1)){ # each SNP
      if(my_DH_AF[i,j] != my_DH_AF[i+1, j]){
        my_rec <- j
        DH_recombination <- c(DH_recombination, my_rec)
      }
    }
  }
  # number of crossovers / DH on this chr
  my_crossovers <- table(DH_recombination) %>% as.data.frame()
  my_crossovers_mean <- mean(my_crossovers$Freq)
  # mean # of crossovers on each chr
  total_recombination <- c(total_recombination, my_crossovers_mean)
}
print("map scale:") 
print(map_scale)
print("mean crossovers, 21 chr:")
print(total_recombination)

# for WT map scale:
#[1] 1.556962 1.521127 1.618421 1.389831 1.537313 1.415385 2.144444 1.523810 1.600000 1.873563 1.366197
#[12] 1.018182 1.661972 1.906977 1.752809 1.373134 1.507692 1.561644 1.779221 1.703704 1.595506

# for map_scale 2:
#[1] 2.565217 2.489362 2.718750 2.101124 2.181818 2.024096 3.309278 2.526882 1.902174 3.484536 2.159091
#[12] 1.196429 2.276596 3.268041 2.597938 2.177778 2.066667 2.368421 2.831579 2.626374 2.718750

# for map_scale 20:
#[1] 14.280000 16.590000 12.130000 11.710000 15.340000  9.170000 16.950000 14.270000  4.878788 18.770000
#[11] 13.310000  2.438776 12.330000 18.140000  7.830000 11.880000 12.630000 12.650000 22.680000 19.690000
#[21] 18.900000

