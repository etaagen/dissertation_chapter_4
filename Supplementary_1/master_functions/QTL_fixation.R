# QTL_fixation
# Function for tracking QTL fixation, called upon by ASR simulation wrapper functions

QTL_fixation <- function(trait, pop, simParam, pop_name){
  my_QTL_effects <- trait@addEff
  which_allele_is_better <- NULL
  for(i in 1:length(my_QTL_effects)){
    if(my_QTL_effects[i] < 0){ # if QTL effect is negative, favors minor
      which_allele_is_better <- c(which_allele_is_better, 0)
    } else if(my_QTL_effects[i] > 0){ # if QTL effect is positive, favors major
      which_allele_is_better <- c(which_allele_is_better, 2)
    }
  }
  which_allele_is_worse <- NULL
  for(i in 1:length(my_QTL_effects)){
    if(my_QTL_effects[i] < 0){ # if QTL effect is negative, favors minor
      which_allele_is_worse <- c(which_allele_is_worse, 2)
    } else if(my_QTL_effects[i] > 0){ # if QTL effect is positive, favors major
      which_allele_is_worse <- c(which_allele_is_worse, 0)
    }
  }
  
  # QTL fixation ratio - want a single value, how many columns have only 1 genotype
  # each column sum will be 0 or 2*nrow if fixed
  QTL_geno <- pullQtlGeno(pop = pop, simParam = simParam) %>% as.data.frame()
  fixed_QTL <- colSums(QTL_geno) %>% as.data.frame() 
  fixed_QTL <- cbind(fixed_QTL, which_allele_is_better*nrow(QTL_geno), which_allele_is_worse*nrow(QTL_geno))
  
  better_QTL_fix_ratio <- 0
  worse_QTL_fix_ratio <- 0
  for(i in 1:ncol(QTL_geno)){
    if(fixed_QTL[i, 1] == fixed_QTL[i, 2]){
      better_QTL_fix_ratio <- better_QTL_fix_ratio + 1
    } else if(fixed_QTL[i, 1] == fixed_QTL[i, 3]){
      worse_QTL_fix_ratio <- worse_QTL_fix_ratio + 1
    }
  }
  better_QTL_fix_ratio <- better_QTL_fix_ratio / ncol(QTL_geno)
  worse_QTL_fix_ratio <- worse_QTL_fix_ratio / ncol(QTL_geno)
  
  QTL_fixation <- as.data.frame(cbind(better_QTL_fix_ratio, worse_QTL_fix_ratio, pop_name))
  QTL_fixation$better_QTL_fix_ratio <- as.numeric(QTL_fixation$better_QTL_fix_ratio)
  QTL_fixation$worse_QTL_fix_ratio <- as.numeric(QTL_fixation$worse_QTL_fix_ratio)
  
  assign(paste0(pop_name, "_QTL"), QTL_fixation, envir = .GlobalEnv)
}

