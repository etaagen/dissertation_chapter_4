# my_QTL_AF_change

# track AF change of QTL in ASR simulation wrapper functions

my_QTL_AF_change <- function(trait, pop_1, pop_2, simParam){ #burnin and gsPop gives entry for GS and end of GS
  # pull QTL genotypes 
  start_QTL_geno <- pullQtlGeno(pop = pop_1, simParam = simParam) %>% as.data.frame()
  end_QTL_geno <- pullQtlGeno(pop = pop_2, simParam = simParam) %>% as.data.frame()
  QTL_colnames <- colnames(start_QTL_geno)
  # pop size
  burnin_size <- pop_1@nInd
  end_size <- pop_2@nInd
  
  # which QTL allele is + / - 
  # need to consider based on the QTL_effect_ratio if the 2 (major) or 0 (minor) genotype is -
  # if the effect is +, the 0 allele is - variant
  # if the effect is -, the 2 allele is - variant 
  my_effects <- trait@addEff
  which_allele_is_DV <- NULL
  for(i in 1:length(my_effects)){
    if(my_effects[i] < 0){
      which_allele_is_DV <- c(which_allele_is_DV, 2)
    } else if(my_effects[i] > 0){
      which_allele_is_DV <- c(which_allele_is_DV, 0)
    }
  }
  # which_allele_is_DV, now we know which allele is -
  # now let's find out if small, medium or large effect
  
  # we're going on effect size, so taking absolute value of QTL effect makes this easier
  my_QTL_effects <- abs(trait@addEff) 
  my_QTL_effects_summary <-data.frame(unclass(summary(my_QTL_effects)))
  
  QTL_effect_sizes <- NULL
  QTL_index <- NULL
  for(i in 1:length(my_QTL_effects)){
    if(my_QTL_effects[i] <= my_QTL_effects_summary[2,]){ # 1st qu.
      QTL_effect_sizes <- c(QTL_effect_sizes, "small")
      QTL_index <- c(QTL_index, QTL_colnames[i])
    }else if (my_QTL_effects[i] > my_QTL_effects_summary[2,] & my_QTL_effects[i] < my_QTL_effects_summary[5,]){ # 3rd qu.
      QTL_effect_sizes <- c(QTL_effect_sizes, "medium")
      QTL_index <- c(QTL_index, QTL_colnames[i])
    }else if (my_QTL_effects[i] >= my_QTL_effects_summary[5,]){ # 3rd qu.
      QTL_effect_sizes <- c(QTL_effect_sizes, "large")
      QTL_index <- c(QTL_index, QTL_colnames[i])
    }
  }
  
  # get start and end AF
  SNP_geno <- c(0, 2)
  burnin_AF <- sapply(start_QTL_geno, function(x) table(factor(x, levels = SNP_geno, ordered = TRUE))) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame()
  
  end_AF <- sapply(end_QTL_geno, function(x) table(factor(x, levels = SNP_geno, ordered = TRUE))) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame()

  
  test_start <- cbind(QTL_index, QTL_effect_sizes, which_allele_is_DV, burnin_AF/burnin_size)
  # define DV AF column
  for(i in 1:nrow(test_start)){
    if(test_start$which_allele_is_DV[i] == 0){
      test_start$AF_neg[i] <- test_start$`0`[i]
      test_start$AF_pos[i] <- test_start$`2`[i]
    }else{
      test_start$AF_neg[i] <- test_start$`2`[i]
      test_start$AF_pos[i] <- test_start$`0`[i]
    }
  }
 
  test_end <- cbind(QTL_index, QTL_effect_sizes, which_allele_is_DV, end_AF/end_size)
  for(i in 1:nrow(test_end)){
    if(test_end$which_allele_is_DV[i] == 0){
      test_end$AF_neg[i] <- test_end$`0`[i]
      test_end$AF_pos[i] <- test_end$`2`[i]
    }else{
      test_end$AF_neg[i] <- test_end$`2`[i]
      test_end$AF_pos[i] <- test_end$`0`[i]
    }
  }
  # how did the AF change for each category
  # large / -, large / +
  # medium / -, medium / +
  # small / -, small / - 
  
  # filter for 6 categories AF change (end - start)
  # positive values means increased freq
  # negative values means decreased freq 
  large_neg <- (test_end %>% filter(QTL_effect_sizes == "large") %>% select(AF_neg)) - (test_start %>% filter(QTL_effect_sizes == "large") %>% select(AF_neg))
  large_pos <- (test_end %>% filter(QTL_effect_sizes == "large") %>% select(AF_pos)) - (test_start %>% filter(QTL_effect_sizes == "large") %>% select(AF_pos))
  med_neg <- (test_end %>% filter(QTL_effect_sizes == "medium") %>% select(AF_neg)) - (test_start %>% filter(QTL_effect_sizes == "medium") %>% select(AF_neg))
  med_pos <- (test_end %>% filter(QTL_effect_sizes == "medium") %>% select(AF_pos)) - (test_start %>% filter(QTL_effect_sizes == "medium") %>% select(AF_pos))
  small_neg <- (test_end %>% filter(QTL_effect_sizes == "small") %>% select(AF_neg)) - (test_start %>% filter(QTL_effect_sizes == "small") %>% select(AF_neg))
  small_pos <- (test_end %>% filter(QTL_effect_sizes == "small") %>% select(AF_pos)) - (test_start %>% filter(QTL_effect_sizes == "small") %>% select(AF_pos))
  
  QTL_AF_change <- rbind(mean(large_neg$AF_neg), mean(large_pos$AF_pos),
        mean(med_neg$AF_neg),mean(med_pos$AF_pos),
        mean(small_neg$AF_neg),mean(small_pos$AF_pos)) %>% as.data.frame()
  Type <- c("large_neg", "large_pos", "med_neg", "med_pos", "small_neg", "small_pos")
  QTL_AF_change <- cbind(Type, QTL_AF_change)
  colnames(QTL_AF_change) <- c("mean_AF_change", "QTL_type")
  QTL_AF_change <<- QTL_AF_change
}
