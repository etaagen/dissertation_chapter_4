### 3 functions
# my_ASR_csv_bp and my_ASR_csv summarize output from parallel rep of ASR wrapper functions
# my_ASR_summary used to plot results 

# for a given my_ASR_output list, create df to summarize 
my_ASR_csv_bp <- function(my_ASR_output, suffix, SNP){
  # genetic gain is the first element of each list of lists 
  genetic_gain <- data.frame()
  for(j in 1:length(my_ASR_output)){
    gg <- as.data.frame((my_ASR_output[[j]][[1]]))
    genetic_gain <- rbind(genetic_gain, gg)
  }
  names(genetic_gain)[4] <- "value"
  genetic_gain$measurement <- "genetic gain"
  genetic_gain$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("gg_", suffix), genetic_gain, envir = .GlobalEnv)
  
  # genetic variance is the second element of each list of lists
  genetic_variance <- data.frame()
  for(j in 1:length(my_ASR_output)){
    gv <- as.data.frame((my_ASR_output[[j]][[2]]))
    genetic_variance <- rbind(genetic_variance, gv)
  }
  # want to save this for ANOVA 
  names(genetic_variance)[3] <- "value"
  genetic_variance$measurement <- "genetic variance"
  genetic_variance$simulation <- suffix
  assign(paste0("gv_", suffix), genetic_variance, envir = .GlobalEnv)
  
  # prediction accuracy is the third element of each list of lists
  prediction_accuracy <- data.frame()
  for(j in 1:length(my_ASR_output)){
    pa <- as.data.frame((my_ASR_output[[j]][[3]]))
    prediction_accuracy <- rbind(prediction_accuracy, pa)
  }
  # want to save this for ANOVA 
  names(prediction_accuracy)[3] <- "value"
  prediction_accuracy$measurement <- "prediction accuracy"
  prediction_accuracy$simulation <- suffix
  assign(paste0("pa_", suffix), prediction_accuracy, envir = .GlobalEnv)
  
  # bulmer effect is the fourth element of each list of lists
  bulmer_effect <- data.frame()
  for(j in 1:length(my_ASR_output)){
    be <- as.data.frame((my_ASR_output[[j]][[4]]))
    bulmer_effect <- rbind(bulmer_effect, be)
  }
  names(bulmer_effect)[3] <- "value"
  bulmer_effect$measurement <- "bulmer effect"
  bulmer_effect$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("be_", suffix), bulmer_effect, envir = .GlobalEnv)
  
  
  # if(SNP == TRUE){
  #   # LD of QTL and SNP is the fifth element of each list of lists
  #   LD_summary <- data.frame()
  #   for(j in 1:length(my_ASR_output)){
  #     LD <- as.data.frame((my_ASR_output[[j]][[5]]))
  #     LD_summary <- rbind(LD_summary, LD)
  #   }
  #   names(LD_summary)[3] <- "value"
  #   LD_summary$measurement <- "LD"
  #   LD_summary$simulation <- suffix
  #   # want to save this for ANOVA 
  #   assign(paste0("ld_", suffix), LD_summary, envir = .GlobalEnv)
  #   
  #   LD_summary_2 <- LD_summary %>% 
  #     group_by(legend, cycle) %>%
  #     summarise(LD = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
  #   LD_summary_2$measurement <- "LD"
  #   LD_summary_2$simulation <- suffix
  #   names(LD_summary_2)[3] <- "value"
  #   # want to save for plots
  #   assign(paste0("ld_summary_", suffix), LD_summary_2, envir = .GlobalEnv)
  # }
  
  # Fixed QTL ratio
  if(SNP == TRUE){ # because changes [[5 or 6 depending]] # NOW BOTH 5 because excluding LD 
    fixed_QTL <- data.frame()
    for(j in 1:length(my_ASR_output)){
      qtl <- as.data.frame((my_ASR_output[[j]][[5]]))
      fixed_QTL <- rbind(fixed_QTL, qtl)
    }
  }else{
    fixed_QTL <- data.frame()
    for(j in 1:length(my_ASR_output)){
      qtl <- as.data.frame((my_ASR_output[[j]][[5]]))
      fixed_QTL <- rbind(fixed_QTL, qtl)
    }
  }
  # need as.int.
  better <- fixed_QTL %>% select(better_QTL_fix_ratio) 
  colnames(better) <- "value"
  better$measurement <- "positive_QTL"
  worse <- fixed_QTL %>% select(worse_QTL_fix_ratio)
  colnames(worse) <- "value"
  worse$measurement <- "negative_QTL"
  
  fixed_QTL <- fixed_QTL %>% select(`rep`, cycle, legend)
  fixed_QTL <- rbind(cbind(fixed_QTL, better), cbind(fixed_QTL, worse))
  fixed_QTL$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("qtl_", suffix), fixed_QTL, envir = .GlobalEnv)
  
  
  
  # QTL AF is the 6th element of each list of lists
  my_qtl_af <- data.frame()
  for(j in 1:length(my_ASR_output)){
    my_qtl <- as.data.frame((my_ASR_output[[j]][[6]]))
    my_qtl_af <- rbind(my_qtl_af, my_qtl)
  }
  names(my_qtl_af)[2] <- "value"
  names(my_qtl_af)[1] <- "type"
  my_qtl_af$measurement <- "mean_QTL_AF_change"
  my_qtl_af$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("qtl_af_", suffix), my_qtl_af, envir = .GlobalEnv)
}

my_ASR_csv <- function(my_ASR_output, suffix, SNP){
  # genetic gain is the first element of each list of lists 
  genetic_gain <- data.frame()
  for(j in 1:length(my_ASR_output)){
    gg <- as.data.frame((my_ASR_output[[j]][[1]]))
    genetic_gain <- rbind(genetic_gain, gg)
  }
  names(genetic_gain)[4] <- "value"
  genetic_gain$measurement <- "genetic gain"
  genetic_gain$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("gg_", suffix), genetic_gain, envir = .GlobalEnv)
  
  # genetic variance is the second element of each list of lists
  genetic_variance <- data.frame()
  for(j in 1:length(my_ASR_output)){
    gv <- as.data.frame((my_ASR_output[[j]][[2]]))
    genetic_variance <- rbind(genetic_variance, gv)
  }
  # want to save this for ANOVA 
  names(genetic_variance)[3] <- "value"
  genetic_variance$measurement <- "genetic variance"
  genetic_variance$simulation <- suffix
  assign(paste0("gv_", suffix), genetic_variance, envir = .GlobalEnv)

  # prediction accuracy is the third element of each list of lists
  prediction_accuracy <- data.frame()
  for(j in 1:length(my_ASR_output)){
    pa <- as.data.frame((my_ASR_output[[j]][[3]]))
    prediction_accuracy <- rbind(prediction_accuracy, pa)
  }
  # want to save this for ANOVA 
  names(prediction_accuracy)[3] <- "value"
  prediction_accuracy$measurement <- "prediction accuracy"
  prediction_accuracy$simulation <- suffix
  assign(paste0("pa_", suffix), prediction_accuracy, envir = .GlobalEnv)
  
  # bulmer effect is the fourth element of each list of lists
  bulmer_effect <- data.frame()
  for(j in 1:length(my_ASR_output)){
    be <- as.data.frame((my_ASR_output[[j]][[4]]))
    bulmer_effect <- rbind(bulmer_effect, be)
  }
  names(bulmer_effect)[3] <- "value"
  bulmer_effect$measurement <- "bulmer effect"
  bulmer_effect$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("be_", suffix), bulmer_effect, envir = .GlobalEnv)
  
  
  # if(SNP == TRUE){
  #   # LD of QTL and SNP is the fifth element of each list of lists
  #   LD_summary <- data.frame()
  #   for(j in 1:length(my_ASR_output)){
  #     LD <- as.data.frame((my_ASR_output[[j]][[5]]))
  #     LD_summary <- rbind(LD_summary, LD)
  #   }
  #   names(LD_summary)[3] <- "value"
  #   LD_summary$measurement <- "LD"
  #   LD_summary$simulation <- suffix
  #   # want to save this for ANOVA 
  #   assign(paste0("ld_", suffix), LD_summary, envir = .GlobalEnv)
  #   
  #   LD_summary_2 <- LD_summary %>% 
  #     group_by(legend, cycle) %>%
  #     summarise(LD = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
  #   LD_summary_2$measurement <- "LD"
  #   LD_summary_2$simulation <- suffix
  #   names(LD_summary_2)[3] <- "value"
  #   # want to save for plots
  #   assign(paste0("ld_summary_", suffix), LD_summary_2, envir = .GlobalEnv)
  # }
  
  # Fixed QTL ratio
  if(SNP == TRUE){ # because changes [[5 or 6 depending]] # NOW BOTH 5 because excluding LD 
    fixed_QTL <- data.frame()
    for(j in 1:length(my_ASR_output)){
      qtl <- as.data.frame((my_ASR_output[[j]][[5]]))
      fixed_QTL <- rbind(fixed_QTL, qtl)
    }
  }else{
    fixed_QTL <- data.frame()
    for(j in 1:length(my_ASR_output)){
      qtl <- as.data.frame((my_ASR_output[[j]][[5]]))
      fixed_QTL <- rbind(fixed_QTL, qtl)
    }
  }
  # need as.int.
  better <- fixed_QTL %>% select(better_QTL_fix_ratio) 
  colnames(better) <- "value"
  better$measurement <- "positive_QTL"
  worse <- fixed_QTL %>% select(worse_QTL_fix_ratio)
  colnames(worse) <- "value"
  worse$measurement <- "negative_QTL"
  
  fixed_QTL <- fixed_QTL %>% select(`rep`, cycle, legend)
  fixed_QTL <- rbind(cbind(fixed_QTL, better), cbind(fixed_QTL, worse))
  fixed_QTL$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("qtl_", suffix), fixed_QTL, envir = .GlobalEnv)
  
  # QTL AF is the 6th element of each list of lists
  my_qtl_af <- data.frame()
  for(j in 1:length(my_ASR_output)){
    my_qtl <- as.data.frame((my_ASR_output[[j]][[6]]))
    my_qtl_af <- rbind(my_qtl_af, my_qtl)
  }
  names(my_qtl_af)[2] <- "value"
  names(my_qtl_af)[1] <- "type"
  my_qtl_af$measurement <- "mean_QTL_AF_change"
  my_qtl_af$simulation <- suffix
  # want to save this for ANOVA 
  assign(paste0("qtl_af_", suffix), my_qtl_af, envir = .GlobalEnv)

}


my_ASR_summary <- function(measurement_df, measurement, df_name){
  
  if(df_name == "qtl"){
    measurement_df_summary <- measurement_df %>% 
      group_by(legend, cycle, Matrix, measurement, Recombination) %>%
      summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
    measurement_df_summary <- measurement_df_summary %>% relocate(measurement, .after = se)
    names(measurement_df_summary)[5] <- "value"
    assign(paste0(df_name, "_summary"), measurement_df_summary, envir = .GlobalEnv)
  } 
  else if(df_name == "qtl_af"){
    measurement_df_summary <- measurement_df %>% 
      group_by(legend, type, Matrix, Recombination) %>%
      summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
    measurement_df_summary$measurement <- measurement
    names(measurement_df_summary)[5] <- "value"
    assign(paste0(df_name, "_summary"), measurement_df_summary, envir = .GlobalEnv)
  } 
  else if(df_name == "dv"){
    measurement_df_summary <- measurement_df %>% 
      group_by(legend, rep, Recombination) %>%
      summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
    measurement_df_summary$measurement <- measurement
    names(measurement_df_summary)[4] <- "value"
    # want to save for plots
    assign(paste0(df_name, "_summary"), measurement_df_summary, envir = .GlobalEnv)
  }
  else {measurement_df_summary <- measurement_df %>% 
    group_by(legend, cycle, Matrix, Recombination) %>%
    summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value))) 
  measurement_df_summary$measurement <- measurement
  names(measurement_df_summary)[5] <- "value"
  
  assign(paste0(df_name, "_summary"), measurement_df_summary, envir = .GlobalEnv)
  }
}
