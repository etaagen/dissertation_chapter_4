### 2/8/22 
# my_anova_within() compares response variable within given parameters (set H2, QTL count, SnpChip density, relationship matrix)
# my_anova_across() compares response variable across different parameters (different H2, QTL count, SnpChip density, relationship matrix)

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

# come back to this! do you needt to subset by suffix first?
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

recombination_anova <- function(measurement_df, pop, df_name){
  # pop = "bp" (biparental) or "founder" (26 cultivars)
  info <- str_split_fixed(measurement_df$legend, "_", 3) %>% as.data.frame()
  colnames(info) <- c("type", "map", "condition")
  # change name of response variable to generic so function can be used broadly
  measurement_df <- cbind(measurement_df, info)
  measurement_df$map <- ordered(measurement_df$map,
                                levels = c("WT", "Peri", "Chr"))
  measurement_df$simulation[measurement_df$simulation == 1] <- "GW"
  measurement_df$simulation[measurement_df$simulation == 2] <- "CV"
  # does map make a difference for mean?
  # map within simulation, type, condition, cycle (ie CV, R, cycle 10, condition 3 WT vs Peri vs Chr)
  measurement_df$split_map <- paste0(measurement_df$simulation, "_", measurement_df$type, "_", measurement_df$condition, "_", measurement_df$cycle)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_map)
  # anova on each list 
  anova_map <- list()
  for(i in 1:length(measurement_df_split)){
    my_anova <- aov(value ~ map, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_map <- measurement_df_split[[i]]$split_map[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_map[[i]] <- my_anova
  }
  # convert to df
  anova_map <- do.call(rbind, lapply(anova_map, as.data.frame))
  
  # across simulation (ie. CV R cycle 10, condtion 3, WT vs GW R cycle 10, condtion 3, WT )
  measurement_df$split_sim <- paste0(measurement_df$map, "_", measurement_df$type, "_", measurement_df$condition, "_", measurement_df$cycle)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_sim)
  # anova on each list 
  anova_sim <- list()
  for(i in 1:length(measurement_df_split)){
    my_anova <- aov(value ~ simulation, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_sim <- measurement_df_split[[i]]$split_sim[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_sim[[i]] <- my_anova
  }
  # convert to df
  anova_sim <- do.call(rbind, lapply(anova_sim, as.data.frame))
  
  ### IF biparental, only one type (random)
  if(pop == "founder"){
    # across type (ie.  R CV, cycle 10, condtion 3 vs DV CV, cycle 10, condtion 3)
    measurement_df$split_type <- paste0(measurement_df$map, "_", measurement_df$simulation, "_", measurement_df$condition, "_", measurement_df$cycle)
    # split on condition
    measurement_df_split <- split(measurement_df, measurement_df$split_type)
    # anova on each list 
    anova_type <- list()
    for(i in 1:length(measurement_df_split)){
      my_anova <- aov(value ~ type, data = measurement_df_split[[i]])
      my_anova <- TukeyHSD(my_anova)
      my_anova <- as.data.frame(my_anova[[1]])
      my_anova$split_type <- measurement_df_split[[i]]$split_type[1]
      my_anova <- rownames_to_column(my_anova)
      # filter for signficant p-value
      my_anova <- my_anova %>% filter(`p adj` < 0.05)
      #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
      anova_type[[i]] <- my_anova
    }
    # convert to df
    anova_type <- do.call(rbind, lapply(anova_type, as.data.frame))
  }
  
  # combine output
  # add column to define which anova condition
  if(nrow(anova_map) == 0){
    anova_map[1,] <- NA
  } 
  anova_map$anova <- "recombination"
  anova_map$model <- paste0(df_name, " ~ recombination")
  names(anova_map)[6] <- "split"
  names(anova_map)[1] <- "TukeyHSD groups"
  anova_map$split_key <- "Matrix_QTL_Cond_Cycle"
  names(anova_map)[5] <- "P-adj"
  anova_map <- anova_map %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_sim) == 0){
    anova_sim[1,] <- NA
  } 
  anova_sim$anova <- "relationship matrix"
  anova_sim$model <- paste0(df_name, " ~ relationship matrix")
  names(anova_sim)[6] <- "split"
  names(anova_sim)[1] <- "TukeyHSD groups"
  anova_sim$split_key <- "Rec_QTL_Cond_Cycle"
  names(anova_sim)[5] <- "P-adj"
  anova_sim <- anova_sim %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(pop == "founder"){
    if(nrow(anova_type) == 0){
      anova_type[1,] <- NA
    } 
    anova_type$anova <- "QTL type"
    anova_type$model <- paste0(df_name, " ~ QTL type")
    names(anova_type)[6] <- "split"
    names(anova_type)[1] <- "TukeyHSD groups"
    anova_type$split_key <- "Rec_Matrix_Cond_Cycle"
    names(anova_type)[5] <- "P-adj"
    anova_type <- anova_type %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  }
  
  # combine
  if(pop == "founder"){
    total_anova <- rbind(anova_map, anova_sim, anova_type)
  } else{
    total_anova <- rbind(anova_map, anova_sim)
  }
  # remove rownames
  rownames(total_anova)<-NULL
  
  # add info columns
  info <- str_split_fixed(df_name, "_", 2) %>% as.data.frame()
  total_anova$measurement <- info[1,1]
  total_anova <- total_anova %>% 
    relocate(measurement, .before = anova) %>% 
    relocate(`TukeyHSD groups`, .before = split_key)
  #save
  assign(paste0("anova_", df_name), total_anova, envir = .GlobalEnv)
}

recombination_QTL_AF_anova <- function(measurement_df, pop, df_name){
  # pop = "bp" (biparental) or "founder" (26 cultivars)
  colnames(measurement_df) <- c("QTL_type", "value", "legend", "rep", "measurement", "simulation")
  info <- str_split_fixed(measurement_df$legend, "_", 3) %>% as.data.frame()
  colnames(info) <- c("type", "map", "condition")
  # change name of response variable to generic so function can be used broadly
  measurement_df <- cbind(measurement_df, info)
  measurement_df$map <- ordered(measurement_df$map,
                                levels = c("WT", "Peri", "Chr"))
  measurement_df$simulation[measurement_df$simulation == 1] <- "GW"
  measurement_df$simulation[measurement_df$simulation == 2] <- "CV"
  
  # does map make a difference for mean?
  # map within simulation, type, condition, cycle (ie CV, R, cycle 10, condition 3 WT vs Peri vs Chr)
  measurement_df$split_map <- paste0(measurement_df$simulation, "_", measurement_df$type, "_", measurement_df$condition, "_", measurement_df$QTL_type)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_map)
  # anova on each list 
  anova_map <- list()
  for(i in 1:length(measurement_df_split)){
    my_anova <- aov(value ~ map, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_map <- measurement_df_split[[i]]$split_map[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_map[[i]] <- my_anova
  }
  # convert to df
  anova_map <- do.call(rbind, lapply(anova_map, as.data.frame))
  
  # across simulation (ie. CV R cycle 10, condtion 3, WT vs GW R cycle 10, condtion 3, WT )
  measurement_df$split_sim <- paste0(measurement_df$map, "_", measurement_df$type, "_", measurement_df$condition, "_", measurement_df$QTL_type)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_sim)
  # anova on each list 
  anova_sim <- list()
  for(i in 1:length(measurement_df_split)){
    my_anova <- aov(value ~ simulation, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_sim <- measurement_df_split[[i]]$split_sim[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_sim[[i]] <- my_anova
  }
  # convert to df
  anova_sim <- do.call(rbind, lapply(anova_sim, as.data.frame))
  
  ### IF biparental, only one type (random)
  if(pop == "founder"){
    # across type (ie.  R CV, cycle 10, condtion 3 vs DV CV, cycle 10, condtion 3)
    measurement_df$split_type <- paste0(measurement_df$map, "_", measurement_df$simulation, "_", measurement_df$condition, "_", measurement_df$cycle)
    # split on condition
    measurement_df_split <- split(measurement_df, measurement_df$split_type)
    # anova on each list 
    anova_type <- list()
    for(i in 1:length(measurement_df_split)){
      my_anova <- aov(value ~ type, data = measurement_df_split[[i]])
      my_anova <- TukeyHSD(my_anova)
      my_anova <- as.data.frame(my_anova[[1]])
      my_anova$split_type <- measurement_df_split[[i]]$split_type[1]
      my_anova <- rownames_to_column(my_anova)
      # filter for signficant p-value
      my_anova <- my_anova %>% filter(`p adj` < 0.05)
      #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
      anova_type[[i]] <- my_anova
    }
    # convert to df
    anova_type <- do.call(rbind, lapply(anova_type, as.data.frame))
  }
  
  # combine output
  # add column to define which anova condition
  if(nrow(anova_map) == 0){
    anova_map[1,] <- NA
  } 
  anova_map$anova <- "recombination"
  anova_map$model <- paste0(df_name, " ~ recombination")
  names(anova_map)[6] <- "split"
  names(anova_map)[1] <- "TukeyHSD groups"
  anova_map$split_key <- "Matrix_QTL_Cond_Type"
  names(anova_map)[5] <- "P-adj"
  anova_map <- anova_map %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_sim) == 0){
    anova_sim[1,] <- NA
  } 
  anova_sim$anova <- "relationship matrix"
  anova_sim$model <- paste0(df_name, " ~ relationship matrix")
  names(anova_sim)[6] <- "split"
  names(anova_sim)[1] <- "TukeyHSD groups"
  anova_sim$split_key <- "Rec_QTL_Cond_Type"
  names(anova_sim)[5] <- "P-adj"
  anova_sim <- anova_sim %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(pop == "founder"){
    if(nrow(anova_type) == 0){
      anova_type[1,] <- NA
    } 
    anova_type$anova <- "QTL type"
    anova_type$model <- paste0(df_name, " ~ QTL type")
    names(anova_type)[6] <- "split"
    names(anova_type)[1] <- "TukeyHSD groups"
    anova_type$split_key <- "Rec_Matrix_Cond_Cycle"
    names(anova_type)[5] <- "P-adj"
    anova_type <- anova_type %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  }
  
  # combine
  if(pop == "founder"){
    total_anova <- rbind(anova_map, anova_sim, anova_type)
  } else{
    total_anova <- rbind(anova_map, anova_sim)
  }
  # remove rownames
  rownames(total_anova)<-NULL
  
  # add info columns
  #info <- str_split_fixed(df_name, "_", 2) %>% as.data.frame()
  total_anova$measurement <- df_name
  total_anova <- total_anova %>% 
    relocate(measurement, .before = anova) %>% 
    relocate(`TukeyHSD groups`, .before = split_key)
  #save
  assign(paste0("anova_", df_name), total_anova, envir = .GlobalEnv)
}









my_anova_within <- function(measurement_df, sim_rep, df_name){
  # Null hypothesis: the means of the different groups are the same
  # Alternative hypothesis: At least one sample mean is not equal to the others.
  # define the different groups 
  # type (R/DV), condition(1:5), cylce(1:6), map (WT/10X)
  bf = 0.05/sim_rep
  
  info <- str_split_fixed(measurement_df$legend, "_", 3) %>% as.data.frame()
  colnames(info) <- c("type", "map", "condition")
  # change name of response variable to generic so function can be used broadly
  measurement_df <- cbind(measurement_df, info)
  measurement_df$map <- ordered(measurement_df$map,
                                levels = c("WT", "Peri", "Chr"))
  measurement_df$cycle <- ordered(measurement_df$cycle,
                                  levels = c(1, 2, 3, 4, 5, 6))
  measurement_df$condition <- ordered(measurement_df$condition,
                                      levels = c("1", "2", "3", "4", "5"))
  measurement_df$type <- ordered(measurement_df$type,
                                 levels = c("R", "DV"))
  
  # within a given type, condition, and cycle, compare map size impact on value
  measurement_df$split_map <- paste0(measurement_df$cycle, "_", measurement_df$type, "_", measurement_df$condition)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_map)
  # anova on each list 
  anova_map <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ map, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ map, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_map <- measurement_df_split[[i]]$split_map[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_map[[i]] <- my_anova
  }
  # convert to df
  anova_map <- do.call(rbind, lapply(anova_map, as.data.frame))
  
  # within in a type, cycle, and map, compare condition impact on value  
  measurement_df$split_condition <- paste0(measurement_df$cycle, "_", measurement_df$type, "_", measurement_df$map)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_condition)
  # anova on each list 
  anova_condition <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ condition, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ condition, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_condition <- measurement_df_split[[i]]$split_condition[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_condition[[i]] <- my_anova
  }
  # convert to df
  anova_condition <- do.call(rbind, lapply(anova_condition, as.data.frame))
  
  # within a given condition, map size, and cycle, compare type impact on value
  measurement_df$split_type <- paste0(measurement_df$cycle, "_", measurement_df$condition, "_", measurement_df$map)
  # split on type
  measurement_df_split <- split(measurement_df, measurement_df$split_type)
  # anova on each list 
  anova_type <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ type, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ type, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_type <- measurement_df_split[[i]]$split_type[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_type[[i]] <- my_anova
  }
  # convert to df
  anova_type <- do.call(rbind, lapply(anova_type, as.data.frame))
  
  # within a given condition, map size, and type, compare cycle impact on value
  measurement_df$split_cycle <- paste0(measurement_df$type, "_", measurement_df$condition, "_", measurement_df$map)
  # split on cycle
  measurement_df_split <- split(measurement_df, measurement_df$split_cycle)
  # anova on each list 
  anova_cycle <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ cycle, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ cycle, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_cycle <- measurement_df_split[[i]]$split_cycle[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_cycle[[i]] <- my_anova
  }
  # convert to df
  anova_cycle <- do.call(rbind, lapply(anova_cycle, as.data.frame))
  
  
  # add column to define which anova condition
  if(nrow(anova_map) == 0){
    anova_map[1,] <- NA
  } 
  anova_map$anova <- "recombination"
  anova_map$model <- paste0(df_name, " ~ recombination")
  names(anova_map)[6] <- "split"
  names(anova_map)[1] <- "TukeyHSD groups"
  anova_map$split_key <- "cycle_type_condition"
  names(anova_map)[5] <- "P-adj"
  anova_map <- anova_map %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_condition) == 0){
    anova_condition[1,] <- NA
  }
  anova_condition$anova <- "condition"
  anova_condition$model <- paste0(df_name, " ~ condition")
  names(anova_condition)[6] <- "split"
  names(anova_condition)[1] <- "TukeyHSD groups"
  anova_condition$split_key <- "cycle_type_recombination"
  names(anova_condition)[5] <- "P-adj"
  anova_condition <- anova_condition %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_type) == 0){
    anova_type[1,] <- NA
  }
  anova_type$anova <- "type"
  anova_type$model <- paste0(df_name, " ~ type")
  names(anova_type)[6] <- "split"
  names(anova_type)[1] <- "TukeyHSD groups"
  anova_type$split_key <- "cycle_condition_recombination"
  names(anova_type)[5] <- "P-adj"
  anova_type <- anova_type %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_cycle) == 0){
    anova_cycle[1,] <- NA
  }
  anova_cycle$anova <- "cycle"
  anova_cycle$model <- paste0(df_name, " ~ cycle")
  names(anova_cycle)[6] <- "split"
  names(anova_cycle)[1] <- "TukeyHSD groups"
  anova_cycle$split_key <- "type_condition_recombination"
  names(anova_cycle)[5] <- "P-adj"
  anova_cycle <- anova_cycle %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  # how to combine without NA
  total_anova <- rbind(anova_map, anova_condition, anova_type, anova_cycle)
  rownames(total_anova)<-NULL
  # remove rownames
  # add info columns
  info <- str_split_fixed(df_name, "_", 2) %>% as.data.frame()
  total_anova$measurement <- info[1,1]
  total_anova$simulation <- info[1,2]
  total_anova <- total_anova %>% 
    relocate(measurement, .before = anova) %>% 
    relocate(simulation, .before = measurement) %>% 
    relocate(`TukeyHSD groups`, .before = split_key)
  #save
  assign(paste0("anova_", df_name), total_anova, envir = .GlobalEnv)
}

#### new anova function for measurements that are not taken every cycle 
my_anova_within_QTL_AF <- function(measurement_df, sim_rep, df_name){
  #### new anova function for measurements that are not taken every cycle 
  # specifically for AF change
  # QTL allele frequency change by effect size and Deleterious variant allele frequency change  
  # Null hypothesis: the means of the different groups are the same
  # Alternative hypothesis: At least one sample mean is not equal to the others.
  # define the different groups 
  # type (R/DV), condition(1:5), cylce(1:6), map (WT/10X)
  bf = 0.05/sim_rep
  
  info <- str_split_fixed(measurement_df$legend, "_", 3) %>% as.data.frame()
  colnames(info) <- c("type", "map", "condition")
  # change name of response variable to generic so function can be used broadly
  measurement_df <- cbind(measurement_df, info)
  measurement_df$map <- ordered(measurement_df$map,
                                levels = c("WT", "Peri", "Chr"))
  measurement_df$condition <- ordered(measurement_df$condition,
                                      levels = c("1", "2", "3", "4", "5"))
  measurement_df$type <- ordered(measurement_df$type,
                                 levels = c("R", "DV"))
  measurement_df$QTL_effect_size <- ordered(measurement_df$QTL_effect_size,
                                            levels = c("small", "large"))
  # within a given type, condition, and QTL_effect_size, compare map size impact on value
  measurement_df$split_map <- paste0(measurement_df$QTL_effect_size, "_", measurement_df$type, "_", measurement_df$condition)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_map)
  # anova on each list 
  anova_map <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ map, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ map, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_map <- measurement_df_split[[i]]$split_map[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_map[[i]] <- my_anova
  }
  # convert to df
  anova_map <- do.call(rbind, lapply(anova_map, as.data.frame))
  
  # within in a type, QTL_effect_size, and map, compare condition impact on value  
  measurement_df$split_condition <- paste0(measurement_df$QTL_effect_size, "_", measurement_df$type, "_", measurement_df$map)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_condition)
  # anova on each list 
  anova_condition <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ condition, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ condition, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_condition <- measurement_df_split[[i]]$split_condition[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_condition[[i]] <- my_anova
  }
  # convert to df
  anova_condition <- do.call(rbind, lapply(anova_condition, as.data.frame))
  
  # within a given condition, map size, and QTL_effect_size, compare type impact on value
  measurement_df$split_type <- paste0(measurement_df$QTL_effect_size, "_", measurement_df$condition, "_", measurement_df$map)
  # split on type
  measurement_df_split <- split(measurement_df, measurement_df$split_type)
  # anova on each list 
  anova_type <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ type, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ type, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_type <- measurement_df_split[[i]]$split_type[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_type[[i]] <- my_anova
  }
  # convert to df
  anova_type <- do.call(rbind, lapply(anova_type, as.data.frame))
  
  # within a given condition, map size, and type , compare QTL_effect_size impact on value
  measurement_df$split_ES <- paste0(measurement_df$type, "_", measurement_df$condition, "_", measurement_df$map)
  # split on QTL_effect_size
  measurement_df_split <- split(measurement_df, measurement_df$split_ES)
  # anova on each list 
  anova_ES <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ QTL_effect_size, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ QTL_effect_size, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_ES <- measurement_df_split[[i]]$split_ES[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_ES[[i]] <- my_anova
  }
  # convert to df
  anova_ES <- do.call(rbind, lapply(anova_ES, as.data.frame))
  
  # add column to define which anova condition
  if(nrow(anova_map) == 0){
    anova_map[1,] <- NA
  } 
  anova_map$anova <- "recombination"
  anova_map$model <- paste0(df_name, " ~ recombination")
  names(anova_map)[6] <- "split"
  names(anova_map)[1] <- "TukeyHSD groups"
  anova_map$split_key <- "QTL-ES_type_condition"
  names(anova_map)[5] <- "P-adj"
  anova_map <- anova_map %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_condition) == 0){
    anova_condition[1,] <- NA
  }
  anova_condition$anova <- "condition"
  anova_condition$model <- paste0(df_name, " ~ condition")
  names(anova_condition)[6] <- "split"
  names(anova_condition)[1] <- "TukeyHSD groups"
  anova_condition$split_key <- "QTL-ES_type_recombination"
  names(anova_condition)[5] <- "P-adj"
  anova_condition <- anova_condition %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_type) == 0){
    anova_type[1,] <- NA
  }
  anova_type$anova <- "type"
  anova_type$model <- paste0(df_name, " ~ type")
  names(anova_type)[6] <- "split"
  names(anova_type)[1] <- "TukeyHSD groups"
  anova_type$split_key <- "QTL-ES_condition_recombination"
  names(anova_type)[5] <- "P-adj"
  anova_type <- anova_type %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_ES) == 0){
    anova_ES[1,] <- NA
  }
  anova_ES$anova <- "QTL ES"
  anova_ES$model <- paste0(df_name, " ~ QTL ES")
  names(anova_ES)[6] <- "split"
  names(anova_ES)[1] <- "TukeyHSD groups"
  anova_ES$split_key <- "type_condition_recombination"
  names(anova_ES)[5] <- "P-adj"
  anova_ES <- anova_ES %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  # how to combine without NA
  total_anova <- rbind(anova_map, anova_condition, anova_type, anova_ES)
  rownames(total_anova)<-NULL
  # remove rownames
  # add info columns
  info <- str_split_fixed(df_name, "_", 2) %>% as.data.frame()
  total_anova$measurement <- info[1,1]
  total_anova$simulation <- info[1,2]
  total_anova <- total_anova %>% 
    relocate(measurement, .before = anova) %>% 
    relocate(simulation, .before = measurement) %>% 
    relocate(`TukeyHSD groups`, .before = split_key)
  #save
  assign(paste0("anova_", df_name), total_anova, envir = .GlobalEnv)
}

my_anova_within_bp <- function(measurement_df, sim_rep, df_name){
  # Null hypothesis: the means of the different groups are the same
  # Alternative hypothesis: At least one sample mean is not equal to the others.
  # define the different groups 
  # condition(1:5), cylce(1:6), map (WT/10X) (no type! remove, only R)
  bf = 0.05/sim_rep
  
  info <- str_split_fixed(measurement_df$legend, "_", 3) %>% as.data.frame()
  colnames(info) <- c("type", "map", "condition")
  # change name of response variable to generic so function can be used broadly
  measurement_df <- cbind(measurement_df, info)
  measurement_df$map <- ordered(measurement_df$map,
                                levels = c("WT", "Peri", "Chr"))
  measurement_df$cycle <- ordered(measurement_df$cycle,
                                  levels = c(1, 2, 3, 4, 5, 6))
  measurement_df$condition <- ordered(measurement_df$condition,
                                      levels = c("1", "2", "3", "4", "5"))
  # measurement_df$type <- ordered(measurement_df$type,
  #                                levels = c("R", "DV"))
  
  # within a given type, condition, and cycle, compare map size impact on value
  measurement_df$split_map <- paste0(measurement_df$cycle, "_", measurement_df$type, "_", measurement_df$condition)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_map)
  # anova on each list 
  anova_map <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ map, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ map, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_map <- measurement_df_split[[i]]$split_map[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_map[[i]] <- my_anova
  }
  # convert to df
  anova_map <- do.call(rbind, lapply(anova_map, as.data.frame))
  
  # within in a type, cycle, and map, compare condition impact on value  
  measurement_df$split_condition <- paste0(measurement_df$cycle, "_", measurement_df$type, "_", measurement_df$map)
  # split on condition
  measurement_df_split <- split(measurement_df, measurement_df$split_condition)
  # anova on each list 
  anova_condition <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ condition, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ condition, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_condition <- measurement_df_split[[i]]$split_condition[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_condition[[i]] <- my_anova
  }
  # convert to df
  anova_condition <- do.call(rbind, lapply(anova_condition, as.data.frame))
  
  # # within a given condition, map size, and cycle, compare type impact on value
  # measurement_df$split_type <- paste0(measurement_df$cycle, "_", measurement_df$condition, "_", measurement_df$map)
  # # split on type
  # measurement_df_split <- split(measurement_df, measurement_df$split_type)
  # # anova on each list 
  # anova_type <- list()
  # for(i in 1:length(measurement_df_split)){
  #   #my_anova <- anova(lm(value ~ type, measurement_df_split[[i]])) %>% as.data.frame()
  #   my_anova <- aov(value ~ type, data = measurement_df_split[[i]])
  #   my_anova <- TukeyHSD(my_anova)
  #   my_anova <- as.data.frame(my_anova[[1]])
  #   my_anova$split_type <- measurement_df_split[[i]]$split_type[1]
  #   my_anova <- rownames_to_column(my_anova)
  #   # filter for signficant p-value
  #   my_anova <- my_anova %>% filter(`p adj` < 0.05)
  #   #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
  #   anova_type[[i]] <- my_anova
  # }
  # # convert to df
  # anova_type <- do.call(rbind, lapply(anova_type, as.data.frame))
  
  # within a given condition, map size, and type, compare cycle impact on value
  measurement_df$split_cycle <- paste0(measurement_df$type, "_", measurement_df$condition, "_", measurement_df$map)
  # split on cycle
  measurement_df_split <- split(measurement_df, measurement_df$split_cycle)
  # anova on each list 
  anova_cycle <- list()
  for(i in 1:length(measurement_df_split)){
    #my_anova <- anova(lm(value ~ cycle, measurement_df_split[[i]])) %>% as.data.frame()
    my_anova <- aov(value ~ cycle, data = measurement_df_split[[i]])
    my_anova <- TukeyHSD(my_anova)
    my_anova <- as.data.frame(my_anova[[1]])
    my_anova$split_cycle <- measurement_df_split[[i]]$split_cycle[1]
    my_anova <- rownames_to_column(my_anova)
    # filter for signficant p-value
    my_anova <- my_anova %>% filter(`p adj` < 0.05)
    #my_anova <- my_anova %>% filter(`Pr(>F)` < bf)
    anova_cycle[[i]] <- my_anova
  }
  # convert to df
  anova_cycle <- do.call(rbind, lapply(anova_cycle, as.data.frame))
  
  
  # add column to define which anova condition
  if(nrow(anova_map) == 0){
    anova_map[1,] <- NA
  } 
  anova_map$anova <- "recombination"
  anova_map$model <- paste0(df_name, " ~ recombination")
  names(anova_map)[6] <- "split"
  names(anova_map)[1] <- "TukeyHSD groups"
  anova_map$split_key <- "cycle_type_condition"
  names(anova_map)[5] <- "P-adj"
  anova_map <- anova_map %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_condition) == 0){
    anova_condition[1,] <- NA
  }
  anova_condition$anova <- "condition"
  anova_condition$model <- paste0(df_name, " ~ condition")
  names(anova_condition)[6] <- "split"
  names(anova_condition)[1] <- "TukeyHSD groups"
  anova_condition$split_key <- "cycle_type_recombination"
  names(anova_condition)[5] <- "P-adj"
  anova_condition <- anova_condition %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  # 
  # if(nrow(anova_type) == 0){
  #   anova_type[1,] <- NA
  # }
  # anova_type$anova <- "type"
  # anova_type$model <- paste0(df_name, " ~ type")
  # names(anova_type)[6] <- "split"
  # names(anova_type)[1] <- "TukeyHSD groups"
  # anova_type$split_key <- "cycle_condition_recombination"
  # names(anova_type)[5] <- "P-adj"
  # anova_type <- anova_type %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  if(nrow(anova_cycle) == 0){
    anova_cycle[1,] <- NA
  }
  anova_cycle$anova <- "cycle"
  anova_cycle$model <- paste0(df_name, " ~ cycle")
  names(anova_cycle)[6] <- "split"
  names(anova_cycle)[1] <- "TukeyHSD groups"
  anova_cycle$split_key <- "type_condition_recombination"
  names(anova_cycle)[5] <- "P-adj"
  anova_cycle <- anova_cycle %>% select(`TukeyHSD groups`, anova, model, split_key, split, `P-adj`)
  
  # how to combine without NA
  total_anova <- rbind(anova_map, anova_condition, anova_cycle) # anova_type, 
  rownames(total_anova)<-NULL
  # remove rownames
  # add info columns
  info <- str_split_fixed(df_name, "_", 2) %>% as.data.frame()
  total_anova$measurement <- info[1,1]
  total_anova$simulation <- info[1,2]
  total_anova <- total_anova %>% 
    relocate(measurement, .before = anova) %>% 
    relocate(simulation, .before = measurement) %>% 
    relocate(`TukeyHSD groups`, .before = split_key)
  #save
  assign(paste0("anova_", df_name), total_anova, envir = .GlobalEnv)
}