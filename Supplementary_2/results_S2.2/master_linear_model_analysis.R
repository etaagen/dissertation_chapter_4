#### significance test methods #### 
# create mixed model with fixed interaction of all reasonable factors and random effect rep
# run Anova on all 2-way interactiosn (only 2 variables)
# if significant, investigate the relationship of the facto levels / all other significant factors using emmeans()

# to run this on pop (Full 26 vs BP, need to add pop as factor to model)

library(tidyverse)
library(data.table)
library(lme4)
library(emmeans)
library(car)

### my_data can ge gg, gv, pa, be, qtl, qtl_af
my_data <- fread("/Users/ellietaagen/Desktop/github/cr_simulation/Results/BP_F_master_results_csv/f_gg.csv") %>% as.data.frame()

# note qtl and qtl_af have additional factors!
response = "Genetic gain"
df_name = "f_gg"
my_cycle = 10 # 6, 10, "total"
QTL = FALSE
QTL_AF = FALSE
Full_founder = TRUE
# filter for cycle 10 (might filter for cycle 6 later)
if(QTL_AF == FALSE){
  my_data <- my_data %>% filter(cycle == my_cycle)
}
if(QTL_AF == TRUE){
  my_data <- my_data %>% filter(type %in% c("large_neg", "med_neg", "small_neg"))
}

# set factors
my_data$rep <- as.factor(my_data$rep)
my_data$Map_type <- as.factor(my_data$Map_type)
my_data$Recombination<- as.factor(my_data$Recombination)
my_data$QTL<- as.factor(my_data$QTL)
my_data$H2<- as.factor(my_data$H2)
my_data$Repulsion<- as.factor(my_data$Repulsion)
my_data$Matrix<- as.factor(my_data$Matrix)
my_data$QTL_type <- as.factor(my_data$QTL_type)
#my_data$cycle <- as.factor(my_data$cycle)
my_data$Pop <- as.factor(my_data$Pop)
# extras for QTL
if(QTL == TRUE){
  names(my_data)[5] <- "Allele"
  my_data$Allele <- as.factor(my_data$Allele)
}
if(QTL_AF == TRUE){
  my_data <- my_data %>% filter(type %in% c("small_neg", "med_neg", "large_neg"))
  names(my_data)[1] <- "Allele"
  my_data$Allele <- as.factor(my_data$Allele)
}


###### make mixed model
# if full founder, QTL
# else if full founder
# else if bp, QTL
# else if bp
if(Full_founder == TRUE & (QTL == TRUE | QTL_AF == TRUE)){
  my_lmer <- lmer(value ~ Map_type*Recombination*QTL*H2*Repulsion*Matrix*QTL_type*Allele + (1|rep), data = my_data)
}else if(Full_founder == TRUE & QTL == FALSE &  QTL_AF == FALSE){
  my_lmer <- lmer(value ~ Map_type*Recombination*QTL*H2*Repulsion*Matrix*QTL_type + (1|rep), data = my_data)
} else if(Full_founder == FALSE & (QTL == TRUE | QTL_AF == TRUE)){
  my_lmer <- lmer(value ~ Map_type*Recombination*QTL*H2*Repulsion*Matrix*Allele + (1|rep), data = my_data)
}else{
  my_lmer <- lmer(value ~ Map_type*Recombination*QTL*H2*Repulsion*Matrix + (1|rep), data = my_data)
}

##### decide which pairwise interactions to investigate 
my_Anova <- Anova(my_lmer) 
my_Anova[1:28,]

##### post-hoc comparision of pairwise interactions

if(Full_founder == TRUE & (QTL == TRUE | QTL_AF == TRUE)){
  my_emmeans <- emmeans(my_lmer, pairwise ~ Map_type|Recombination|QTL|H2|Repulsion|Matrix|Allele|QTL_type)
  Map_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Recombination|Map_type|QTL|H2|Repulsion|Matrix|Allele|QTL_type)
  Recombination_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Matrix|Recombination|Map_type|QTL|H2|Repulsion|Allele|QTL_type)
  Matrix_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL_type|Matrix|Recombination|Map_type|QTL|H2|Repulsion|Allele)
  QTL_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Repulsion|Matrix|Recombination|Map_type|QTL|H2|Allele|QTL_type)
  Repulsion_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL|Map_type|Recombination|H2|Repulsion|Matrix|Allele|QTL_type)
  QTL_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ H2|Recombination|Map_type|QTL|Repulsion|Matrix|Allele|QTL_type)
  H2_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Allele|H2|Recombination|Map_type|QTL|Repulsion|Matrix|QTL_type)
  Allele_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
}else if(Full_founder == TRUE & QTL == FALSE &  QTL_AF == FALSE){
  my_emmeans <- emmeans(my_lmer, pairwise ~ Map_type|Recombination|QTL|H2|Repulsion|Matrix|QTL_type)
  Map_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Recombination|Map_type|QTL|H2|Repulsion|Matrix|QTL_type)
  Recombination_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Matrix|Recombination|Map_type|QTL|H2|Repulsion|QTL_type)
  Matrix_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL_type|Matrix|Recombination|Map_type|QTL|H2|Repulsion)
  QTL_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Repulsion|Matrix|Recombination|Map_type|QTL|H2|QTL_type)
  Repulsion_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL|Map_type|Recombination|H2|Repulsion|Matrix|QTL_type)
  QTL_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ H2|Recombination|Map_type|QTL|Repulsion|Matrix|QTL_type)
  H2_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  
} else if(Full_founder == FALSE & (QTL == TRUE | QTL_AF == TRUE)){
  my_emmeans <- emmeans(my_lmer, pairwise ~ Map_type|Recombination|QTL|H2|Repulsion|Matrix|Allele)
  Map_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Recombination|Map_type|QTL|H2|Repulsion|Matrix|Allele)
  Recombination_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Matrix|Recombination|Map_type|QTL|H2|Repulsion|Allele)
  Matrix_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Repulsion|Matrix|Recombination|Map_type|QTL|H2|Allele)
  Repulsion_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL|Map_type|Recombination|H2|Repulsion|Matrix|Allele)
  QTL_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ H2|Recombination|Map_type|QTL|Repulsion|Matrix|Allele)
  H2_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Allele|H2|Recombination|Map_type|QTL|Repulsion|Matrix)
  Allele_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
}else{
  my_emmeans <- emmeans(my_lmer, pairwise ~ Map_type|Recombination|QTL|H2|Repulsion|Matrix)
  Map_type_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Recombination|Map_type|QTL|H2|Repulsion|Matrix)
  Recombination_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Matrix|Recombination|Map_type|QTL|H2|Repulsion)
  Matrix_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ Repulsion|Matrix|Recombination|Map_type|QTL|H2)
  Repulsion_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ QTL|Map_type|Recombination|H2|Repulsion|Matrix)
  QTL_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05) 
  
  my_emmeans <- emmeans(my_lmer, pairwise ~ H2|Recombination|Map_type|QTL|Repulsion|Matrix)
  H2_summary <- as.data.frame(summary(my_emmeans$contrasts)) %>% filter(p.value < 0.05)
}


Map_type_summary <- cbind(response, Map_type_summary)
Recombination_summary <- cbind(response, Recombination_summary)
Matrix_summary <- cbind(response, Matrix_summary)
Repulsion_summary <- cbind(response, Repulsion_summary)
QTL_summary <- cbind(response, QTL_summary)
H2_summary <- cbind(response, H2_summary)
if(Full_founder == TRUE){
  QTL_type_summary <- cbind(response, QTL_type_summary)
}
if(QTL == TRUE | QTL_AF == TRUE){
  Allele_summary <- cbind(response, Allele_summary)
}

##### save output 
fwrite(Map_type_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                    df_name, "_", my_cycle, "_map_type.csv"))
fwrite(Recombination_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                df_name, "_", my_cycle, "_recombination.csv"))
fwrite(Matrix_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                     df_name, "_", my_cycle, "_matrix.csv"))
fwrite(Repulsion_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                              df_name, "_", my_cycle, "_repulsion.csv"))
fwrite(QTL_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                 df_name, "_", my_cycle, "_qtl.csv"))
fwrite(H2_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                           df_name, "_", my_cycle, "_H2.csv"))
if(Full_founder == TRUE){
  fwrite(QTL_type_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                  df_name, "_", my_cycle, "_qtl_type.csv"))
}
if(QTL == TRUE | QTL_AF == TRUE){
  fwrite(Allele_summary, paste0("/Users/ellietaagen/Desktop/github/cr_simulation/Results/master_linear_model_analysis_csv/full_founder/",
                                  df_name, "_", my_cycle, "_allele.csv"))
}

