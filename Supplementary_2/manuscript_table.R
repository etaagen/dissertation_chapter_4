### Table 1

library(tidyverse)
library(data.table)
library(kableExtra)
library(ggsci)
library(gt)
library(lme4)
library(emmeans)
library(car)
library(broom.mixed)

gv <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gv.csv") %>% as.data.frame()
gg <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gg.csv") %>% as.data.frame()

## genetiv variance 
my_data <- gv
my_data$rep <- as.factor(my_data$rep)
my_data$Map_type <- as.factor(my_data$Map_type)
my_data$Recombination<- as.factor(my_data$Recombination)
my_data$QTL<- as.factor(my_data$QTL)
my_data$H2<- as.factor(my_data$H2)
my_data$Repulsion<- as.factor(my_data$Repulsion)
my_data$Matrix<- as.factor(my_data$Matrix)
my_data$QTL_type <- as.factor(my_data$QTL_type)

# filter to 6th and 10th cycle
my_data_6 <- my_data %>% filter(cycle == 6)
my_data_10 <- my_data %>% filter(cycle == 10)

# run mixed model 
my_lmer_6 <- lmer(value ~ (Map_type + Recombination+QTL+H2+Repulsion+Matrix+QTL_type)^2 + (1|rep), data = my_data_6)
my_lmer_10 <- lmer(value ~ (Map_type + Recombination+QTL+H2+Repulsion+Matrix+QTL_type)^2 + (1|rep), data = my_data_10)

# run Anova
my_Anova_6 <- anova(my_lmer_6) 
my_Anova_10 <- anova(my_lmer_10) 

afss6 <- my_Anova_6$"Sum Sq"
afss10 <- my_Anova_10$"Sum Sq"

cycle6 <- cbind(my_Anova_6,PctExp=afss6/sum(afss6)*100)
cycle6 <-rownames_to_column(cycle6)
cycle10 <- cbind(my_Anova_10,PctExp=afss10/sum(afss10)*100)
cycle10 <-rownames_to_column(cycle10)

write_csv(cycle6, "/Users/ellietaagen/Desktop/gv_6.csv")
write_csv(cycle10, "/Users/ellietaagen/Desktop/gv_10.csv")


### genetic gain
my_data <- gg
my_data$rep <- as.factor(my_data$rep)
my_data$Map_type <- as.factor(my_data$Map_type)
my_data$Recombination<- as.factor(my_data$Recombination)
my_data$QTL<- as.factor(my_data$QTL)
my_data$H2<- as.factor(my_data$H2)
my_data$Repulsion<- as.factor(my_data$Repulsion)
my_data$Matrix<- as.factor(my_data$Matrix)
my_data$QTL_type <- as.factor(my_data$QTL_type)

# filter to 6th and 10th cycle
my_data_6 <- my_data %>% filter(cycle == 6)
my_data_10 <- my_data %>% filter(cycle == 10)

# run mixed model 
my_lmer_6 <- lmer(value ~ (Map_type + Recombination+QTL+H2+Repulsion+Matrix+QTL_type)^2 + (1|rep), data = my_data_6)
my_lmer_10 <- lmer(value ~ (Map_type + Recombination+QTL+H2+Repulsion+Matrix+QTL_type)^2 + (1|rep), data = my_data_10)

# run Anova
my_Anova_6 <- anova(my_lmer_6) 
my_Anova_10 <- anova(my_lmer_10) 

afss6 <- my_Anova_6$"Sum Sq"
afss10 <- my_Anova_10$"Sum Sq"

cycle6 <- cbind(my_Anova_6,PctExp=afss6/sum(afss6)*100)
cycle6 <-rownames_to_column(cycle6)
cycle10 <- cbind(my_Anova_10,PctExp=afss10/sum(afss10)*100)
cycle10 <-rownames_to_column(cycle10)

write_csv(cycle6, "/Users/ellietaagen/Desktop/gg_6.csv")
write_csv(cycle10, "/Users/ellietaagen/Desktop/gg_10.csv")
