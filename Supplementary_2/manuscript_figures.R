library(tidyverse)
library(data.table)
library(kableExtra)
library(ggsci)

source("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_1/master_functions/summary_stats.R")
gv <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gv.csv") %>% as.data.frame()
gg <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_gg.csv") %>% as.data.frame()
pa <- fread("https://raw.githubusercontent.com/etaagen/dissertation_chapter_4/main/Supplementary_2/results_S2.1/f_pa.csv")  %>% as.data.frame()
my_ASR_summary(measurement_df = (gg %>% filter(QTL == 200 & H2 == 8)), measurement = "genetic gain", df_name = "gg")


##### genetic variance, figure 3
# filter for 20X recombination, 200 QTL, H2 0.8, coupling / repulsion conditions 1 and 4

gv_summary <- gv %>% 
  filter(Recombination == "20X", QTL == 200, H2 == 8, Repulsion %in% c(1, 4)) %>% 
  group_by(cycle, Map_type, Matrix, Repulsion, QTL_type) %>% 
  summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
gv_summary$Matrix <- factor(gv_summary$Matrix, levels = c("GW", "CV"))
gv_summary$Map_type <- factor(gv_summary$Map_type, levels = c("WT", "Peri", "Chr"))
gv_summary$Repulsion <- factor(gv_summary$Repulsion, levels = c("1", "4"))
gv_summary$QTL_type <- factor(gv_summary$QTL_type, levels = c("R", "DV"))

gv_summary <- gv_summary %>% filter(!cycle == "founder")
gv_summary$cycle[gv_summary$cycle == "burnin"] <- 0
gv_summary$cycle <- as.numeric(gv_summary$cycle)
# jitter the cycle for easier viewing on plot
gv_summary$cycle <- as.numeric(gv_summary$cycle)
for(i in 1:nrow(gv_summary)){
  if(gv_summary$Map_type[i] == "WT"){
    gv_summary$cycle[i] <- gv_summary$cycle[i] - 0.05
  }else if(gv_summary$Map_type[i] == "Chr"){
    gv_summary$cycle[i] <- gv_summary$cycle[i] + 0.05
  }
}
gv_summary$cycle <- factor(gv_summary$cycle, 
                           levels = c(-0.05, 0, 0.05, 0.95,  1.0,  1.05,  1.95,  2.0,  2.05,  2.95,  3.0,  3.05,
                                      3.95,  4.0,  4.05,  4.95,  5.0,  5.05,  5.95,  6.0, 6.05,
                                      6.95,  7.0,  7.05,  7.95,  8.0,  8.05,  8.95,  9.0,  9.05,
                                      9.95, 10.0, 10.05))
# plots
gv_summary %>% 
  ggplot(aes(x = cycle, y = my_value, color = Map_type, shape = Repulsion)) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se, width=0.15)) +
  scale_color_jco(alpha = 0.6) + 
  scale_x_discrete(breaks=0:10)+
  #ggtitle("Genetic variance, R QTL") +
  ylab("") +
  xlab("") + 
  theme_minimal() + #labs(color='repulsion') +
  theme(text=element_text(size=18, face= "bold"), legend.position="bottom") +
  facet_grid(Matrix  ~ QTL_type)

ggsave("gv.svg", path= "/Users/ellietaagen/Desktop/", width = 7, height = 5)

##### genetic gain, figure 4
gg_summary <- gg %>% 
  filter(Recombination == "20X", QTL == 200, H2 == 8, Repulsion %in% c(1, 4)) %>% 
  group_by(cycle, Map_type, Matrix, Repulsion, QTL_type) %>%
  summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
gg_summary$Matrix <- factor(gg_summary$Matrix, levels = c("GW", "CV"))
gg_summary$Map_type <- factor(gg_summary$Map_type, levels = c("WT", "Peri", "Chr"))
gg_summary$Repulsion <- factor(gg_summary$Repulsion, levels = c("1", "4"))
gg_summary$cycle <- as.numeric(gg_summary$cycle)
gg_summary$QTL_type <- factor(gg_summary$QTL_type, levels = c("R", "DV"))
# jitter the cycle for easier viewing on plot
for(i in 1:nrow(gg_summary)){
  if(gg_summary$Map_type[i] == "WT"){
    gg_summary$cycle[i] <- gg_summary$cycle[i] - 0.05
  }else if(gg_summary$Map_type[i] == "Chr"){
    gg_summary$cycle[i] <- gg_summary$cycle[i] + 0.05
  }
}
gg_summary$cycle <- factor(gg_summary$cycle, 
                           levels = c(0.95,  1.0,  1.05,  1.95,  2.0,  2.05,  2.95,  3.0,  3.05,
                                      3.95,  4.0,  4.05,  4.95,  5.0,  5.05,  5.95,  6.0, 6.05,
                                      6.95,  7.0,  7.05,  7.95,  8.0,  8.05,  8.95,  9.0,  9.05,
                                      9.95, 10.0, 10.05))

gg_summary %>% 
  ggplot(aes(x = cycle, y = my_value, color = Map_type, shape = Repulsion)) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se, width=0.15)) +
  scale_color_jco(alpha = 0.6) + 
  scale_x_discrete(breaks=1:10)+
  # ggtitle("Genetic gain, R QTL") +
  ylab("") +
  xlab("") + 
  theme_minimal() + #+ labs(color='repulsion') + labs(color='repulsion') +
  theme(text=element_text(size=18, face= "bold"), legend.position="bottom") +
  facet_grid(Matrix  ~ QTL_type)

ggsave("gg.svg", path= "/Users/ellietaagen/Desktop/", width = 7, height = 5)




### prediction accuracy, figure 5
pa_summary <- pa %>% 
  filter(Recombination == "20X", QTL == 200, H2 == 8, Repulsion %in% c(1, 4)) %>% 
  group_by(cycle, Map_type, Matrix, Repulsion, QTL_type) %>%
  summarise(my_value = mean(value), sd = sd(value), se = sd(value)/sqrt(length(value)))
pa_summary$Matrix <- factor(pa_summary$Matrix, levels = c("GW", "CV"))
pa_summary$Map_type <- factor(pa_summary$Map_type, levels = c("WT", "Peri", "Chr"))
pa_summary$Repulsion <- factor(pa_summary$Repulsion, levels = c("1", "4"))
pa_summary$cycle <- as.numeric(pa_summary$cycle)
pa_summary$QTL_type <- factor(pa_summary$QTL_type, levels = c("R", "DV"))
# jitter the cycle for easier viewing on plot
for(i in 1:nrow(pa_summary)){
  if(pa_summary$Map_type[i] == "WT"){
    pa_summary$cycle[i] <- pa_summary$cycle[i] - 0.05
  }else if(pa_summary$Map_type[i] == "Chr"){
    pa_summary$cycle[i] <- pa_summary$cycle[i] + 0.05
  }
}
pa_summary$cycle <- factor(pa_summary$cycle, 
                           levels = c(0.95,  1.0,  1.05,  1.95,  2.0,  2.05,  2.95,  3.0,  3.05,
                                      3.95,  4.0,  4.05,  4.95,  5.0,  5.05,  5.95,  6.0, 6.05,
                                      6.95,  7.0,  7.05,  7.95,  8.0,  8.05,  8.95,  9.0,  9.05,
                                      9.95, 10.0, 10.05))

pa_summary %>% 
  ggplot(aes(x = cycle, y = my_value, color = Map_type, shape = Repulsion)) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin=value-se, ymax=value+se, width=0.15)) +
  scale_color_jco(alpha = 0.6) + 
  scale_x_discrete(breaks=1:10)+
  # ggtitle("Genetic gain, R QTL") +
  ylab("") +
  xlab("") + 
  theme_minimal() + #+ labs(color='repulsion') + labs(color='repulsion') +
  theme(text=element_text(size=18, face= "bold"), legend.position="bottom") +
  facet_grid(Matrix  ~ QTL_type)

ggsave("pa.svg", path= "/Users/ellietaagen/Desktop/", width = 7, height = 5)