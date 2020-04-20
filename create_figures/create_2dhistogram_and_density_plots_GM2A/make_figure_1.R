library(tidyverse)
library(reshape2)

setwd("BIO465/Correlations")

#THESE ARE THE RESULTS OF STEP1 CALCULATE PMR
en_can <- read.csv("en_cancer.csv", header = T)
en_con <- read.csv("en_control.csv", header = T) 

ren_con <- read.csv("cr_control.csv", header = T)
ren_can <- read.csv("cr_cancer.csv", header = T)

lung_con <- read.csv("lung_control.csv", header = T)
lung_can <- read.csv("lung_cancer.csv", header = T)

####melting gene columns into single categorical variable and combining sample and gene into single unqiue id to join
melt_can <- melt(ren_can, id.vars = "X") %>% rename(value.can = value) %>% unite(ID, X, variable, sep = "_")
melt_con <- melt(ren_con, id.vars = "X") %>% rename(value.con = value) %>% unite(ID, X, variable, sep = "_")

####join by unique ID, separate back into sample and gene, and filter NA's (unpaired)
joined <- full_join(melt_can, melt_con, by = "ID")
joined_final <- joined %>% separate(ID, c("sample", "gene"), sep = "_")
joined_filtered <- joined_final %>% filter(!is.na(value.con) & !is.na(value.can))

####summarizes mean for control and cancer for each gene
joined_filtered %>% filter(!is.na(value.can) | !is.na(value.con)) %>% group_by(gene) %>% mutate(n = n()) %>% filter(n > 10) -> trimmed  


master_graph <- trimmed %>% group_by(gene) %>% dplyr::summarize(av_can = mean(value.can, na.rm = T), av_con = (mean(value.con, na.rm = T)))




###scatter plot
ggplot(master_graph, aes(av_con, av_can)) + geom_point() + theme_bw()
ggsave("raw_scatter.png", height = 5, width = 5)

###2d histgram ("binned" scatterplot)
ggplot(master_graph, aes(av_con, av_can)) + geom_bin2d(bins = 100) + geom_smooth(method = "lm", se = T, col = "red", size = .75, linetype = "dashed") + theme_bw() + theme(text = element_text(size = 18)) + labs(x = "Average Control Ratio", y = "Average Cancer Ratio") + xlim(.5, 2) + ylim(.5, 2) #+geom_segment(aes(x = .5, y = .5, xend = 1.7, yend = 1.7), col = "red", size = .75) + geom_vline(xintercept = 1.67) + geom_hline(yintercept = .905)
ggsave("ren_2d_post_shift.png", height = 5, width = 7)

#### determine the linear regression formula
summary(lm(av_can~av_con, master_graph))

#### tests whether regression line is significantly different from y = x
summary(lm(av_can~av_con, master_graph, offset = 1*av_con))
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.10476    0.00259   40.45   <2e-16 ***
#  av_con       0.85977    0.00353  243.58   <2e-16 ***
  ---
ggsave("2d_hist_fit.png", height = 5, width = 7)

###eyeballed cuttoffs by looked at gene in lower right quadrant, easilt disinguishable from other genes
eg <- master_graph %>% filter(av_con > 1.5 & av_can < 1)
head(eg)
 ### extractable gene @ "EBP (1.49, .617)"


###returning to original melted frames to graph distribution of single gene from cancer and control
con_ex <- melt_con %>% separate(ID, c("sample", "gene"), sep = "_") %>%  filter(gene == "GM2A")
can_ex <- melt_can %>% separate(ID, c("sample", "gene"), sep = "_") %>%  filter(gene == "GM2A")


###hist control
###vertical lines are mean values
###custom xlim coerces both graphs to use same scale to easily compare
ggplot(con_ex, aes(value.con)) + geom_density() + theme_bw() + theme(text = element_text(size = 18), plot.title = element_text(hjust = .5)) + labs(x = "mRNA/Protein Ratio", y = "Density")+ ggtitle("Control") + geom_vline(xintercept = 1.67, linetype = "dashed", color = "red") + xlim(.5, 2) + ylim(0, 3.6)
ggsave("con_GM2A_density.png", height = 5, width = 5) 

###hist cancer
ggplot(can_ex, aes(value.can)) + geom_density() + theme_bw() + theme(text = element_text(size = 18), plot.title = element_text(hjust = .5)) + labs(x = "mRNA/Protein Ratio", y = "Density")+ ggtitle("Cancer") + geom_vline(xintercept = .905, linetype = "dashed", color = "red") + xlim(.5, 2) + ylim(0, 3.6)
ggsave("can_GM2A_density.png", height = 5, width = 5) 

