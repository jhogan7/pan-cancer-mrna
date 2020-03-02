library(tidyverse)
library(reshape2)

setwd("BIO465/Correlations")

en_can <- read.csv("en_cancer.csv", header = T)
en_con <- read.csv("en_control.csv", header = T) 

####melting gene columns into single categorical variable and combining sample and gene into single unqiue id to join
melt_can <- melt(en_can, id.vars = "X") %>% rename(value.can = value) %>% unite(ID, X, variable, sep = "_")
melt_con <- melt(en_con, id.vars = "X") %>% rename(value.con = value) %>% unite(ID, X, variable, sep = "_")

####join by unique ID, separate back into sample and gene, and filter NA's (unpaired)
joined <- full_join(melt_can, melt_con, by = "ID")
joined_final <- joined %>% separate(ID, c("sample", "gene"), sep = "_")
joined_filtered <- joined_final %>% filter(!is.na(value.con) & !is.na(value.can))

####summarizes mean for control and cancer for each gene
master_graph <- joined_filtered %>% group_by(gene) %>% dplyr::summarize(av_can = mean(value.can, na.rm = T), av_con = (mean(value.con, na.rm = T)))


###scatter plot
ggplot(master_graph, aes(av_con, av_can)) + geom_point() + theme_bw()
ggsave("raw_scatter.png", height = 5, width = 5)

###2d histgram ("binned" scatterplot)
ggplot(master_graph, aes(av_con, av_can)) + geom_bin2d(bins = 100) + geom_smooth(method = "lm", se = T, col = "red", size = .75, linetype = "dashed") + theme_bw() + theme(text = element_text(size = 18)) + labs(x = "Average Control Ratio", y = "Average Cancer Ratio") + geom_segment(aes(x = 0, y = 0, xend = 1.3, yend = 1.3), col = "red", size = .75)

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
eg <- master_graph %>% filter(av_con > .75, av_can < .45)
head(eg)
 ### extractable gene @ "DPT"


###returning to original melted frames to graph distribution of single gene from cancer and control
con_ex <- melt_con %>% separate(ID, c("sample", "gene"), sep = "_") %>%  filter(gene == "DPT")
can_ex <- melt_can %>% separate(ID, c("sample", "gene"), sep = "_") %>%  filter(gene == "DPT")


###hist control
###vertical lines are mean values
###custom xlim coerces both graphs to use same scale to easily compare
ggplot(con_ex, aes(value.con)) + geom_density() + theme_bw() + theme(text = element_text(size = 18), plot.title = element_text(hjust = .5)) + labs(x = "mRNA/Protein Ratio", y = "Density")+ ggtitle("Control") + xlim(0, .9) + ylim(0, 7)  + geom_vline(xintercept = .772, linetype = "dashed", color = "red")
ggsave("con_DPT_density.png", height = 5, width = 5) 

###hist cancer
ggplot(can_ex, aes(value.can)) + geom_density() + theme_bw() + theme(text = element_text(size = 18), plot.title = element_text(hjust = .5)) + labs(x = "mRNA/Protein Ratio", y = "Density")+ ggtitle("Cancer") + geom_vline(xintercept = .441, linetype = "dashed", color = "red") + xlim(0, .9) + ylim(0, 7)
ggsave("can_DPT_density.png", height = 5, width = 5) 

###t.test for DPT cancer vs control
t.test(con_ex$value.con, can_ex$value.can, paired = T)
#t = 6.0375, df = 13, p-value = 4.183e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.2125384 0.4493941
#sample estimates:
#  mean of the differences: 0.3309663
