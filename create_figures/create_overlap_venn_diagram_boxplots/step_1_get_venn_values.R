#input: ttest csv for each cancer type with a p value, effect size, and fdr correction 
#        for all paired comparisons between protein/mRNA ratios for each gene from paired_t.R
#output: counts of gene numbers shared by different cancer types ordered to be input in
#         venn_diagrams.ipynb for Figure 2
library(tidyverse)
library(reshape2)

# get lists of all genes with valid data present in each cancer types
setwd("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's")

en_can <- read.csv("en_cancer.csv", header = T)
en_con <- read.csv("en_control.csv", header = T)
lung_can <- read.csv("lung_cancer.csv", header = T)
lung_con <- read.csv("lung_control.csv", header = T)
cr_can <- read.csv("cr_cancer.csv", header = T)
cr_con <- read.csv("cr_control.csv", header = T)


####melting gene columns into single categorical variable and combining sample and gene into single unqiue id to join
melt_can <- melt(en_can, id.vars = "X") %>% rename(value.can = value) %>% unite(ID, X, variable, sep = "_")
melt_con <- melt(en_con, id.vars = "X") %>% rename(value.con = value) %>% unite(ID, X, variable, sep = "_")

####join by unique ID, separate back into sample and gene, and filter NA's (unpaired)
joined <- full_join(melt_can, melt_con, by = "ID")
joined_final <- joined %>% separate(ID, c("sample", "gene"), sep = "_")
en_joined_filtered <- joined_final %>% filter(!is.na(value.con) & !is.na(value.can))

####melting gene columns into single categorical variable and combining sample and gene into single unqiue id to join
melt_can <- melt(lung_can, id.vars = "X") %>% rename(value.can = value) %>% unite(ID, X, variable, sep = "_")
melt_con <- melt(lung_con, id.vars = "X") %>% rename(value.con = value) %>% unite(ID, X, variable, sep = "_")

####join by unique ID, separate back into sample and gene, and filter NA's (unpaired)
joined <- full_join(melt_can, melt_con, by = "ID")
joined_final <- joined %>% separate(ID, c("sample", "gene"), sep = "_")
lung_joined_filtered <- joined_final %>% filter(!is.na(value.con) & !is.na(value.can))
lung_joined_filtered$value.can <- as.numeric(lung_joined_filtered$value.can)
lung_joined_filtered$value.con <- as.numeric(lung_joined_filtered$value.con)

####melting gene columns into single categorical variable and combining sample and gene into single unqiue id to join
melt_can <- melt(cr_can, id.vars = "X") %>% rename(value.can = value) %>% unite(ID, X, variable, sep = "_")
melt_con <- melt(cr_con, id.vars = "X") %>% rename(value.con = value) %>% unite(ID, X, variable, sep = "_")

####join by unique ID, separate back into sample and gene, and filter NA's (unpaired)
joined <- full_join(melt_can, melt_con, by = "ID")
joined_final <- joined %>% separate(ID, c("sample", "gene"), sep = "_")
cr_joined_filtered <- joined_final %>% filter(!is.na(value.con) & !is.na(value.can))
cr_joined_filtered$value.can <- as.numeric(as.character(cr_joined_filtered$value.can))
cr_joined_filtered$value.con <- as.numeric(as.character(cr_joined_filtered$value.con))



get_gene_list <- function(table) { 
  gene_list_from_table <- unique(table$gene)
  used_genes = c()
  
  for (g in gene_list_from_table) {
    sub <- table %>% subset(gene == g)
    
    normal <- na.omit(sub$value.con)
    cancer <- na.omit(sub$value.can)
    
    
    if (length(normal) <= 10 || length(cancer) <= 10) {
      # print(paste("Skipping: ",g))
      next
    }
    
    used_genes = c(used_genes, g)
  }
  
  return (used_genes)
}


en_gene_list <- get_gene_list(en_joined_filtered)
cr_gene_list <- get_gene_list(cr_joined_filtered)
lung_gene_list <- get_gene_list(lung_joined_filtered)


# make list of genes with valid data for all three cancer types
all_genes = c(cr_gene_list, en_gene_list, lung_gene_list)
gene_list = unique(all_genes[all_genes %in% cr_gene_list & all_genes %in% en_gene_list & all_genes %in% lung_gene_list])

minef = 0
# create upregulated lists of genes for each cancer type
ucr <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\cr_ttests.csv", header = T) %>% filter(fdr < .05, effect > minef))$gene
ucr <- ucr[ucr %in% gene_list]
uen <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\en_ttests.csv", header = T) %>% filter(fdr < .05, effect > minef))$gene
uen <- uen[uen %in% gene_list]
ulung <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\lung_ttests.csv", header = T) %>% filter(fdr < .05, effect > minef))$gene
ulung <- ulung[ulung %in% gene_list]


# create downregulated lists of genes for each cancer type
dcr <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\cr_ttests.csv", header = T) %>% filter(fdr < .05, effect < minef))$gene
dcr <- dcr[dcr %in% gene_list]
den <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\en_ttests.csv", header = T) %>% filter(fdr < .05, effect < minef))$gene
den <- den[den %in% gene_list]
dlung <- (read.csv("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's\\lung_ttests.csv", header = T) %>% filter(fdr < .05, effect < minef))$gene
dlung <- dlung[dlung %in% gene_list]

# find upregulated lists of unique and shared genes by cancer types
uall = intersect(intersect(ucr, uen), ulung)
ucr_lung = setdiff(intersect(ucr, ulung), uall)
ucr_en = setdiff(intersect(ucr, uen), uall)
uen_lung = setdiff(intersect(uen, ulung), uall)
uen_genes = setdiff(setdiff(uen, ucr), ulung)
ucr_genes = setdiff(setdiff(ucr, uen), ulung)
ulung_genes = setdiff(setdiff(ulung, ucr), uen)

# find downregulated lists of unique and shared genes by cancer types
dall = intersect(intersect(dcr, den), dlung)
dcr_lung = setdiff(intersect(dcr, dlung), dall)
dcr_en = setdiff(intersect(dcr, den), dall)
den_lung = setdiff(intersect(den, dlung), dall)
den_genes = setdiff(setdiff(den, dcr), dlung)
dcr_genes = setdiff(setdiff(dcr, den), dlung)
dlung_genes = setdiff(setdiff(dlung, dcr), den)

# print numbers of genes 
print(c(length(ucr_genes), length(ulung_genes), length(ucr_lung), length(uen_genes), length(ucr_en), length(uen_lung), length(uall)))
print(c(length(dcr_genes), length(dlung_genes), length(dcr_lung), length(den_genes), length(dcr_en), length(den_lung), length(dall)))

