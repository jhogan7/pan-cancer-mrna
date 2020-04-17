# input: cancer and control ratio CSV for en, cccrcc, and lung cancer created 
#         in step 1 folder
# output: file with original ratio values and effect sizes for reference

library(tidyverse)
library(reshape2)

# read in CSV's
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

# make list of genes with valid data for all three cancer types
all_genes = c(cr_gene_list, en_gene_list, lung_gene_list)
gene_list = unique(all_genes[all_genes %in% cr_gene_list & all_genes %in% en_gene_list & all_genes %in% lung_gene_list])

# create rows for every gene with at least 10 samples in all cancer types
get_data_frame <- function(e_table, c_table, l_table) { 
  e_gene_list <- unique(e_table$gene)
  c_gene_list <- unique(c_table$gene)
  l_gene_list <- unique(l_table$gene)
  
  print_table <- c("gene", "e_n", "e_c", "e_effect", "l_n", "l_c", "l_effect", "c_n", "c_c", "c_effect")
  for (g in gene_list) {
    print(g)
    e_sub <- e_table %>% subset(gene == g)
    c_sub <- c_table %>% subset(gene == g)
    l_sub <- l_table %>% subset(gene == g)
    
    e_normal <- na.omit(e_sub$value.con)
    e_cancer <- na.omit(e_sub$value.can)
    c_normal <- na.omit(c_sub$value.con)
    c_cancer <- na.omit(c_sub$value.can)
    l_normal <- na.omit(l_sub$value.con)
    l_cancer <- na.omit(l_sub$value.can)
    
    
    if (length(e_normal) <= 10 || length(e_cancer) <= 10) {
      print(paste("Skipping: ",g))
      next
    }
    
    e_e <- mean(e_cancer) - mean(e_normal)
    l_e <- mean(l_cancer) - mean(l_normal)
    c_e <- mean(c_cancer) - mean(c_normal)
    
    new_row <- c(g, mean(e_normal), mean(e_cancer), e_e, mean(l_normal), mean(l_cancer), l_e, mean(c_normal), mean(c_cancer), c_e)
    
    print_table <- rbind(print_table, new_row)
  }
  
  return (print_table)
}


# create a column of absolute value of all effect sizes combined
table_complete <- get_data_frame(en_joined_filtered, cr_joined_filtered, lung_joined_filtered)
table_names = table_complete[1,]
table_effect = as_tibble(table_complete[-1,])
names(table_effect) = table_names
effect_final = table_effect %>% mutate(total_effect = abs(as.numeric(e_effect)) + abs(as.numeric(l_effect)) + abs(as.numeric(c_effect)))

write.csv(effect_final, "ratio_values.csv", row.names = T)

