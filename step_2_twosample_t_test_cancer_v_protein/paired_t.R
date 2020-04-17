# input: cancer and control ratio CSV for en, cccrcc, and lung cancer created 
#         in step 1 folder
# output: t-test result p-value, effect size, and fdr corrected p-value in 3 ttests.csv files


library(tidyverse)
library(reshape2)

# file path to CSV's from step 1
setwd("C:\\Users\\taylo\\Google Drive\\Capstone\\mRNA_Protein_Correlation\\csv's")

# read csv data
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



unique_genes <- unique(en_joined_filtered$gene)
test_sub <- unique_genes[1:50]
test_run <- en_joined_filtered %>% filter(gene %in% test_sub)


# function to run t-tests on newly tidy data
run_paired_t <- function(table) { 
  gene_list <- unique(table$gene)
  
  print_table <- c("gene", "p.value", "effect")
  for (g in gene_list) {
    print(g)
    sub <- table %>% subset(gene == g)
  
    normal <- na.omit(sub$value.con)
    cancer <- na.omit(sub$value.can)
    
  
    if (length(normal) <= 10 || length(cancer) <= 10) {
      #print(paste("Skipping: ",g))
      next
    }
    p <- t.test(normal, cancer, paired = T, na.rm = T)$p.value
    # positive m means cancer is upregulated- produces more protein
    m <- mean(cancer) - mean(normal)
  
    new_row <- c(g, p, m)
  
    print_table <- rbind(print_table, new_row)
  }
  
  write_matrix_fde <- print_table %>% as_data_frame()
  write_matrix_fde$fdr <- p.adjust(write_matrix_fde$V2, method = "fdr")
  final_write <- write_matrix_fde %>% filter(V1 != "gene") %>% rename(gene = V1, pval = V2, effect = V3)
  return (final_write)
}

test <- run_paired_t(test_run)

# run t-test for 3 cancer types and write to outfile 
en_complete <- run_paired_t(en_joined_filtered)
cr_complete <- run_paired_t(cr_joined_filtered)
lung_complete <- run_paired_t(lung_joined_filtered)

write.csv(en_complete, "en_ttests.csv", row.names = F)
write.csv(cr_complete, "cr_ttests.csv", row.names = F)
write.csv(lung_complete, "lung_ttests.csv", row.names = F)
