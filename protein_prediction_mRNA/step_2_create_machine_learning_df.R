library(tidyverse)
library(ggplot2)

##############CREATE ML DATAFRAME#################
#input: multi-indexed dataframe of dimensions (n_patients * n_genes, 3)
  #rows: 1 for each patient_id and gene
  #columns: patient_id, protein, mRNA
#output: dataframe of dimensions (n_patients*n_signif_genes, sum(n_rbp_genes+3))
  #rows: 1 for each patient_id and gene of interest
  #columns: patient_id, protein values for all RNA Binding Proteins (RBP), mRNA value gene x, protein value gene x
  #x is a gene for each gene with significantly different mRNA ratios between cancer and control 
  #NOTE: the last column (protein value gene x) is our y vector, all other values are X (explanatory variable matrix) 

##############CREATE COLON DATAFRAME#################
colon_df = read_csv('colon_dataset_df.csv') %>% select(-X1)

#TODO select only genes that are RBP
out_df = pivot_wider(colon_df, values_from=c(Gene_Name, protein), names_from=Patient_ID)

endometrial_df = read_csv('endometrial_dataset_df.csv')

