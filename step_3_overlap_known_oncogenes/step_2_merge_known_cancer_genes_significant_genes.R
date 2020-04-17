library(tidyverse)
library(varhandle)
library(gsubfn)
#STK11

#############LOAD FILE################
df_upregulated = read_csv('upregulation.csv')
df_downregulated = read_csv('downregulation.csv')

ifname = 'known_cancer_genes.csv'

out_df = data.frame(n_signif=numeric(), n_known=numeric(), n_overlap=numeric())

#get list of cancer-specific tsgs and oncogenes
df_cancer_genes= read_csv(ifname) 
df_oncogenes = filter(df_cancer_genes, tumorsupressor_oncogene %in% c('oncogene', 'possible oncogene'))#uknown not included
df_tsg = filter(df_cancer_genes, tumorsupressor_oncogene %in% c('tsg', 'possible tsg'))

###############UPREGULATED################

#get ALL lung genes upreg
lung = cbind(stack(select(df_upregulated, c(Uall, Ucr_lung, Uen_lung, Ulung_genes))))
lung = unique(drop_na(lung))$values
known_Ulung = filter(df_oncogenes, Cancer=='LUAD')$Gene
lung_oncogenes = lung[lung %in% known_Ulung]#1 gene
sprintf('lung oncogenes. signif: %d, known: %d, overlap :%d', length(lung), length(known_Ulung), length(lung_oncogenes))
out_df['lung_oncogenes',] = c(length(lung), length(known_Ulung), length(lung_oncogenes))

#get ALL en genes
en = cbind(stack(select(df_upregulated, c(Uall, Ucr_en, Uen_lung, Uen_genes))))
en = unique(drop_na(en))$values
known_Uen = filter(df_oncogenes, Cancer=='UCEC')$Gene
en_oncogenes = en[en %in% known_Uen]#2 genes
sprintf('en oncogenes. signif: %d, known: %d, overlap :%d', length(en), length(known_Uen), length(en_oncogenes))
out_df['en_oncogenes',] = c(length(en), length(known_Uen), length(en_oncogenes))

#get ALL cr genes
cr = cbind(stack(select(df_upregulated, c(Uall, Ucr_en, Ucr_lung, Ucr_genes))))
cr = unique(drop_na(cr))$values
known_Ucr = filter(df_oncogenes, Cancer=='KIRC')$Gene
cr_oncogenes = cr[cr %in% known_Ucr]#1 gene
sprintf('cr oncogenes. signif: %d, known: %d, overlap :%d', length(cr), length(known_Ucr), length(cr_oncogenes))
out_df['cr_oncogenes',] = c(length(cr), length(known_Ucr), length(cr_oncogenes))

###############DOWNREGULATED################
#get ALL lung genes 
lungD = cbind(stack(select(df_downregulated, c(Dall, Dcr_lung, Den_lung, Dlung_genes))))
lungD = unique(drop_na(lungD))$values
known_Dlung = filter(df_tsg, Cancer=='LUAD')$Gene
lung_tsg = lungD[lungD %in% known_Dlung]#1 gene
sprintf('lung tsg. signif: %d, known: %d, overlap :%d', length(lungD), length(known_Dlung), length(lung_tsg))
out_df['lung_tsg',] = c(length(lungD), length(known_Dlung), length(lung_tsg))

#get ALL en genes
enD = cbind(stack(select(df_downregulated, c(Dall, Dcr_en, Den_lung, Den_genes))))
enD = unique(drop_na(enD))$values
known_Den = filter(df_tsg, Cancer=='UCEC')$Gene
en_tsg = enD[enD %in% known_Den]#1 gene
sprintf('en tsg. signif: %d, known: %d, overlap :%d', length(enD), length(known_Den), length(en_tsg))
out_df['en_tsg',] = c(length(enD), length(known_Den), length(en_tsg))

#get ALL cr genes
crD = cbind(stack(select(df_downregulated, c(Dall, Dcr_en, Dcr_lung, Dcr_genes))))
crD = unique(drop_na(crD))$values
known_Dcr = filter(df_tsg, Cancer=='KIRC')$Gene
cr_tsg = crD[crD %in% known_Dcr]#1 gene
sprintf('en tsg. signif: %d, known: %d, overlap :%d', length(crD), length(known_Dcr), length(cr_tsg))
out_df['cr_tsg',] = c(length(crD), length(known_Dcr), length(cr_tsg))

gene_names = list(lung_oncogenes, en_oncogenes, cr_oncogenes, lung_tsg, en_tsg, cr_tsg)
name_df = sapply(gene_names, '[', seq(max(lengths(gene_names))))
name_df = as.data.frame(name_df)
colnames(name_df) = c('lung_oncogenes','en_oncogenes','cr_oncogenes', 'lung_tsg','en_tsg','cr_tsg')

###########TO CSV##############
write.csv(out_df, 'number_signif_genes.csv')
write.csv(name_df, 'overlap_gene_names.csv')

