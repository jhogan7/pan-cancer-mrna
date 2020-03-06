#input: csv of dimensions (n_signif_genes, n_cancer_types)
#output: heatmap using the first 20 columns
library(tidyverse)
library(vroom)
library(RColorBrewer)

test=T
if (test){
  tbl = gen_tbl(rows=20, cols=3, col_types='ddd')
  tbl = cbind(sprintf('gene%d',1:20), tbl)
  colnames(tbl) = c('gene','colon','renal','endometrial') 
} else{
  filename='TODO_replaceme.csv'
  tbl = read_csv(filename)
}

####################MELT DF####################
#dim (n_genes*n_cancer_types, 2)
cancer_types = c('colon','renal','endometrial')
df_melt = pivot_longer(data=tbl, cols=cancer_types, names_to='cancer_type', values_to='mRNA_to_prot')

top2=round(min(df_melt$mRNA_to_prot)*2)/2
bot2=round(max(df_melt$mRNA_to_prot)*2)/2
scale=(bot2:top2)
ggplot(df_melt, aes(x=gene, y=cancer_type, fill=mRNA_to_prot)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient2('Difference in \nRNA/Protein ratio', breaks=scale, low='blue', mid='white', high='red') +
  ylab('Cancer Type') + 
  xlab('Gene') +
  theme(axis.text.x = element_text(angle = 45, size=8), legend.title.align=0.5) #+
  #title('mRNA to Protein Difference Cancer ')
