#input: csv of dimensions (n_signif_genes, n_cancer_types)
#output: heatmap using the first 20 columns
library(tidyverse)
library(vroom)
library(RColorBrewer)
library(cowplot)

cr <- read.csv("cr_ttests.csv", header = T) %>% filter(fdr < .05) %>% rename(renal = effect) %>% select(gene, renal)
en <- read.csv("en_ttests.csv", header = T) %>% filter(fdr < .05) %>% rename(endometrial = effect) %>% select(gene, endometrial)
lung <- read.csv("lung_ttests.csv", header = T) %>% filter(fdr < .05) %>% rename(lung = effect) %>% select(gene, lung)

merged <- lung %>% full_join(cr, by = "gene") %>% full_join(en, by = "gene")

####gene selection. Choosing top 20 genes for which there is data for all 3 types, and who have the highest difference in comparison
merged %>% filter(!is.na(lung)& !is.na(renal) & !is.na(endometrial)) %>% mutate(lung_norm = lung/sd(lung), renal_norm = renal/sd(renal), endometrial_norm = endometrial/sd(endometrial), mean = (lung_norm + renal_norm + endometrial_norm)/3) -> all_norm_mean

####find greatest magnitude positive and negative different, but for genes with same direction for all 3
all_norm_mean %>% filter(lung_norm > 0 & renal_norm > 0 & endometrial_norm > 0) %>% arrange(desc(mean)) %>% slice(1:10) -> top_10
all_norm_mean %>% filter(lung_norm < 0 & renal_norm < 0 & endometrial_norm < 0) %>% arrange(mean) %>% slice(1:10) -> bottom_10

top_10 %>% select(gene, lung_norm, renal_norm, endometrial_norm) -> relevant_top
bottom_10 %>% select(gene, lung_norm, renal_norm, endometrial_norm) -> relevant_bottom  

final_list <- rbind(relevant_top, relevant_bottom)

write.csv(final_list, "figure2_list.csv", row.names = F)


test=F
if (test){
  tbl = gen_tbl(rows=20, cols=3, col_types='ddd')
  tbl = cbind(sprintf('gene%d',1:20), tbl)
  colnames(tbl) = c('gene','lung','renal','endometrial') 
} else{
  filename='figure2_list.csv'
  tbl = read.csv(filename , header = T)
}

####################MELT DF####################
#dim (n_genes*n_cancer_types, 2)
tbl <- tbl %>% rename(lung = lung_norm, renal = renal_norm, endometrial = endometrial_norm) %>% filter(!gene %in% c("ALDOB", "VIM"))
cancer_types = c('lung','renal','endometrial')
df_melt = pivot_longer(data=tbl, cols=cancer_types, names_to='cancer_type', values_to='mRNA_to_prot')

##########align legend################
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

top2=round(min(df_melt$mRNA_to_prot)*2)/2
bot2=round(max(df_melt$mRNA_to_prot)*2)/2
scale=(bot2:top2)
p=ggplot(df_melt, aes(x=gene, y=cancer_type, fill=mRNA_to_prot)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2('Difference in \nmRNA/Protein ratio', breaks=scale, low='blue', mid='white', high='red') +
  ylab('Cancer Type') + 
  xlab('Gene') +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, size=8, vjust=-.05), legend.text = element_text(size = 8), legend.title.align=0.5, legend.box.just = "center", axis.ticks = element_blank())
ggdraw(align_legend(p))
#+
  #title('mRNA to Protein Difference Cancer ')

ggsave("figure2.png", height = 5, width = 10)
