library(ReactomePA)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)

add.id <- function(gene.table){
  eg <- bitr(
    geneID = gene.table$gene, 
    fromType = "SYMBOL",
    toType=c("ENTREZID"), 
    OrgDb="org.Hs.eg.db")
  colnames(eg) <- c("gene", "ENTREZID")
  gene.table <- gene.table %>% 
    left_join(eg,by = "gene")
  return(gene.table)
}

degs.df.AD <- read_tsv("./out/deseq_tables/AD/AD_AD_T10_T12_vs_AD_T3_T5.txt")
degs.df.UB <- read_tsv("./out/deseq_tables/UB/UB_UB_T10_T12_vs_UB_T3_T5.txt")

degs.TRRUST.AD <- degs.df.AD %>% filter(stat == "up") %>% select(gene)
write_tsv(degs.TRRUST.AD, "./data/pathway_data/TRRUST/TRRUST_gene_table_AD.txt", col_names = F)

degs.TRRUST.UB <- degs.df.UB %>% filter(stat == "up") %>% select(gene)
write_tsv(degs.TRRUST.UB, "./data/pathway_data/TRRUST/TRRUST_gene_table_UB.txt", col_names = F)

# read TRRUST result

AD.TRRUST <- read_tsv("out/pathways/TRRUST/AD_TRRUST_result.txt") %>% 
  mutate(log10padj = -log10(Q_value)) %>% 
  top_n(10, wt = log10padj) %>% 
  arrange(log10padj)

UB.TRRUST <- read_tsv("out/pathways/TRRUST/UB_TRRUST_result.txt") %>% 
  mutate(log10padj = -log10(Q_value))%>% 
  top_n(10, wt = log10padj) %>% 
  arrange(log10padj)

p1 <- ggplot(AD.TRRUST, aes(as_factor(Key_TF), log10padj)) + 
      geom_bar(stat = "identity") + 
      labs(x = "Key_TF") + 
      coord_flip() + ggtitle("AD") + theme_bw()

p2 <- ggplot(UB.TRRUST, aes(as_factor(Key_TF), log10padj)) + 
      geom_bar(stat = "identity") + 
      labs(x = "Key_TF") + 
      coord_flip() + ggtitle("UB")+ theme_bw()



# plot gene
fpkm <- read_tsv("data/fpkm_gene.txt")
gene_symbol <- read_tsv("data/gene_symbol.txt")
# AD
counts.AD <- fpkm %>% select(starts_with("AD_"), gene_id) %>% left_join(gene_symbol, by = "gene_id")

exp <- counts.AD %>% filter(gene %in% c("CDKN2A", "CCL2","IL6","BIRC3")) %>% 
  gather(key = "group", value = "count", -c(gene, gene_id)) %>% 
  mutate(log10count = log10(count + 1)) %>% 
  mutate(
    group2 = case_when(
      str_detect(group, "T10_T12") ~ "old",
      str_detect(group, "T3_T5") ~ "young"
    )
  ) %>% drop_na()

p.AD <- ggboxplot(exp, x = "group2", y = "log10count",
          fill = "group2", palette = "npg",
          add = "jitter") + facet_grid(~gene, scales = "free_x") + 
  stat_compare_means() + 
  theme_bw()+ggtitle("AD")
write_tsv(exp, "out/pathways/TRRUST/geneexp_AD_table.txt")

# UB

counts.UB <- fpkm %>% select(starts_with("UB_"), gene_id) %>% left_join(gene_symbol, by = "gene_id")

exp <- counts.UB %>% filter(gene %in% c("CDKN2A", "CCL2","IL6","BIRC3")) %>% 
  gather(key = "group", value = "count", -c(gene, gene_id)) %>% 
  mutate(log10count = log10(count + 1)) %>% 
  mutate(
    group2 = case_when(
      str_detect(group, "T10_T12") ~ "old",
      str_detect(group, "T3_T5") ~ "young"
    )
  ) %>% drop_na()

p.UB <- ggboxplot(exp, x = "group2", y = "log10count",
               # 配色方案 ?ggboxplot
               fill = "group2", palette = "npg",
               add = "jitter") + facet_grid(~gene, scales = "free_x") + 
  stat_compare_means() + 
  theme_bw()+ggtitle("UB")



p.trrust <- p1+p2
p.gene.list <- list()
p.gene.list[[1]] <- p.AD
p.gene.list[[2]] <- p.UB

p.gene <- cowplot::plot_grid(plotlist = p.gene.list,nrow = 2)


ggsave(filename = "./out/pathways/TRRUST/geneexp.png",plot = p.gene, height = 8.5, width = 10)
ggsave(filename = "./out/pathways/TRRUST/key_TF.png",plot = p.trrust)
write_tsv(exp, "out/pathways/TRRUST/geneexp_UB_table.txt")

# --GSEA ------------------------------------------------------------------


# ReactomePA
# AD
AD.genes <- degs.df.AD

AD.genes <- add.id(AD.genes)

AD.genes.use <- AD.genes %>% 
  drop_na() %>% 
  with(setNames(.$log2FoldChange, .$ENTREZID))


AD.Reactomeresult <- gsePathway(
  sort(AD.genes.use, decreasing = T), pvalueCutoff = 0.05)



ad.gsep <- enrichplot::gseaplot2(AD.Reactomeresult,
                      AD.Reactomeresult@result[["ID"]][1:3],
                      pvalue_table = T, title = "AD")

# UB

UB.genes <- degs.df.UB 

UB.genes <- add.id(UB.genes)

UB.genes.use <- UB.genes %>% 
  drop_na() %>% 
  with(setNames(.$log2FoldChange, .$ENTREZID))


UB.Reactomeresult <- gsePathway(
  sort(UB.genes.use, decreasing = T), pvalueCutoff = 0.05)



UB.gsep <- enrichplot::gseaplot2(UB.Reactomeresult,
                                 UB.Reactomeresult@result[["ID"]][1:3],
                                 pvalue_table = T, title = "UB")


# SENMAY
# AD
AD.genes.gsea <- AD.genes %>%
  with(setNames(.$log2FoldChange, .$gene))

AD.gsea <- GSEA(
  sort(AD.genes.gsea, decreasing = T),
  TERM2GENE = read.gmt("data/pathway_data/SAUL_SEN_MAYO.v2023.1.Hs.gmt")
)

AD.gsea.p <- enrichplot::gseaplot2(AD.gsea,
                      AD.gsea@result[["ID"]],
                      pvalue_table = T, title = "AD")



# UB
UB.genes.gsea <- UB.genes %>%
  with(setNames(.$log2FoldChange, .$gene))

UB.gsea <- GSEA(
  sort(UB.genes.gsea, decreasing = T),
  TERM2GENE = read.gmt("data/pathway_data/SAUL_SEN_MAYO.v2023.1.Hs.gmt")
)

UB.gsea.p <- enrichplot::gseaplot2(UB.gsea,
                                   UB.gsea@result[["ID"]],
                                   pvalue_table = T, title = "UB",)

# out
# ggsave(filename = "out/pathways/SENMAYO/SENMAYO_AD.png", plot = AD.gsea.p,device = "png" )
# ggsave(plot = UB.gsea.p, filename = "out/pathways/SENMAYO/SENMAYO_UB.png")

# ggsave(plot = ad.gsep, filename = "out/pathways/ReactomePA/ReactomePA_AD.png")
# ggsave(plot = UB.gsep, filename = "out/pathways/ReactomePA/ReactomePA_UB.png")


# out results
write_tsv(UB.Reactomeresult@result, "out/pathways/ReactomePA/ReactomePA_UB_table.txt")
write_tsv(AD.Reactomeresult@result, "out/pathways/ReactomePA/ReactomePA_AD_table.txt")

write_tsv(UB.gsea@result, "out/pathways/SENMAYO//SENMAYO_UB_table.txt")
write_tsv(AD.gsea@result, "out/pathways/SENMAYO//SENMAYO_AD_table.txt")



