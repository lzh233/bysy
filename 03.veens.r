library(ggVennDiagram)
library(tidyverse)

dir.create("out/veen_AD_UB/plot/", recursive = T)
dir.create("out/veen_AD_UB/tables/", recursive = T)

make.degs.vene <- function(f1, f2, name1, name2) {
  df1 <- read_tsv(f1)
  df2 <- read_tsv(f2)
  out.list <- list()
  
  #df1.up <- df1 %>% filter(stat == "up") %>% .$gene
  #df2.up <- df2 %>% filter(stat == "up") %>% .$gene
  
  #df1.down <- df1 %>% filter(stat == "down") %>% .$gene
  #df2.down <- df2 %>% filter(stat == "down") %>% .$gene
  
  df1.all <- df1 %>% mutate(group = name1) %>% filter(stat %in% c("down", "up")) #%>% .$gene
  df2.all <- df2 %>% mutate(group = name2) %>% filter(stat %in% c("down", "up")) #%>% .$gene
  # up veen
  up.list <- list()
  up.list[[name1]] <- df1.all$gene
  up.list[[name2]] <- df2.all$gene
  
  p.all <- ggVennDiagram(up.list,edge_size = 1, set_size = 3)+
    scale_fill_distiller(palette = "RdBu")+
    scale_color_manual(values = c(rep("black",3)))
  
  # save data
  comm.genes <- intersect(df1.all$gene, df2.all$gene)
  print(comm.genes)
  comm.df <- rbind(
    df1.all %>% filter(gene %in% comm.genes),
    df2.all %>% filter(gene %in% comm.genes)
  )
  uniq.df1 <- df1.all %>% filter(!gene %in% comm.genes)
  uniq.df2 <- df2.all %>% filter(!gene %in% comm.genes)
  
  # out
  write_tsv(comm.df, str_glue("./out/veen_AD_UB/tables/common_{name1}_AND_{name2}.txt"))
  write_tsv(uniq.df1, str_glue("./out/veen_AD_UB/tables/unique_{name1}.txt"))
  write_tsv(uniq.df2, str_glue("./out/veen_AD_UB/tables/unique_{name2}.txt"))
  
  return(p.all)
}

p1 <- make.degs.vene(
  "out/deseq_tables/AD/AD_AD_T10_T12_vs_AD_T3_T5.txt", 
  "out/deseq_tables/UB/UB_UB_T10_T12_vs_UB_T3_T5.txt",
  "AD_T10_T12_vs_T3_T5",
  "UB_T10_T12_vs_T3_T5"
)


p2 <- make.degs.vene(
  "out/deseq_tables/AD/AD_AD_T10_T12_vs_AD_T8.txt",
  "out/deseq_tables/UB/UB_UB_T10_T12_vs_UB_T8.txt",
  "AD_T10_T12_vs_T8",
  "UB_T10_T12_vs_T8"
)

p3 <- make.degs.vene(
  "out/deseq_tables/AD/AD_AD_T8_vs_AD_T3_T5.txt",
  "out/deseq_tables/UB/UB_UB_T8_vs_UB_T3_T5.txt",
  "AD_T8_vs_T3_T5",
  "UB_T8_vs_T3_T5"
)



p <- list()
p[[1]] <- p1
p[[2]] <- p2
p[[3]] <- p3
p.out <- cowplot::plot_grid(plotlist = p,nrow = 3)
ggsave(filename = "out/veen_AD_UB/plot/veen.pdf",plot = p.out, width = 8, height = 8)







