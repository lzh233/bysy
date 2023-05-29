library(tidyverse)
library(DESeq2)
library(ggrepel)
#library(ComplexHeatmap)

dir.create(str_glue("out/deseq_tables/"),recursive = T)
dir.create(str_glue("out/deseq_plots/"), recursive = T)

plot_theme <- 
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 13,color = "black"))


df <- read.table("data/raw_gene.txt", header = T, row.names = 1)
gene.ann <- read.table("data/gene_symbol.txt", header = T)
AD.gr <- read_tsv("data/AD_gr.txt")
UB.gr <- read_tsv("data/UB_gr.txt")

p <- 0.05
fc <- 2
lable.num <- 10


for (ori in c("UB")) {
  degs.all <- list()
  # AD or UB
  data.use <- ori

  if (data.use == "AD") {
    use <- AD.gr
  }else if (data.use == "UB") {
    use <- UB.gr
  }
  
  dir.create(str_glue("out/deseq_tables/{data.use}"))
  dir.create(str_glue("out/deseq_plots/{data.use}"))
  
  for (v1 in seq_along(use$col1)) {
    v1.gr <- use[v1, 1]$col1
    v2.gr <- use[v1, 2]$col2
    df.use <- df %>% 
      select(names(.)[str_detect(names(.), v1.gr) | str_detect(names(.), v2.gr)]) %>% 
      filter(rowSums(.) > 0)
    
    group.list <- gsub("\\_\\d+$", "", names(df.use))
    print(str_c(group.list, collapse = "    "))
    
    # deseq
    colData <- data.frame(row.names=colnames(df.use), 
                          group_list=factor(group.list, levels = c(v1.gr, v2.gr)))
    
    print(colData)
    print(str_c(colnames(df.use), collapse = "    "))
    
    dds <- DESeqDataSetFromMatrix(countData = df.use,
                                  colData = colData,
                                  design = ~ group_list)
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds, contrast=c("group_list",v1.gr, v2.gr)))
    
    # drop NA
    res <- res %>% 
      as.data.frame(.) %>% 
      filter(!is.na(padj)) %>% 
      mutate(gene_id = rownames(.)) %>% 
      left_join(gene.ann, by = "gene_id")
      
    # plot VOL
    res <- res %>% mutate(
      stat =  case_when(
        (padj <= p & log2FoldChange >= fc) ~ "up",
        (padj <= p & log2FoldChange <= -fc) ~ "down",
        TRUE ~ "None"
      )
    )
    
    # marks table
    # fc
    label.fc <- res %>% 
      filter(stat %in% c("up", "down")) %>%
      group_by(stat) %>% 
      top_n(round(lable.num/2, 0),wt = abs(log2FoldChange)) %>% 
      ungroup()
    
    label.p <- res %>% 
      filter(stat %in% c("up", "down")) %>%
      group_by(stat) %>% 
      top_n(round(lable.num/2, 0),wt = -padj)%>% 
      ungroup()
    
    lable.df <- rbind(label.p, label.fc) %>% distinct()
  
    p.vol <- ggplot(res, aes(log2FoldChange, -log10(padj), color = stat)) + 
      geom_point(size = 0.7) + 
      geom_label_repel(lable.df, mapping = aes(log2FoldChange, -log10(padj),label = gene),show.legend = F)+
      scale_color_manual(values = c(up = "#F26522", down = "#00AEEF", None ="gray")) +
      geom_hline(aes(yintercept=-log10(p)), linewidth = 1, linetype="dashed", color = "gray70") + 
      geom_vline(aes(xintercept=fc), linewidth = 0.5, linetype="dashed", color = "gray70") + 
      geom_vline(aes(xintercept=-fc), linewidth = 0.5, linetype="dashed", color = "gray70") + 
      ggtitle(str_glue("{v1.gr}.vs.{v2.gr}"))+
      theme_bw() + plot_theme
   
    
    res.out <- res %>% 
      select(gene_id, gene, log2FoldChange, padj, stat) %>% 
      mutate(log2FoldChange = round(log2FoldChange, 4),
             padj  = round(padj, 4))
    
    degs.all[[v1]] <- res %>% mutate(group = str_glue("{v1.gr}_vs_{v2.gr}"))
    
    write_tsv(res.out,
              str_glue("out/deseq_tables/{data.use}/{data.use}_{v1.gr}_vs_{v2.gr}.txt")
      ) 
    
    ggsave(filename = str_glue("out/deseq_plots/{data.use}/{data.use}_{v1.gr}_vs_{v2.gr}.png"), 
           plot = p.vol, 
           device = "png",
           width = 8,
           height = 8)
    
    use.genes <- res %>% filter(stat %in% c("up", "down")) %>% .$gene_id
    count.data <- df.use[rownames(df.use) %in% use.genes, ]
    
    count.data <- count.data %>% 
      mutate(gene_id = rownames(.)) %>% 
      left_join(gene.ann, by = "gene_id")
    
    write_tsv(count.data, 
                str_glue("out/deseq_tables/{data.use}/{data.use}_{v1.gr}_vs_{v2.gr}_count.txt"))
    write_tsv(lable.df, 
              str_glue("out/deseq_tables/{data.use}/{data.use}_{v1.gr}_vs_{v2.gr}_top.txt"))
    
  }
  degs.all.df <- do.call("rbind", degs.all)
  
  degs.count.df <- degs.all.df %>% 
    group_by(group, stat) %>% summarise(ngenes = n()) %>%
    filter(!stat %in% "None")
  
  p <- ggplot(degs.count.df, aes(group, ngenes, fill = stat)) + 
    geom_bar(stat = "identity",
             position = position_dodge(0.75),width = 0.75) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  write_tsv(degs.count.df,
            str_glue("out/deseq_tables/{data.use}/{data.use}_DEGs_count.txt"))
  
  ggsave(filename = str_glue("out/deseq_plots/{data.use}/{data.use}_DEGs.png"), 
                   plot = p, 
                   device = "png",
                   width = 8,
                   height = 5)
  degs.all <- list()
  
}

