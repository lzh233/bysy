library(UpSetR)
library(tidyverse)
library(patchwork)

AD.gr <- read_tsv("data/AD_gr.txt")
UB.gr <- read_tsv("data/UB_gr.txt")
dir.create("out/upsets/tables/",recursive = T)
dir.create("out/upsets/tables/AD",recursive = T)
dir.create("out/upsets/tables/UB",recursive = T)

get.unique.genes <- function(df, colname) {
  df <- df %>% as_tibble()
  
  f <- rowSums(df %>% filter(across(all_of(colname)) == 1) %>% dplyr::select_if(is.numeric)) == 1
  df.out <- df %>% filter(across(all_of(colname)) == 1)
  return(df.out[f,])
}

for (use.gr in c("AD", "UB")) {
  dir.create(str_glue("./out/upsets/{use.gr}/"))

  if (use.gr == "AD") {
    use <- AD.gr
  }else{
    use <- UB.gr
  }
 
  
  exp.all <- list()
  use.cols <- c()
  
  for (v1 in seq_along(use$col1)) {
    v1.gr <- use[v1, 1]$col1
    v2.gr <- use[v1, 2]$col2
    df.use <- read_tsv(str_glue("out/deseq_tables/{use.gr}/{use.gr}_{v1.gr}_vs_{v2.gr}.txt")) 
    
    use.cols <- append(use.cols, str_glue("{v1.gr}_vs_{v2.gr}"))
    
    df.use <- df.use %>% filter(!stat %in% "None") %>%
      mutate(is_use = 1, group = str_glue("{v1.gr}_vs_{v2.gr}"))
    exp.all[[v1]] <- df.use
  }
  all.df <- do.call("rbind", exp.all) %>% 
    dplyr::select(-log2FoldChange, -padj, -gene_id) %>% 
    spread(key = group, value = is_use )
  all.df[is.na(all.df)] <- 0
  
  write_tsv(all.df, str_glue("out/upsets/tables/{use.gr}_DEGs_table.txt"))
  
  all.df$name <- str_c(all.df$gene, "_", all.df$stat)
  all.df <- all.df %>%  column_to_rownames("name")
  
  p <- upset(data = all.df, point.size = 2,
        sets = use.cols,
        main.bar.color = "#00AEEF",
        queries = list(
          list(
            query = elements, params = list("stat", "up"),color = c("#F26522"),active = T
          )
        )
  ) 
  pdf.name <- str_glue("out/upsets/{use.gr}/{use.gr}_DEGs_upsets.pdf")
  pdf(file = pdf.name,height=5,width=10)
  print (p)
  dev.off()
  # unique genes
  for (v1 in seq_along(use$col1)) {
    v1.gr <- use[v1, 1]$col1
    v2.gr <- use[v1, 2]$col2
    
    uniqu.gene.df <- get.unique.genes(all.df, str_glue("{v1.gr}_vs_{v2.gr}"))[,c("gene", "stat")]
    write_tsv(uniqu.gene.df,str_glue("./out/upsets/tables/{use.gr}/{use.gr}_{v1.gr}_vs_{v2.gr}_DEGs_unique.txt"))
    }
}
