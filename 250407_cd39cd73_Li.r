BiocManager::install("GEOquery")
BiocManager::install("limma")
pacman::p_load(GEOquery, tidyverse, ggrepel, survminer, limma, oligo, DT, marray)
pacman::p_load(pheatmap, tidyplots)
pacman::p_load(data.table, biomaRt, arrayQualityMetrics, testit, gplots, AgiMicroRna)
library(reshape2)

# Load the data
expr_data <- read.delim("/Users/denriquez/Downloads/GSE208679_processed_data.txt", 
                       header = TRUE, 
                       check.names = FALSE) 
anno_data <- read.delim("/Users/denriquez/Downloads/Human.GRCh38.p13.annot.tsv", 
                       header = TRUE, 
                       check.names = FALSE) 

expr_data_norm = read.delim("/Users/denriquez/Downloads/GSE208679_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", 
                       header = TRUE, 
                       check.names = FALSE) 


anno_data |>
    filter(if_any(1:ncol(anno_data), ~grepl("NT5E|ENTPD1", ., ignore.case = TRUE)))                 

raji = c("GSM6364108", "GSM6364109", "GSM6364110")

raji_exp = expr_data_norm |>
    select(GeneID, contains(raji)) |>
    filter(GeneID%in%c("4907", "953")) |>
    mutate(GeneID=dplyr::recode(GeneID, "4907" = "NT5E", "953" = "ENTPD1")) |>
    mutate(GeneID = factor(GeneID, levels = c("NT5E", "ENTPD1"))) |>
    pivot_longer(-GeneID, names_to = "Sample", values_to = "FPKM") |>
    mutate(cell="Raji")

expr_data2 <- read.delim("/Users/denriquez/Downloads/GSE108601_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", 
                       header = TRUE, 
                       check.names = FALSE) 
head(expr_data2) 

mt1 =c("GSM2905805", "GSM2905806")
mt1_exp = expr_data2 |>
    select(GeneID, contains(mt1)) |>
    filter(GeneID%in%c("4907", "953")) |>
    mutate(GeneID=dplyr::recode(GeneID, "4907" = "NT5E", "953" = "ENTPD1")) |>
    mutate(GeneID = factor(GeneID, levels = c("NT5E", "ENTPD1"))) |>
    pivot_longer(-GeneID, names_to = "Sample", values_to = "FPKM") |>
    mutate(cell="MT1") 

datasetx = rbind(raji_exp, mt1_exp)
datasetx |>
    tidyplots::tidyplot(x=GeneID, y=FPKM, color=cell) |>
    tidyplots::add_median_bar() |>
    tidyplots::add_data_points_beeswarm() |>
    tidyplots::adjust_x_axis(labels=c("NT5E"="CD73", "ENTPD1"="CD39")) |>
    tidyplots::adjust_y_axis_title("Gene expression (FPKM)") |>
    tidyplots::adjust_x_axis_title("") |>
    tidyplots::adjust_legend_title("") |>
    tidyplots::add_test_asterisks() |>
    tidyplots::remove_caption() |>
    tidyplots::save_plot("/Users/denriquez/Downloads/raji_mt1.png",  bg="transparent")




