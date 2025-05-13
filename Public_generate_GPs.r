## Load packages
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ComplexHeatmap)

## Selected chip labels
selected_chip <- c(
    "P0_FT_1", "P0_OC_2", "P2_OM_z", "P2_FT_2", "P2_ZF_1",
    "P2_YF_1", "P3_PM_z", "P3_ZY_3", "P4_IC_z", "P4_OM_2",
    "P4_AW_z", "P4_CO_2", "P5_ME_z", "P5_OM_2", "P5_FT_1",
    "P5_YF_1", "P5_ZF_1", "P6_OM_2", "P6_LI_1", "P6_FT_1",
    "P6_SOM_", "P6_YF_3", "P6_ZF_1"
)

## Set path
result_dir <- "./result/"

## Find DEGs for each CellTopic
Cell_markers_all <- data.frame()
for(st in 1:length(selected_chip)){
    st_obj <- readRDS(paste0(result_dir, selected_chip[st], ".rds"))
    Idents(st_obj) <- "CellTopic"
    Cell_markers <- FindAllMarkers(st_obj)
    Cell_markers <- subset(Cell_markers, p_val < 0.05)
    Cell_markers$chip <- selected_chip[st]
    rownames(Cell_markers) <- NULL
    Cell_markers_all <- rbind(Cell_markers_all, Cell_markers)
}

## Filter DEGs
Cell_markers_all_select <- Cell_markers_all %>% subset(p_val_adj <= 0.0001 & avg_log2FC > 0.25)

## Load MetaTopic assignment
row_cluster <- read.table(paste0(result_dir, "row_cluster.txt"), sep = "\t")
row_cluster$chip <- gsub("_CellTopic[0-9]*", "", rownames(row_cluster))

## Binarize the result of DEG analysis
gene_table <- table(Cell_markers_all_select$gene, paste0(Cell_markers_all_select$chip, "_", Cell_markers_all_select$cluster))
row_cluster_a <- row_cluster[rownames(row_cluster) %in% colnames(gene_table), ]
gene_table <- gene_table[, rownames(row_cluster_a)[order(row_cluster_a$cluster)]]

## Calculate the frequency of each gene in each MetaTopic
group <- sort(row_cluster_a$cluster)
gene_count_list <- gene_table[, ] %>% apply(1, function(x, group) {
    df <- data.frame(count = x, group = as.character(group))
    df <- df %>%
        group_by(group) %>%
        summarise(sum = sum(count) / n())
    df <- df %>% column_to_rownames(var = "group")
    return(df)
}, group)
gene_count <- Reduce(cbind, gene_count_list)
colnames(gene_count) <- rownames(gene_table)

## Choose top 30 genes
gene <- gene_count[, ] %>% apply(1, function(x) {
    names(x)[(rank(-x, ties.method = "min") <= 30) & (x > 0)]
})

## Filter genes without corresponding ENTREZID
gene <- gene %>% lapply(function(x) {
    a <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% as.data.frame()
    a$SYMBOL
})

## Filter ribosomal genes and mitochondrial genes
gene <- gene %>% unlist() %>% unique()
gene <- gene[!grepl("^RP[SL]", gene)]
gene <- gene[!grepl("^MT-", gene)]

## Obtain the average expression of filtered genes
for (st in 1:23) {
    # Load and process data
    st_obj <- readRDS(paste0(result_dir, selected_chip[st], ".rds"))
    data <- GetAssayData(object = st_obj, layer = "data")
    data <- data %>% t() %>% scale() %>% t()
    # Average the expression levels
    for (i in unique(st_obj$CellTopic)) {
        df <- data[, st_obj$CellTopic == i] %>% apply(1, function(x) { mean(x) }) %>% as.data.frame()
        colnames(df) <- paste0(selected_chip[st], "_", i)
        df <- df %>% rownames_to_column(var = "gene")
        if (i == unique(st_obj$CellTopic)[1]) {
            data_df_all <- df
        } else {
            data_df_all <- merge(data_df_all, df, by = "gene", all = T)
        }
    }
    # Summarize the results
    if (st == 1) {
        data_df_all_chip <- data_df_all
    } else {
        data_df_all_chip <- merge(data_df_all_chip, data_df_all, by = "gene", all = T)
    }
}

## Process matrix for hierarchical clustering 
data_df_all_chip[is.na(data_df_all_chip)] <- -2.58
data_df_all_chip <- column_to_rownames(data_df_all_chip, var = "gene")
data_df_all_chip <- data_df_all_chip[, rownames(row_cluster)]
data_df_all_chip_scale <- data_df_all_chip %>% t() %>% scale() %>% t()
data_df_all_chip_scale[data_df_all_chip_scale > 2] <- 2
data_df_all_chip_scale[data_df_all_chip_scale < -2] <- -2

## Hierarchical clustering 
set.seed(123)
cluster_num_row = 5
cluster_num_column = 8
p <- Heatmap(as.matrix(data_df_all_chip_scale[gene, ]),
    col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    show_row_names = F,
    show_column_names = F,
    column_title = NULL,
    cluster_rows = T,
    cluster_columns = T,
    row_km = cluster_num_row,
    column_km = cluster_num_column
)
row_order_list <- row_order(p)

## Obtain GP1-GP5
heatmap_mtx <- as.matrix(data_df_all_chip_scale[gene, ])
gene_list <- row_order_list %>% lapply(function(x, heatmap_mtx) {
    rownames(heatmap_mtx)[x]
}, heatmap_mtx)
new_row_order <- c(4, 5, 1, 2, 3)
GPs <- gene_list[new_row_order]
names(GPs) <- paste0("GP", 1:5)