## Load packages
library(tidyverse)
library(Seurat)
library(CARD)
library(SpaTopic)
library(compositions)

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

## Extract SpaTopic results
for (st in 1:length(selected_chip)) {
    st_obj <- readRDS(paste0(result_dir, selected_chip[st], ".rds"))
    ct_topic_data <- Misc(object = st_obj, slot = "celltype_topic") %>% as.data.frame()
    ct_topic_data <- column_to_rownames(ct_topic_data, var = "celltype")
    colnames(ct_topic_data) <- paste0(selected_chip[st], "_", colnames(ct_topic_data))
    ct_topic_data <- rownames_to_column(ct_topic_data, var = "celltype")
    if (st == 1) {
        CellTopic_df <- ct_topic_data
    } else {
        CellTopic_df <- merge(
            CellTopic_df, ct_topic_data,
            by = "celltype", all = TRUE
        )
    }
}
CellTopic_df[is.na(CellTopic_df)] <- 0
CellTopic_df <- CellTopic_df %>% column_to_rownames(var = "celltype")

## Remove outlier
CellTopic_df <- CellTopic_df[, !grepl("P6_YF_3_CellTopic5", colnames(CellTopic_df))]

## Perform centered log-ratio transformation
CellTopic_df_clr <- clr(t(as.matrix(CellTopic_df))) %>% t() %>% as.data.frame()

## Calculate correlation
corr <- cor(CellTopic_df_clr)

## Hierarchical clustering 
cluster_num <- 12
p <- pheatmap::pheatmap(corr, 
    border_color = "black", 
    scale = "none", 
    cluster_rows = T, 
    cluster_cols = T, 
    legend = TRUE, 
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("low", " ", "high"), 
    show_rownames = T, 
    show_colnames = T, 
    fontsize = 8,
    cutree_rows = cluster_num,
    cutree_cols = cluster_num,
    silent = TRUE
)

## Obtain MetaTopic assignment
row_cluster <- cutree(p$tree_row, k = cluster_num) %>% as.data.frame()

## Save MetaTopic assignment
write.table(row_cluster, paste0(result_dir, "row_cluster.txt"), sep = "\t")