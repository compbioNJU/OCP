## Load packages
library(tidyverse)
library(Seurat)
library(CARD)
library(SpaTopic)

## Selected chip labels
selected_chip <- c(
    "P0_FT_1", "P0_OC_2", "P2_OM_z", "P2_FT_2", "P2_ZF_1",
    "P2_YF_1", "P3_PM_z", "P3_ZY_3", "P4_IC_z", "P4_OM_2",
    "P4_AW_z", "P4_CO_2", "P5_ME_z", "P5_OM_2", "P5_FT_1",
    "P5_YF_1", "P5_ZF_1", "P6_OM_2", "P6_LI_1", "P6_FT_1",
    "P6_SOM_", "P6_YF_3", "P6_ZF_1"
)

## File paths of deconvolution results and bin100 data
data_dir <- "./data/"
CARD_name <- list.files(data_dir, pattern = "CARDobj.rds")
CARD_name_minor <- CARD_name[!grepl("major", CARD_name)]
data_name <- list.files(data_dir, pattern = "bin100.rds")

## Only choose selected chips
CARD_name_minor <- CARD_name_minor[grepl(paste(selected_chip, collapse = "|"), CARD_name_minor)]
data_name <- data_name[grepl(paste(selected_chip, collapse = "|"), data_name)]

## Match chip and data
CARD_name_minor <- CARD_name_minor[match(selected_chip, gsub(".CARDobj.rds", "", CARD_name_minor))]
data_name <- data_name[match(selected_chip, gsub("_bin100.rds", "", data_name))]

## Run SpaTopic
result_dir <- "./result/"
cluster_resolution =  readRDS('./cluster_resolution.rds')

for (st in 1:length(selected_chip)) {

    # Load deconvolution result
    CARD_obj <- readRDS(paste0(data_dir, CARD_name_minor[st]))
    spot_celltype <- CARD_obj@Proportion_CARD

    # Load and process data
    st_obj <- readRDS(paste0(data_dir, data_name[st]))
    st_obj <- RunPCA(st_obj, assay = "SCT", verbose = FALSE)
    st_obj <- FindNeighbors(st_obj, reduction = "pca", dims = 1:30)
    st_obj <- FindClusters(st_obj, verbose = FALSE, resolution = cluster_resolution[st])

    # SpaTopic
    spot_clusters <- st_obj@meta.data
    SpaTopic_res <- SpaTopic::CellTopic(
        spot_celltype,
        spot_clusters,
        cluster = "seurat_clusters",
        num_topics = 15,
        percent = 0.7)
    CellT <- SpaTopic_res[["CellTopic"]]["CellTopic"]
    ct_topic_data <- SpaTopic_res[["celltype_topic"]]
    Misc(object = st_obj, slot = "celltype_topic") <- ct_topic_data %>% rownames_to_column(var = "celltype")
    st_obj <- AddMetaData(st_obj, CellT)

    # Save seurat object
    saveRDS(st_obj, paste0(result_dir, selected_chip[st], ".rds"))
}