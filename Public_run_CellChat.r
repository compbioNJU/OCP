## Load packages
library(tidyverse)
library(Seurat)
library(CellChat)

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

## Run CellChat
options(future.globals.maxSize= 4*1024*1024^2)
for (st in 1:length(selected_chip)) {
    # Load data
    st_obj <- readRDS(paste0(result_dir, selected_chip[st], ".rds"))
    # Prepare CellChat input
    Idents(st_obj) <- st_obj$MetaTopic
    data.input <- Seurat::GetAssayData(st_obj, layer = "data", assay = "SCT")
    meta <- data.frame(labels = Idents(st_obj), row.names = names(Idents(st_obj)))
    spatial.locs <- Seurat::GetTissueCoordinates(st_obj, scale = NULL, cols = c("imagerow", "imagecol"))
    spatial.factors <- data.frame(ratio = 10, tol = 50)
    # Create CellChat object
    cellchat <- createCellChat(data.input, meta = meta, datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
    # Select ligand receptor database
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    # Run CellChat
    future::plan("multisession", workers = 8)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat,
        type = "truncatedMean", trim = 0.1,
        distance.use = TRUE, interaction.range = 400, scale.distance = 0.01,
        contact.dependent = TRUE, contact.range = 200
    )
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    cellchat <- aggregateNet(cellchat)
    # Save CellChat result
    saveRDS(cellchat, paste0(result_dir, selected_chip[st], "_Meta_cellchat.rds"))
}