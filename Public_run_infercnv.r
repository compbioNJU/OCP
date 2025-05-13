## Load packages
library(infercnv)
library(Seurat)

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
gene_order_file_path <- './gene_pos_symbol.txt'

## Run infercnv
for (st in 1:length(selected_chip)) {
    # Load data
    sc_obj_patient <- readRDS(paste0(result_dir, selected_chip[st], ".rds"))
    ref_list <- readRDS(paste0(result_dir, "infer_cnv_ref_list.rds"))

    # Set reference
    chosen_ref <- ref_list[[st]]

    # Prepare input
    counts <- sc_obj_patient@assays$Spatial@counts
    anno <- sc_obj_patient[["CellTopic"]]

    # Create infercnv object
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts,
                                            annotations_file = anno,
                                            delim="\t",
                                            gene_order_file = gene_order_file_path,   
                                            ref_group_names = chosen_ref
                                            )

    # Run infercnv
    infercnv_obj <- infercnv::run(infercnv_obj,
                                    cutoff=0.1,  
                                    out_dir=paste0(result_dir,selected_chip[st],'/'),
                                    cluster_by_groups=T, 
                                    denoise=T,
                                    HMM=T,no_plot = T,analysis_mode = 'samples')
}


