
# Open Data ---------------------------------------------------------------
load(file = file.path(data_dir,"Object_R/AllSample_v3.Robj"))
Sample_list_sc <- gsub("_", ".", Sample_list)

# AllSample_subset <- subset(AllSample, orig.ident %in% c(Sample_list_sc, "FL5.PB1", "FL6.PB1"))

{
  Metadata_All <- AllSample@meta.data %>%
    select(orig.ident, integrated_snn_res.0.3, type, Lymphomatype) %>%
    mutate(Sample_code = gsub("\\.", "_", orig.ident)) %>%
    mutate(cell_id = gsub("\\.", "_", rownames(.))) %>%
    mutate(
      Sample_code = ifelse(Sample_code == "FL5_PB1", "FL5_PB1_3", Sample_code),
      Sample_code = ifelse(Sample_code == "FL6_PB1", "FL6_PB1_2", Sample_code)
    ) %>%
    mutate(
      Sample_code_paper = Sample_code,
      Clusters = integrated_snn_res.0.3
    ) %>%
    mutate(
      Sample_code_paper = recode(Sample_code_paper, !!!Replacement_table),
      Clusters = recode(Clusters, !!!Replacement_table_Sc_all),
      Lymphomatype = factor(Lymphomatype, levels = c("FL","tFL","DLBCL"))
    ) %>%
    mutate(Clusters = factor(Clusters, levels = Replacement_table_Sc_all)) %>%
    mutate(
      Sample_code_paper = factor(Sample_code_paper, levels = Replacement_table),
      Sample_code = factor(Sample_code, levels = Sample_order)
    ) %>%
    mutate(type = if_else(type == "PB", "PBMC", type)) %>%
    mutate(
      Cell_type = "Others",
      Cell_type = (if_else(grepl("CD4", Clusters) == TRUE, "CD4 T-Cells", "Others")),
      Cell_type = (if_else(grepl("CD8", Clusters) == TRUE, "CD8 T-Cells", Cell_type)),
      Cell_type = (if_else(grepl("NK", Clusters) == TRUE, "NK Cells", Cell_type)),
      Cell_type = (if_else(grepl("Bcells", integrated_snn_res.0.3) == TRUE, "B-cells", Cell_type)),
      # Cell_type = (if_else(grepl("Healthy", Clusters) == TRUE, "Normal B Cells", Cell_type)),
      Cell_type = (if_else(grepl("myeloid", integrated_snn_res.0.3) == TRUE, "Monocytes", Cell_type))
    ) %>%
    mutate(Cell_type = factor(Cell_type, levels = c("B-cells", "Monocytes", "CD4 T-Cells", "CD8 T-Cells", "NK Cells", "Others"))) %>%
    left_join(data_talyies_depletion)
  
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$orig.ident, col.name = "orig.ident")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Cell_type, col.name = "Cell_type")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Sample_code_paper, col.name = "Sample_code_paper")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Sample_code, col.name = "Sample_code")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Clusters, col.name = "Clusters")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Response_type_total, col.name = "Response_type_total")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$type, col.name = "Origin")
  AllSample <- AddMetaData(AllSample, metadata = Metadata_All$Lymphomatype, col.name = "Disease")
}

Metadata_All <- as_tibble(AllSample@meta.data)


