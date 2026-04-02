options(future.globals.maxSize = 10 * 1024^3)

# Open Data ---------------------------------------------------------------
library(here)
library(gtools)
library(Seurat)
library(tidyverse)
source(here('src/00_Open_metadata.R'))
source(here('src/01_Open_sc_files.R'))

#Filter Monocytes et PBMC
AllSample_subset_monocytes <- subset(AllSample, Cell_type %in% "Monocytes" & Origin == "PBMC")
rm(AllSample)
FeaturePlot(AllSample_subset_monocytes, features =(c("CD4","CD3E","CD14")))
length(rownames(AllSample_subset_monocytes@meta.data))

# Define bounds
x_min <- 0
y_min <- 7.5
x_max <- 10
y_max <- 20

# Extract UMAP coordinates from the reduction object
umap_coords <- AllSample_subset_monocytes[["umap_all"]]@cell.embeddings # Or use @reductions$umap@data if using older versions

# Create logical vector
keep_cells <- (umap_coords[, 1] >= x_min) & 
  (umap_coords[, 1] <= x_max) & 
  (umap_coords[, 2] >= y_min) & 
  (umap_coords[, 2] <= y_max)

# Subset the Seurat object
AllSample_subset_monocytes_filtered <- subset(
  AllSample_subset_monocytes, 
  cells = rownames(umap_coords)[keep_cells]  
)
length(rownames(AllSample_subset_monocytes_filtered@meta.data))


AllSample_subset_monocytes_filtered <- subset(AllSample_subset_monocytes_filtered, CD3E<0.35 & MS4A1<0.27 & CD19<0.27 & CD22<0.27 & CD3D < 0.27)
length(rownames(AllSample_subset_monocytes_filtered@meta.data))
FeaturePlot(AllSample_subset_monocytes_filtered, features =(c("CD4","CD3E","CD14",
                                                              "CD3D","CD8A","CD48","CD38","BCL6","CD19", "MS4A1", "CD22", "ITGAM", "FCGR3A", "CD3E", "CD4", "CD8A", "NCAM1","IL3RA")))


pochi::MetaDataPlot(AllSample_subset_monocytes_filtered, group.by = "Cell_type", split.by = "Sample_code", as.freq = F) +
  labs(
    title = "",
    x = "", # Éviter de répéter des informations évidentes
    y = "Frequency",
    fill = "Cell Type"
  ) +
  labs(fill = "Cell Type") +
  scale_fill_brewer(palette = "GnBu", direction = -1) +
  theme_blood() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold")) # Inclinaison des labels



# Data Integration Seurat V5 ----------------------------------------------
# AllSample_subset_monocytes_filtered <- LoadSeuratRds(here("data/Seurat_obj/AllSample_subset_monocytes_filtered.rds"))
AllSample_subset_monocytes_filtered <- LoadSeuratRds(here("data/Seurat_obj/AllSample_subset_monocytes_cca_integrated.rds"))

# AllSample_subset_monocytes_filtered <- subset(AllSample_subset_monocytes_filtered, cells = cells_to_remove, invert = TRUE)

AllSample_subset_monocytes_filtered[["percent.mt"]] <- PercentageFeatureSet(AllSample_subset_monocytes_filtered, pattern = "^MT-")

VlnPlot(AllSample_subset_monocytes_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AllSample_subset_monocytes_filtered <- subset(AllSample_subset_monocytes_filtered, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 7)

plot1 <- FeatureScatter(AllSample_subset_monocytes_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AllSample_subset_monocytes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

cells_to_remove <- WhichCells(AllSample_subset_monocytes_filtered, idents = c(5))
AllSample_subset_monocytes_filtered <- subset(AllSample_subset_monocytes_filtered, cells = cells_to_remove, invert = TRUE)


AllSample_subset_monocytes_filtered[["RNA"]] <- split(AllSample_subset_monocytes_filtered[["RNA"]], f = AllSample_subset_monocytes_filtered$Sample_code)
AllSample_subset_monocytes_filtered <- NormalizeData(AllSample_subset_monocytes_filtered)
AllSample_subset_monocytes_filtered <- FindVariableFeatures(AllSample_subset_monocytes_filtered)
AllSample_subset_monocytes_filtered <- ScaleData(AllSample_subset_monocytes_filtered)
AllSample_subset_monocytes_filtered <- RunPCA(AllSample_subset_monocytes_filtered)

AllSample_subset_monocytes_filtered <- FindNeighbors(AllSample_subset_monocytes_filtered, dims = 1:30, reduction = "pca")
AllSample_subset_monocytes_filtered <- FindClusters(AllSample_subset_monocytes_filtered, resolution = 2, cluster.name = "unintegrated_clusters")

AllSample_subset_monocytes_filtered <- RunUMAP(AllSample_subset_monocytes_filtered, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(AllSample_subset_monocytes_filtered, reduction = "umap.unintegrated", group.by = c("Disease","Sample_code"))

AllSample_subset_monocytes_filtered <- IntegrateLayers(
  object = AllSample_subset_monocytes_filtered, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

AllSample_subset_monocytes_filtered <- IntegrateLayers(
  object = AllSample_subset_monocytes_filtered, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  k.weight = 80,
  verbose = FALSE
)

AllSample_subset_monocytes_filtered <- FindNeighbors(AllSample_subset_monocytes_filtered, reduction = "integrated.cca", dims = 1:30)
AllSample_subset_monocytes_filtered <- FindNeighbors(AllSample_subset_monocytes_filtered, reduction = "integrated.rpca", dims = 1:30)


AllSample_subset_monocytes_filtered <- FindClusters(AllSample_subset_monocytes_filtered, resolution = 0.5, cluster.name = "cca_clusters")
AllSample_subset_monocytes_filtered <- FindClusters(AllSample_subset_monocytes_filtered, resolution = 0.5, cluster.name = "rpca_clusters")

AllSample_subset_monocytes_filtered <- RunUMAP(AllSample_subset_monocytes_filtered, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
AllSample_subset_monocytes_filtered <- RunUMAP(AllSample_subset_monocytes_filtered, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

DimPlot(
  AllSample_subset_monocytes_filtered,
  reduction = "umap.cca",
  group.by = c("Disease", "Sample_code", "cca_clusters"),
  combine = , label.size = 2
)

DimPlot(
  AllSample_subset_monocytes_filtered,
  reduction = "umap.rpca",
  group.by = c("Disease", "Sample_code", "rpca_clusters"),
  combine = , label.size = 2
)

DimPlot(
  AllSample_subset_monocytes_filtered,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  combine = , label.size = 2
)

DimPlot(
  AllSample_subset_monocytes_filtered,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters"),
  combine = , label.size = 5,
  label = TRUE
)

AllSample_subset_monocytes_filtered <- AllSample_subset_monocytes

Idents(object = AllSample_subset_monocytes_filtered) <- "rpca_clusters"

AllSample_subset_monocytes_filtered <- JoinLayers(AllSample_subset_monocytes_filtered)
AllSample_subset_monocytes_filtered.markers <- FindAllMarkers(AllSample_subset_monocytes_filtered, logfc.threshold = 0.25, min.pct = 0.25)

AllSample_subset_monocytes_filtered.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 17) %>%
  ungroup() -> top10

rio::export(top10, "test.tsv")

# 2. Rename your clusters for the plot
new_ident_names <- c(
  "0" = "Resident Macrophages",
  "1" = "Activated Classical Monocytes",
  "2" = "Mature Classical Monocytes",
  "3" = "CD16+ Non-classical Monocytes",
  "4" = "Early Classical Monocytes",
  "5" = "CD16+ transitional Monocytes",
  "6" = "cDC2",
  "7" = "IFN-stimulated Monocytes",
  "8" = "CD11c+ Non-classical Monocytes",
  "9" = "Proliferating Monocytes",
  "10" = "cDC1"
)

# Apply names to the object
Idents(object = AllSample_subset_monocytes_filtered) <- "rpca_clusters"
AllSample_subset_monocytes_filtered <- RenameIdents(AllSample_subset_monocytes_filtered, new_ident_names)

# Define your desired biological order
biological_order <- c(
  "Early Classical Monocytes", 
  "Proliferating Monocytes",
  "Mature Classical Monocytes", 
  "Activated Classical Monocytes", 
  "IFN-stimulated Monocytes",
  "Resident Macrophages", 
  "CD16+ transitional Monocytes",
  "CD16+ Non-classical Monocytes",
  "CD11c+ Non-classical Monocytes",
  "cDC1",
  "cDC2"

)

# Apply the order to the object's active identities
levels(AllSample_subset_monocytes_filtered) <- rev(biological_order)


# 1. Define the canonical markers based on your DE results
myeloid_markers <- c(
  # Monocyte Lineage
  "CD14","CCR2", "VCAN", "S100A9", "S100A8","ANXA1",  # Cluster 4: Early Classical
  "F5",            
  "PADI4", "ALOX5AP", "RETN", "SLC2A3",      # Cluster 2: Mature Classical

  "IL1B", "CXCL8", "NLRP3", "TREM1",         # Cluster 1: Inflammatory Mac/Mono
  
  "ISG15", "IFIT1", "MX1", "SIGLEC1",        # Cluster 7: ISG-High (Interferon)
  


  # Macrophage

  "VSIG4", "CD9", "SMAD7", "OLR1",           # Cluster 0: Resident Mac
  "C1QA", "C1QB", "LIPA",            # Cluster 5: C1q+ Scavenging Mac
  
  # Non Classical
  "CDKN1C", "FCGR3A", "HES4", "RHOC",        # Cluster 3: Non-classical
  "SAT1", "HSPA6", "ITGAX", "LILRB2",        # Cluster 8: Activated Non-classical
  
  
  # Dendritic Cells
  "CLEC9A","C1orf54", "CADM1", "XCR1",       # Cluster 10: cDC1 (Villani et al. 2019)
  "CD1C", "FCER1A", "CLEC10A",               # Cluster 6: cDC2 (Villani et al. 2019)

  # State
  "MKI67", "TOP2A", "BIRC5", "CENPF"     # Cluster 9: Cycling Myeloid
  
  
)

# 3. Generate the DotPlot
DotPlot(AllSample_subset_monocytes_filtered, features = myeloid_markers, dot.scale = 8) + 
  RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Myeloid Subset Canonical Markers", x = "Genes", y = "Cluster")

DimPlot(
  AllSample_subset_monocytes_filtered,
  reduction = "umap.rpca",
  combine = , label.size = 5,
  label = TRUE
)
VlnPlot(AllSample_subset_monocytes_filtered, features = "LGALS9",pt.size = 0,split.by = "Disease")

VlnPlot(AllSample_subset_monocytes_filtered, features = "LGALS9",pt.size = 0)

pochi::MetaDataPlot(AllSample_subset_monocytes_filtered, group.by = "Idents", split.by = "Disease", as.freq = T) +
  labs(
    title = "",
    x = "", # Éviter de répéter des informations évidentes
    y = "Frequency",
    fill = "Cell Type"
  ) +
  labs(fill = "Cell Type") +
  # scale_fill_brewer(palette = "GnBu", direction = -1) +
  # theme_blood() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold")) # Inclinaison des labels


# # Data Integration with Anchors -------------------------------------------
# split the object by dataset
samples.list <- SplitObject(AllSample_subset_monocytes_filtered, split.by = "Sample_code")

  # perform standard preprocessing on each object
  for (i in 1:length(samples.list)) {
    samples.list[[i]] <- NormalizeData(samples.list[[i]], verbose = FALSE)
    samples.list[[i]] <- FindVariableFeatures(
      samples.list[[i]], selection.method = "vst",
      nfeatures = 3000, verbose = FALSE
    )
  }

  # find anchors
  anchors <- FindIntegrationAnchors(object.list = samples.list)

  # integrate data
  integrated <- IntegrateData(anchorset = anchors)


# # Find Clusters + Umap ----------------------------------------------------
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
ElbowPlot(integrated)

integrated <- FindNeighbors(integrated, dims = 1:11)
integrated <- FindClusters(integrated, resolution = 0.9)
integrated <- RunUMAP(integrated, dims = 1:11)
DimPlot(integrated, reduction = "umap")

FeaturePlot(integrated, features =(c("CD4","LGALS9","CD14",
                                                              "CD3D","CD8A","CD48","CD38","BCL6","CD19", "MS4A1", "CD22", "ITGAM", "FCGR3A", "CD3E", "CD4", "CD8A", "NCAM1","IL3RA")))


integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE)

integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# library(rio)
# export(top10, "test.tsv")

DoHeatmap(integrated, features = top10$gene) + NoLegend()
SaveSeuratRds(integrated, file = here("data/Seurat_obj/AllSample_subset_monocytes_integrated.rds"))




# SaveSeuratRds(AllSample_subset_monocytes_filtered, file = here("data/Seurat_obj/AllSample_subset_monocytes_cleaned_cca_rpca_integrated.rds"))
