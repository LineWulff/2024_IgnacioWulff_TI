library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)        # handy for %>%

#### RNA INTERATION

# Read in each dataset
ha107_lps <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/HA107_LPS_pbmc_multiome.rds")
pbs_lps   <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM PBS LPS 8 8wk/outs/PBS_LPS_pbmc_multiome.rds")
ha107_pbs <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-PBS-8wk/outs/HA107_PBS_pbmc_multiome.rds")
tpbs_pbs  <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-PBS-PBS-8wk/outs/TPBS_PBS_pbmc_multiome.rds")

## Fragment file paths
Fragments(ha107_lps@assays$ATAC) <- NULL
fragments <- CreateFragmentObject(path = "/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/atac_fragments.tsv.gz", cells = colnames(ha107_lps), validate.fragments = TRUE)
Fragments(ha107_lps@assays$ATAC) <- fragments

Fragments(pbs_lps @assays$ATAC) <- NULL
fragments <- CreateFragmentObject(path = "/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM PBS LPS 8 8wk/outs/atac_fragments.tsv.gz", cells = colnames(pbs_lps), validate.fragments = TRUE)
Fragments(pbs_lps @assays$ATAC) <- fragments

Fragments(ha107_pbs@assays$ATAC) <- NULL
fragments <- CreateFragmentObject(path = "/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-PBS-8wk/outs/atac_fragments.tsv.gz", cells = colnames(ha107_pbs), validate.fragments = TRUE)
Fragments(ha107_pbs@assays$ATAC) <- fragments

Fragments(tpbs_pbs@assays$ATAC) <- NULL
fragments <- CreateFragmentObject(path = "/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-PBS-PBS-8wk/outs/atac_fragments.tsv.gz", cells = colnames(tpbs_pbs), validate.fragments = TRUE)
Fragments(tpbs_pbs@assays$ATAC) <- fragments

# Assign orig.ident labels
ha107_lps$orig.ident <- "BM-HA107-LPS-8wk"
pbs_lps$orig.ident   <- "BM-PBS-LPS-8wk"
ha107_pbs$orig.ident <- "BM-HA107-PBS-8wk"
tpbs_pbs$orig.ident  <- "BM-PBS-PBS-8wk"

# ATAC and GEX barcodes seperate
ha107_lps_bar <- read.csv("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-8wk/per_barcode_metrics.csv", row.names = 1)
pbs_lps_bar   <- read.csv("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-8wk/per_barcode_metrics.csv", row.names = 1)
ha107_pbs_bar <- read.csv("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-8wk/per_barcode_metrics.csv", row.names = 1)
tpbs_pbs_bar  <- read.csv("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-8wk/per_barcode_metrics.csv", row.names = 1)

ha107_lps$gex_barcode <- ha107_lps_bar[Cells(ha107_lps),]$gex_barcode
ha107_lps$atac_barcode <- ha107_lps_bar[Cells(ha107_lps),]$atac_barcode
pbs_lps$gex_barcode <- pbs_lps_bar[Cells(pbs_lps),]$gex_barcode
pbs_lps$atac_barcode <- pbs_lps_bar[Cells(pbs_lps),]$atac_barcode
ha107_pbs$gex_barcode <- ha107_pbs_bar[Cells(ha107_pbs),]$gex_barcode
ha107_pbs$atac_barcode <- ha107_pbs_bar[Cells(ha107_pbs),]$atac_barcode
tpbs_pbs$gex_barcode <- tpbs_pbs_bar[Cells(tpbs_pbs),]$gex_barcode
tpbs_pbs$atac_barcode <- tpbs_pbs_bar[Cells(tpbs_pbs),]$atac_barcode

# Ensure unique cell names
ha107_lps <- RenameCells(ha107_lps, new.names = paste(ha107_lps$orig.ident,ha107_lps$atac_barcode,sep="_"))
pbs_lps   <- RenameCells(pbs_lps,   new.names = paste(pbs_lps$orig.ident,pbs_lps$atac_barcode,sep="_"))
ha107_pbs <- RenameCells(ha107_pbs, new.names = paste(ha107_pbs$orig.ident,ha107_pbs$atac_barcode,sep="_"))
tpbs_pbs  <- RenameCells(tpbs_pbs,  new.names = paste(tpbs_pbs$orig.ident,tpbs_pbs$atac_barcode,sep="_"))

# Set RNA as the default assay before integration
DefaultAssay(ha107_lps) <- "RNA"
DefaultAssay(pbs_lps)   <- "RNA"
DefaultAssay(ha107_pbs) <- "RNA"
DefaultAssay(tpbs_pbs)  <- "RNA"

# Apply SCTransform to each object separately
ha107_lps <- SCTransform(ha107_lps, verbose = FALSE)
pbs_lps   <- SCTransform(pbs_lps, verbose = FALSE)
ha107_pbs <- SCTransform(ha107_pbs, verbose = FALSE)
tpbs_pbs  <- SCTransform(tpbs_pbs, verbose = FALSE)

# Create a list of objects for RNA integration
rna_list <- list(ha107_lps, pbs_lps, ha107_pbs, tpbs_pbs)

# Select features for integration
features <- SelectIntegrationFeatures(object.list = rna_list, nfeatures = 3000)

rna_list <- PrepSCTIntegration(object.list = rna_list, anchor.features = features)


rna_anchors <- FindIntegrationAnchors(
  object.list = rna_list, 
  normalization.method = "SCT", 
  anchor.features = features
)
rna_integrated <- IntegrateData(
  anchorset = rna_anchors, 
  normalization.method = "SCT"
)

saveRDS(rna_integrated, "2510_HA107PBSLPSPBS-8wk_RNAdata.rds")

# Run PCA on integrated RNA
rna_integrated <- RunPCA(rna_integrated)

# Run UMAP for visualization
rna_integrated <- RunUMAP(rna_integrated, reduction = "pca", dims = 1:30, reduction.name = "rna_umap", reduction.key = "RNA_UMAP_")

# Visualize RNA Integration
DimPlot(rna_integrated, reduction = "rna_umap", group.by = "dataset") + ggtitle("RNA Integration (Corrected)")


rna_integrated <- FindNeighbors(rna_integrated, dims = 1:30, reduction = "pca")
rna_integrated <- FindClusters(rna_integrated, resolution = 0.4)  # Adjust resolution if needed

## cluster
DimPlot(
  rna_integrated, 
  reduction = "rna_umap", 
  group.by = "seurat_clusters", 
  label = TRUE, 
  repel = TRUE
) + ggtitle("Clusters on RNA Integration (Corrected)")

 # Change features to gene names
FeaturePlot(
  rna_integrated, 
  features = "Il7r",  
  reduction = "rna_umap"
) + ggtitle("Il7r Expression in RNA")
DefaultAssay(rna_integrated) <- "SCT"

FeaturePlot(
  rna_integrated, 
  features = "Siglecf",  
  reduction = "rna_umap",  
  split.by = "dataset"
) + ggtitle("Cd34 Expression in RNA (SCT) - Split by Dataset")


#atac
#Assign dataset labels
ha107_lps$dataset <- "HA107_LPS"
pbs_lps$dataset   <- "PBS_LPS"
ha107_pbs$dataset <- "HA107_PBS"
tpbs_pbs$dataset  <- "TPBS_PBS"

# Ensure unique cell names
ha107_lps <- RenameCells(ha107_lps, add.cell.id = "HA107LPS")
pbs_lps   <- RenameCells(pbs_lps,   add.cell.id = "PBSLPS")
ha107_pbs <- RenameCells(ha107_pbs, add.cell.id = "HA107PBS")
tpbs_pbs  <- RenameCells(tpbs_pbs,  add.cell.id = "TPBSPBS")


DefaultAssay(ha107_lps) <- "ATAC"
DefaultAssay(pbs_lps)   <- "ATAC"
DefaultAssay(ha107_pbs) <- "ATAC"
DefaultAssay(tpbs_pbs)  <- "ATAC"

# Unify peak sets across all datasets
all.peaks <- reduce(x = c(
  granges(ha107_lps),
  granges(pbs_lps),
  granges(ha107_pbs),
  granges(tpbs_pbs)
))

# Rebuild the ATAC assay in each object to use the same peak set
updateATAC <- function(obj) {
  counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features  = all.peaks,
    cells     = colnames(obj)
  )
  obj[["ATAC"]] <- CreateChromatinAssay(counts = counts, fragments = Fragments(obj))
  DefaultAssay(obj) <- "ATAC"
  return(obj)
}

ha107_lps <- updateATAC(ha107_lps)
pbs_lps   <- updateATAC(pbs_lps)
ha107_pbs <- updateATAC(ha107_pbs)
tpbs_pbs  <- updateATAC(tpbs_pbs)

# Preprocessing function for ATAC
prepATAC <- function(obj) {
  obj <- FindTopFeatures(obj, min.cutoff = 10)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj, reduction.name = "lsi", n = 50)
  return(obj)
}

ha107_lps <- prepATAC(ha107_lps)
pbs_lps   <- prepATAC(pbs_lps)
ha107_pbs <- prepATAC(ha107_pbs)
tpbs_pbs  <- prepATAC(tpbs_pbs)
pbmc.combined <- merge(
  x = ha107_lps,
  y = list(pbs_lps, ha107_pbs, tpbs_pbs)
)

# Re-run LSI on the merged object (fixes missing LSI error)
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined, reduction.name = "lsi", n = 50)
peak.features <- rownames(pbmc.combined)  #



pbmc.combined <- RunUMAP(
  pbmc.combined,
  reduction      = "lsi",
  dims           = 2:30,
  reduction.name = "umap_lsi",
  reduction.key  = "uLSI_"
)

DimPlot(pbmc.combined, reduction = "umap_lsi", group.by = "dataset") + ggtitle("Merged (Uncorrected)")
obj.list <- list(ha107_lps, pbs_lps, ha107_pbs, tpbs_pbs)

anchors <- FindIntegrationAnchors(
  object.list     = obj.list,
  anchor.features = peak.features,  
  reduction       = "rlsi",  
  dims            = 2:30
)
pbmc.integrated <- IntegrateEmbeddings(
  anchorset          = anchors,
  reductions         = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate  = 1:30
)

pbmc.integrated <- RunUMAP(
  pbmc.integrated,
  reduction      = "integrated_lsi",
  dims           = 2:30,
  reduction.name = "intUMAP",
  reduction.key  = "intUMAP_"
)

DimPlot(pbmc.integrated, reduction = "intUMAP", group.by = "dataset") + ggtitle("Integrated (Corrected)")

DimPlot(
  pbmc.integrated, 
  reduction = "intUMAP", 
  group.by = "seurat_clusters",
  label = TRUE,  # Show cluster numbers on the plot
  repel = TRUE   # Avoid overlapping labels
) + ggtitle("Clusters on ATAC Integration (Corrected)")


pbmc.integrated <- FindNeighbors(
  pbmc.integrated, 
  reduction = "integrated_lsi", 
  dims = 2:30, 
  graph.name = "integrated_snn"  # Explicitly define graph name
)

pbmc.integrated <- FindClusters(
  pbmc.integrated, 
  graph.name = "integrated_snn",  # Use the correct graph name
  resolution = 0.4
)

DimPlot(
  pbmc.integrated, 
  reduction = "intUMAP", 
  group.by = "seurat_clusters",
  label = TRUE,  
  repel = TRUE   
) + ggtitle("Corrected Seurat Clusters on ATAC Integration")







# 1. Create a named vector of ATAC cluster identities
atac_clusters <- pbmc.integrated$seurat_clusters
names(atac_clusters) <- colnames(pbmc.integrated)

# 2. Find overlapping cell barcodes in your RNA integrated object
common_cells <- intersect(colnames(rna_integrated), names(atac_clusters))

# 3. Make a new metadata column in the RNA object to store ATAC cluster IDs
rna_integrated$ATAC_clusters <- NA
rna_integrated$ATAC_clusters[common_cells] <- atac_clusters[common_cells]

# 4. (Optional) Set identity to that new field
Idents(rna_integrated) <- "ATAC_clusters"








### Transfer ATc to RNA clusters:
# Ensure cell names are consistent between RNA and ATAC
rna_integrated <- RenameCells(rna_integrated, new.names = Cells(pbmc.integrated))


rna_integrated$seurat_clusters <- pbmc.integrated$seurat_clusters
Idents(rna_integrated) <- "seurat_clusters"

DimPlot(
  rna_integrated, 
  reduction = "rna_umap", 
  group.by = "seurat_clusters",
  label = TRUE,  
  repel = TRUE   
) + ggtitle("Seurat Clusters on RNA Integration")


library(patchwork)  # Needed for side-by-side plots

# ATAC Cluster Plot
p1 <- DimPlot(
  pbmc.integrated, 
  reduction = "intUMAP", 
  group.by = "seurat_clusters",
  label = TRUE,  
  repel = TRUE
) + ggtitle("Corrected Seurat Clusters on ATAC")

# RNA Cluster Plot
p2 <- DimPlot(
  rna_integrated, 
  reduction = "rna_umap", 
  group.by = "seurat_clusters",
  label = TRUE,  
  repel = TRUE
) + ggtitle("Corrected Seurat Clusters on RNA")

# Combine plots side by side
p1 | p2



####### UNSUPERVISED CLUSTERING  FOR RNA AND ATAC - CLUSTER 5 AND 8

# Subset clusters 5 & 8
rna_subset_5_8 <- subset(rna_integrated, idents = c("5"))

# Re-run SCTransform on the subset
DefaultAssay(rna_subset_5_8) <- "RNA"
rna_subset_5_8 <- SCTransform(rna_subset_5_8, verbose = FALSE)

# Now run PCA & UMAP
rna_subset_5_8 <- RunPCA(rna_subset_5_8, verbose = FALSE)
rna_subset_5_8 <- RunUMAP(rna_subset_5_8, dims = 1:30)

DimPlot(rna_subset_5_8, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Fresh SCT + PCA/UMAP for clusters 5 & 8")

######### FOR ATAC

# 1) Subset to clusters 5 & 8
Idents(pbmc.integrated) <- "seurat_clusters"
atac_subset_5_8 <- subset(pbmc.integrated, idents = c("5"))

# 2) Set the default assay to "ATAC"
DefaultAssay(atac_subset_5_8) <- "ATAC"

# 3) Re-run the standard Signac steps for this smaller subset
#    (You can adjust min.cutoff or n if you like.)
atac_subset_5_8 <- FindTopFeatures(atac_subset_5_8, min.cutoff = 10)
atac_subset_5_8 <- RunTFIDF(atac_subset_5_8)
atac_subset_5_8 <- RunSVD(atac_subset_5_8, n = 30, reduction.name = "lsi_subset")
# "Cd79a", "Csf3r", "Mpo", "Lyz2", 
#"Gypa", "Adgre1", "Csf1r"

# 4) Run a fresh UMAP on the new 'lsi_subset' reduction
atac_subset_5_8 <- RunUMAP(
  atac_subset_5_8,
  reduction      = "lsi_subset",
  dims           = 2:30,  # or 2:20, etc.
  reduction.name = "umap_subset"
)

# 5) Plot the new subset-based UMAP
DimPlot(
  object   = atac_subset_5_8,
  reduction = "umap_subset",
  group.by  = "seurat_clusters",
  label     = TRUE
) + ggtitle("Fresh UMAP for Clusters 5 & 8 - ATAC Only")


#### Naming Clusters:  (Done fter Feature Plots etc ) and Plotting ATAC and RNA 

# For the RNA integrated object
Idents(rna_integrated) <- "seurat_clusters"
levels(rna_integrated)

Idents(pbmc.integrated) <- "seurat_clusters"
levels(pbmc.integrated)
new.cluster.ids <- c(
  "0" = "Monocytes",
  "1" = "Monocytes",
  "2" = "Neutrophils",
  "3" = "Monocytes",
  "4" = "Monocytes",
  "5" = "LSK",
  "6" = "Monocytes",
  "7" = "ILC",
  "8" = "   ",
  "9" = "NK cells"
  
)
#  active identity is 'seurat_clusters'
Idents(rna_integrated) <- "seurat_clusters"

# Apply RenameIdents() using your new labels
rna_integrated <- RenameIdents(rna_integrated, new.cluster.ids)

# (Optional) Store these new labels in metadata, e.g. 'celltype'
rna_integrated$celltype <- Idents(rna_integrated)

levels(rna_integrated)

Idents(pbmc.integrated) <- "seurat_clusters"
pbmc.integrated <- RenameIdents(pbmc.integrated, new.cluster.ids)
pbmc.integrated$celltype <- Idents(pbmc.integrated)
# RNA UMAP with new labels
p_rna <- DimPlot(
  rna_integrated,
  reduction = "rna_umap",
  group.by = "celltype",  # or just omit if celltype is the Idents
  label = TRUE,
  repel = TRUE
) + ggtitle("RNA UMAP with Named Clusters")

# ATAC UMAP with new labels
p_atac <- DimPlot(
  pbmc.integrated,
  reduction = "intUMAP",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("ATAC UMAP with Named Clusters")

# Side-by-side
library(patchwork)
p_rna | p_atac

## PLOT BY DATSET
# If your cell type names are in 'celltype' metadata (or just use 'seurat_clusters')
DimPlot(
  rna_integrated,
  reduction = "rna_umap",
  group.by = "celltype",      # or "seurat_clusters"
  split.by = "dataset",       # Splits the plot by your "dataset" metadata
  label = TRUE,
  repel = TRUE
) + ggtitle("RNA UMAP, Split by Dataset")

DimPlot(
  pbmc.integrated,
  reduction = "intUMAP",
  group.by = "celltype",      # or "seurat_clusters"
  split.by = "dataset",
  label = TRUE,
  repel = TRUE
) + ggtitle("ATAC UMAP, Split by Dataset")


##### CHROMATIN ACCESSBILITY ANALYSIS - LR 


# Example: LR, with min.pct
da_peaks_pbs <- FindMarkers(
  object          = pbmc.integrated,
  ident.1         = "HA107_LPS",
  ident.2         = "PBS_LPS",
  test.use        = "LR",
  latent.vars     = "nCount_ATAC",
  min.pct         = 0.05,
  logfc.threshold = 0.25
)
da_peaks_pbs$peak <- rownames(da_peaks_pbs)
peak_annot_pbs <- ClosestFeature(pbmc.integrated, da_peaks_pbs$peak)
da_annot_pbs <- merge(da_peaks_pbs, peak_annot_pbs, by.x="peak", by.y="query_region", all.x=TRUE)
write.csv(da_annot_pbs, "HA107LPS_vs_PBSLPS_DA_annot.csv", row.names=FALSE)


###############################################################################
# 0) Setup
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)  # for mouse GO



###############################################################################
# 1) Volcano Plot for HA107_PBS vs TPBS_PBS

res_df <- da_annot_pbs  # rename for convenience

# Define thresholds
p_cutoff  <- 0.05
fc_cutoff <- 0.5  # absolute log2FC threshold 

# Classify Up/Down/NotSig
res_df <- res_df %>%
  mutate(
    logP = -log10(p_val_adj + 1e-300),
    diffexpressed = case_when(
      p_val_adj < p_cutoff & avg_log2FC >=  fc_cutoff  ~ "Up",
      p_val_adj < p_cutoff & avg_log2FC <= -fc_cutoff  ~ "Down",
      TRUE                                            ~ "NotSig"
    )
  )

# Top 10 "Up" & "Down" for labeling
top_up <- res_df %>%
  filter(diffexpressed == "Up") %>%
  arrange(p_val_adj) %>%
  head(10)

top_down <- res_df %>%
  filter(diffexpressed == "Down") %>%
  arrange(p_val_adj) %>%
  head(10)

top_labelled <- bind_rows(top_up, top_down)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = avg_log2FC, y = logP, color = diffexpressed)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "royalblue", "NotSig" = "grey60")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_labelled,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 50
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano: HA107_PBS vs TPBS_PBS",
    x     = "Log2 Fold Change",
    y     = expression(-log[10](Adjusted~p))
  )

print(volcano_plot)
ggsave("Volcano_HA107PBS_vs_TPBSPBS.pdf", plot = volcano_plot, width = 5, height = 5)

###############################################################################
# 2) Extract Up & Down Genes + GO Enrichment
###############################################################################
sig_up <- subset(res_df, p_val_adj < 0.05 & avg_log2FC >  0.5)
sig_down <- subset(res_df, p_val_adj < 0.05 & avg_log2FC < -0.5)

genes_up   <- unique(sig_up$gene_name[!is.na(sig_up$gene_name)])
genes_down <- unique(sig_down$gene_name[!is.na(sig_down$gene_name)])

# GO on Up Genes
ego_up <- enrichGO(
  gene          = genes_up,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# GO on Down Genes
ego_down <- enrichGO(
  gene          = genes_down,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Write to CSV
write.csv(as.data.frame(ego_up),   "GO_Up_HA107PBS_vs_TPBSPBS.csv",   row.names = FALSE)
write.csv(as.data.frame(ego_down), "GO_Down_HA107PBS_vs_TPBSPBS.csv", row.names = FALSE)

# Dotplot
dotplot(ego_up, showCategory = 15) + ggtitle("GO: Up in HA107_PBS")
dotplot(ego_down, showCategory = 15) + ggtitle("GO: Down in HA107_PBS")














####.  FOR CHROMATIN ACCESSIBILTY - CELLTYPE

# 1) Subset to Monocytes, Neutrophils, LSK

Idents(pbmc.integrated) <- "celltype"

sub_obj <- subset(
  pbmc.integrated, 
  idents = c("LSK", "Monocytes", "Neutrohpils")
)

# Confirm the subset:
table(Idents(sub_obj))

# Also ensure identity is set to 'dataset' so we can do HA107_LPS vs PBS_LPS:
Idents(sub_obj) <- "dataset"
table(Idents(sub_obj))


# 2) Run FindMarkers on the subset object

da_markers_lps <- FindMarkers(
  object        = sub_obj,        #
  ident.1       = "HA107_LPS",
  ident.2       = "PBS_LPS",
  assay         = "ATAC",
  test.use      = "LR",
  latent.vars   = "nCount_ATAC",
  min.pct       = 0,              # 
  logfc.threshold = 0             # 
)

cat("\nDone! #Peaks tested:", nrow(da_markers_lps), "\n")

# annotate peaks:
da_markers_lps$peak <- rownames(da_markers_lps)
peak_annot_lps <- ClosestFeature(sub_obj, regions = da_markers_lps$peak)
da_annot_lps <- merge(
  da_markers_lps,
  peak_annot_lps,
  by.x="peak",
  by.y="query_region",
  all.x=TRUE
)

write.csv(da_annot_lps, "DA_HA107LPS_vs_PBSLPS_MonoNeutroLSK_LR_noCutoff.csv", row.names=FALSE)
cat("Saved CSV: DA_HA107LPS_vs_PBSLPS_MonoNeutroLSK_LR_noCutoff.csv\n")


###############################################################################
# 2) Run FindMarkers on the subset object
###############################################################################
da_markers_pps <- FindMarkers(
  object        = sub_obj,        # <--- the subset object
  ident.1       = "HA107_PBS",
  ident.2       = "TPBS_PBS",
  assay         = "ATAC",
  test.use      = "LR",
  latent.vars   = "nCount_ATAC",
  min.pct       = 0,              # no coverage filter
  logfc.threshold = 0             # no fold-change filter
)

cat("\nDone! #Peaks tested:", nrow(da_markers_pps), "\n")

#  annotate peaks:
da_markers_pps$peak <- rownames(da_markers_pps)
peak_annot_pps <- ClosestFeature(sub_obj, regions = da_markers_pps$peak)
da_annot_pps <- merge(
  da_markers_pps,
  peak_annot_pps,
  by.x="peak",
  by.y="query_region",
  all.x=TRUE
)

write.csv(da_annot_pps, "DA_HA107PBS_vs_PBSPBS_MonoNeutroLSK_LR_noCutoff.csv", row.names=FALSE)
cat("Saved CSV: DA_HA107PBS_vs_PBSPBS_MonoNeutroLSK_LR_noCutoff.csv\n")



### Just Monocytes



###############################################################################
# 1) Subset to Monocytes, Neutrophils, LSK
###############################################################################
Idents(pbmc.integrated) <- "celltype"

sub_obj <- subset(
  pbmc.integrated, 
  idents = c("LSK", "Monocytes", "Neutrophils")
)

# Confirm the subset:
table(Idents(sub_obj))

# Also ensure identity is set to 'dataset' so we can do HA107_LPS vs PBS_LPS:
Idents(sub_obj) <- "dataset"
table(Idents(sub_obj))

###############################################################################
# 2) Run FindMarkers on the subset object
###############################################################################
da_markers_lps <- FindMarkers(
  object        = sub_obj,        # <--- the subset object
  ident.1       = "HA107_LPS",
  ident.2       = "PBS_LPS",
  assay         = "ATAC",
  test.use      = "LR",
  latent.vars   = "nCount_ATAC",
  min.pct       = 0,              # no coverage filter
  logfc.threshold = 0             # no fold-change filter
)

cat("\nDone! #Peaks tested:", nrow(da_markers_lps), "\n")

# If desired, annotate peaks:
da_markers_lps$peak <- rownames(da_markers_lps)
peak_annot_lps <- ClosestFeature(sub_obj, regions = da_markers_lps$peak)
da_annot_lps <- merge(
  da_markers_lps,
  peak_annot_lps,
  by.x="peak",
  by.y="query_region",
  all.x=TRUE
)

write.csv(da_annot_lps, "DA_HA107LPS_vs_PBSLPS_MonoNeutroLSK_LR_noCutoff.csv", row.names=FALSE)
cat("Saved CSV: DA_HA107LPS_vs_PBSLPS_MonoNeutroLSK_LR_noCutoff.csv\n")


###############################################################################
# 2) Run FindMarkers on the subset object
###############################################################################
da_markers_pps <- FindMarkers(
  object        = sub_obj,        # <--- the subset object
  ident.1       = "HA107_PBS",
  ident.2       = "TPBS_PBS",
  assay         = "ATAC",
  test.use      = "LR",
  latent.vars   = "nCount_ATAC",
  min.pct       = 0,              # no coverage filter
  logfc.threshold = 0             # no fold-change filter
)

cat("\nDone! #Peaks tested:", nrow(da_markers_pps), "\n")

# If desired, annotate peaks:
da_markers_pps$peak <- rownames(da_markers_pps)
peak_annot_pps <- ClosestFeature(sub_obj, regions = da_markers_pps$peak)
da_annot_pps <- merge(
  da_markers_pps,
  peak_annot_pps,
  by.x="peak",
  by.y="query_region",
  all.x=TRUE
)

write.csv(da_annot_pps, "DA_HA107PBS_vs_PBSPBS_MonoNeutroLSK_LR_noCutoff.csv", row.names=FALSE)
cat("Saved CSV: DA_HA107PBS_vs_PBSPBS_MonoNeutroLSK_LR_noCutoff.csv\n")

