library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(UCell)

dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
## colouring
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
mycols<- myColorRamp(mycols, seq(1:50))
## Project and sample info
proj_data_dir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/"
## output plots dir
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/output/"

col_cols <- c("PBS"="#F8766D","HA107"="#00BFC4")
stim_cols <- c("PBS"="white","LPS"="lightgrey")


#### ----obj read in plus set levels ---- ####
obj@meta.data$orig.ident <- factor(obj$orig.ident, levels =
                                     c("BM-PBS-PBS-21d","BM-PBS-LPS-21d","BM-HA107-PBS-21d","BM-HA107-LPS-21d",
                                       "BM-PBS-PBS-8wk","BM-PBS-LPS-8wk","BM-HA107-PBS-8wk","BM-HA107-LPS-8wk" ))
obj@meta.data$ID_labs_ext <- factor(obj$ID_labs_ext, levels =
                                      c("Q1","Q2","Q3","Ly6c lo monocytes","Ly6c hi monocytes","LSK", "Neutrophils","Dendritic cells", "NK cells"))
Idents(obj) <- 'ID_labs_ext'

obj_sub <- subset(obj, idents = c("Q1","Q2","Q3"))

# violin colors
obj_sub$colonization <- factor(obj_sub$colonization, levels=c("PBS","HA107"))


#### ---- Hallmark signatures ---- ####
# ── 1. Get Hallmark gene sets for mouse ───────────────────────────────────────
options(timeout = 600)
hallmarks <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

C2 <- msigdbr(species = "Mus musculus", category = "C2") %>% # CGP
  dplyr::select(gs_name, gene_symbol)

C7 <- msigdbr(species = "Mus musculus", category = "C7", subcollection = "IMMUNESIGDB") %>% # CGP
  dplyr::select(gs_name, gene_symbol)


grep("CUI", unique(C7$gs_name), ignore.case = TRUE, value = TRUE)

# ── 2. Extract your two gene sets of interest ─────────────────────────────────
ifn_genes <- hallmarks %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene_symbol)

inflam_genes <- hallmarks %>%
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
  pull(gene_symbol)


# ── 1. Get Hallmark gene sets ─────────────────────────────────────────────────
ifn_genes <- hallmarks %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene_symbol)

ifng_genes <- hallmarks %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(gene_symbol)

inflam_genes <- hallmarks %>%
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
  pull(gene_symbol)

IFNaup_genes <- C2 %>%
  filter(gs_name == "DER_IFN_GAMMA_RESPONSE_UP") %>%
  pull(gene_symbol)

MAPK_genes <- C2 %>%
  filter(gs_name == "KEGG_MAPK_SIGNALING_PATHWAY") %>%
  pull(gene_symbol)

NemLPS_genes <- C2 %>%
  filter(gs_name == "NEMETH_INFLAMMATORY_RESPONSE_LPS_UP") %>%
  pull(gene_symbol)


# ── 2. Score cells in your Seurat object ──────────────────────────────────────
# Replace "condition" with your actual metadata column name
# Replace "IL10" and "LPS" with your actual condition labels

signatures <- list(
  IFN_Alpha     = ifn_genes,
  Inflammatory  = inflam_genes,
  ISG = IFNaup_genes,
  IFN_Gamma = ifng_genes,
  MAPK = MAPK_genes,
  NemLPS = NemLPS_genes)

# UCell is preferred for scRNAseq — rank-based, robust to library size
obj <- AddModuleScore_UCell(obj, features = signatures, assay = "imputed_t2")
obj <- AddModuleScore(obj, features = list(ifn_genes), name = "IFN_Alpha", ctrl = 100, assay = "imputed_t2" )
obj <- AddModuleScore(obj, features = list(inflam_genes), name = "Inflammatory", ctrl = length(inflam_genes), assay = "imputed_t2" )
obj <- AddModuleScore(obj, features = list(IFNaup_genes), name = "ISG", ctrl = length(IFNaup_genes), assay = "imputed_t2" )


VlnPlot(obj, features = c("IFN_Alpha_UCell","IFN_Alpha1"), pt.size = 0, group.by = "orig.ident")

obj_sub <- subset(obj, idents = c("Q1","Q2","Q3"))


# ── 3. Extract scores + condition metadata ────────────────────────────────────
sign <- "IFN_Alpha_UCell" # "IFN_Alpha_UCell" "Inflammatory_UCell" "Inflammatory1" "IFN_Alpha1"
score_df <- obj_sub@meta.data %>%
  dplyr::select(colonization,stimulation, ID_labs_ext, orig.ident,         # your condition column
                sign) %>%
  pivot_longer(cols = sign,
               names_to  = "signature",
               values_to = "score") %>%
  mutate(signature = str_remove(signature, "_UCell"),
         signature = str_replace(signature, "_", " "))
sign <- sign %>% str_remove("_UCell") %>%
                        str_replace("_", " ") %>% str_remove("1")

summary_stats <- score_df %>%
  group_by(orig.ident,colonization,stimulation, ID_labs_ext) %>%
  summarise(
    mean_val = mean(score),
    median_val = median(score)
  )

# ── 4. Visualize — violin plot per condition and signature ────────────────────
pdf(paste0(outdir,dato,"_MSigDB_",sign,"_imp2.pdf"),height = 4, width = 7)
ggplot(score_df, aes(x = orig.ident, y = score, fill = stimulation, color=colonization)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  #geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  facet_grid(~ID_labs_ext) +
  scale_fill_manual(values = stim_cols) +
  scale_color_manual(values = col_cols)+
  geom_point(data = summary_stats, aes(y = median_val), color = "black", size = 2) + # points at mean
  geom_line(data = summary_stats, aes(y = median_val, group = colonization), color = "black", linewidth = 0.5) + # line connecting means
  labs(title = sign,
       x = NULL, y = "UCell score") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()


#### ---- HM of ifn and inflam genes ---- ####
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# ── 1. Combine IFN and Inflammatory genes ─────────────────────────────────────
all_sig_genes <- union(ifn_genes, inflam_genes)

# Keep only genes present in the assay
all_sig_genes <- all_sig_genes[all_sig_genes %in% rownames(obj[["imputed_t2"]])]

# ── 2. Pseudobulk with AverageExpression ──────────────────────────────────────
avg_exp <- AverageExpression(
  obj_sub,
  assays      = "imputed_t2",
  features    = all_sig_genes,
  group.by    = c("colonization", "stimulation", "orig.ident", "ID_labs_ext"),
  return.seurat = FALSE
)$imputed_t2

# ── 3. Extract metadata per pseudobulk sample ─────────────────────────────────
# Column names are "colonization_stimulation_orig.ident_ID_labs_ext"
meta_df <- data.frame(label = colnames(avg_exp)) %>%
  separate(label, into = c("colonization", "stimulation", "orig.ident", "ID_labs_ext"),
           sep = "_", remove = FALSE) %>%
  column_to_rownames("label")

# ── 4. Scale expression matrix (z-score per gene) ────────────────────────────
mat <- avg_exp[all_sig_genes, , drop = FALSE] # t(scale(t(avg_exp)))
mat <- mat[, order(meta_df$ID_labs_ext)]  # order columns by ID_labs_ext
meta_df <- meta_df[order(meta_df$ID_labs_ext), ]

# ── 5. Define colors ──────────────────────────────────────────────────────────
# Adjust these to match your actual levels
col_colonization <- col_cols
col_stimulation <- stim_cols


# ── 6. Row annotation — gene signature membership ─────────────────────────────
row_anno_df <- data.frame(
  signature = case_when(
    rownames(mat) %in% ifn_genes    & rownames(mat) %in% inflam_genes ~ "Both",
    rownames(mat) %in% ifn_genes                                       ~ "IFN Alpha",
    rownames(mat) %in% inflam_genes                                    ~ "Inflammatory"
  ),
  row.names = rownames(mat)
)

row_anno <- rowAnnotation(
  Signature = row_anno_df$signature,
  col = list(Signature = c(
    "IFN Alpha"     = "#7B9BC8",
    "Inflammatory"  = "#E8856A",
    "Both"          = "#9B6BB5"
  )),
  show_legend = TRUE
)

# ── 7. Column annotation — colonization and stimulation color bars ─────────────
col_anno <- HeatmapAnnotation(
  Colonization = meta_df$colonization,
  Stimulation  = meta_df$stimulation,
  col = list(
    Colonization = col_colonization,
    Stimulation  = col_stimulation
  ),
  annotation_name_side = "left"
)

# ── 8. Column split by ID_labs_ext ────────────────────────────────────────────
col_split <- factor(meta_df$ID_labs_ext, levels = unique(meta_df$ID_labs_ext))

# ── 9. Draw heatmap ───────────────────────────────────────────────────────────
ht <- Heatmap(
  as.matrix(mat),
  name                  = "Z-score",
  col                   = colorRamp2(c(min(mat_sig), median(mat_sig), max(mat_sig)),c("#4393C3", "white", "#D6604D")),#colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D")),
  top_annotation        = col_anno,
  left_annotation       = row_anno,
  column_split          = col_split,
  cluster_rows          = TRUE,
  cluster_columns       = FALSE,   # grouped by ID_labs_ext, no clustering within
  cluster_column_slices = FALSE,
  show_column_names     = FALSE,
  row_names_gp          = gpar(fontsize = 7),
  column_title_gp       = gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param  = list(direction = "vertical")
)

pdf(paste0(outdir,dato,"MSigDB_InflamIFNAlphaSignature_HMunscaled.pdf"),height = 25, width = 8 )
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#### ---- separate HM of ifn and inflam genes ---- ####
# ── 1. Heatmap function ───────────────────────────────────────────────────────
plot_sig_heatmap <- function(genes, title) {
  
  # Filter to signature genes present in assay
  sig_genes <- intersect(genes, rownames(avg_exp))
  mat_sig   <- avg_exp[sig_genes, , drop = FALSE]
  #mat_sig   <- t(scale(t(mat_sig)))
  mat_sig   <- mat_sig[, order(meta_df$ID_labs_ext), drop = FALSE]
  
  Heatmap(
    as.matrix(mat_sig),
    name                  = "Z-score",
    #column_title          = title,
    column_title_gp       = gpar(fontsize = 11, fontface = "bold"),
    col                   = colorRamp2(c(min(mat_sig), median(mat_sig), max(mat_sig)),c("#4393C3", "white", "#D6604D")),#colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D")),
    top_annotation        = col_anno,
    column_split          = col_split,
    cluster_rows          = TRUE,
    cluster_columns       = FALSE,
    cluster_column_slices = FALSE,
    show_column_names     = FALSE,
    row_names_gp          = gpar(fontsize = 7),
    heatmap_legend_param  = list(direction = "vertical")
  )
}


# ── 2. Draw each heatmap ──────────────────────────────────────────────────────
ht_ifn <- plot_sig_heatmap(ifn_genes,    "HALLMARK IFN Alpha Response")
ht_inf <- plot_sig_heatmap(inflam_genes, "HALLMARK Inflammatory Response")

pdf(paste0(outdir,dato,"MSigDB_IFNAlphaSignature_HMunscaled.pdf"),height = 15, width = 8 )
draw(ht_ifn, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
pdf(paste0(outdir,dato,"MSigDB_InflamSignature_HMunscaled.pdf"),height = 18, width = 8 )
draw(ht_inf, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

ht_isg <- plot_sig_heatmap(IFNaup_genes, "ISG")
draw(ht_isg, heatmap_legend_side = "right", annotation_legend_side = "right")

ht_lps <- plot_sig_heatmap(NemLPS_genes, "Nemeth LPS up")
draw(ht_lps, heatmap_legend_side = "right", annotation_legend_side = "right")

