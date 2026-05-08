#' R script for saving plots in figure format
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on hhttps://stuartlab.org/signac/articles/pbmc_vignette.html

#### ---- Initiate libraries ---- ####
library(ggplot2)
library(stringr)
library(ggrastr)
library(viridis)
library(scales)
library(Signac)
library(Seurat)
library(biovizBase)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges) 
library(colorRamp2)
library(scales)
library(matrixStats)
library(openxlsx)
library(ggrastr)


#### ---- variables used throughout script ---- ####
rm(list = ls())
projdir <- getwd()
RAID_dir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff"
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
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/output/UMAPs/"
## sample combination
project <- "BM-PBSvsHA107-PBSvsLPS-8wk"

#### ---- read in data ---- ####
obj <- readRDS("25_10_06_PBSHA107PBALPS_8wk_clean.rds")
comb_df <- cbind(obj@meta.data, obj@reductions$umap@cell.embeddings)
head(comb_df)

#### ---- Initial umaps and annotation ---- ####
## With Cell annotations (wide)
pdf(paste(outdir,dato,"_UMAP_CellAnnotation.pdf",sep = ""),height = 3, width = 5)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=ID_labs))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  labs(x="UMAP1", y="UMAP2", colour="Cell annotation")+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  theme()
dev.off()

res2col <- hue_pal()(7)
names(res2col) <- seq(0,6,1)

## Clustering
pdf(paste(outdir,dato,"_UMAP_ClusteringRes.0.2.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=ATAC_snn_res.0.2))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  labs(x="UMAP1", y="UMAP2", colour="res.0.2")+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  theme()
dev.off()

#### Annotation related
## Monocyte separation - feature plots on UMAP
pdf(paste(outdir,dato,"_UMAP_Ly6c2_imp6gene.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=obj@assays$imputed_t6@data["Ly6c2",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ly6c2")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_Ly6c1_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=obj@assays$imputed_t2@data["Ly6c1",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ly6c1")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_Ccr2_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=obj@assays$RNA@data["Ccr2",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ccr2")+
  theme()
dev.off()

## LSK seperation - featuresplots on UMAP
pdf(paste(outdir,dato,"_UMAP_Cd34_imp2gene.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=obj@assays$imputed_t2@data["Cd34",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Cd34")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_cKit_imp2gene.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=obj@assays$imputed_t2@data["Kit",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Kit")+
  theme()
dev.off()


## Modulescore vln plots based on the NicheNetData
# monocytes
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_Monocytes.pdf",sep = ""),height = 3, width = 4)
VlnPlot(obj, features = "Monocytes1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_Monocytes.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Monocytes1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# neutrophils
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_Neutrophils.pdf",sep = ""),height = 3, width = 4)
VlnPlot(obj, features = "Neutrophils1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_Neutrophils.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Neutrophils1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# DCs
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_DendriticCells.pdf",sep = ""),height = 3, width = 4)
VlnPlot(obj, features = "Dendritic.cells1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_DendriticCells.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Dendritic.cells1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# NK cells
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_NKcells.pdf",sep = ""),height = 3, width = 4)
VlnPlot(obj, features = "NK.cells1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_NKcells.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "NK.cells1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

### Imp gene exp - vln plots
pdf(paste(outdir,dato,"_Vln_res.0.2_imp6Ly6c2.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Ly6c2",
        group.by = "ATAC_snn_res.0.2",
        assay = "imputed_t6",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_res.0.2_imp6Kit.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Kit",
        group.by = "ATAC_snn_res.0.2",
        assay = "imputed_t6",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_res.0.2_imp6Cd34.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj, features = "Cd34",
        group.by = "ATAC_snn_res.0.2",
        assay = "imputed_t6",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

#### ---- cell type distribution ---- ####
dist_df <- perc_function_samp("ID_labs", Cells(obj), obj,"orig.ident")
dist_df$samp <- factor(dist_df$samp, 
                       levels = c("BM-PBS-PBS-21d","BM-PBS-PBS-8wk","BM-HA107-PBS-21d","BM-HA107-PBS-8wk","BM-HA107-LPS-21d","BM-HA107-LPS-8wk","BM-PBS-LPS-21d","BM-PBS-LPS-8wk"))

pdf(paste(outdir,dato,"_Distribution_PersampleStacked.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black")+
  theme_classic()+
  labs(x="", y="% of sample")+
  theme(axis.text.x = element_text(angle=90))
dev.off()

pdf(paste(outdir,dato,"_Distribution_PersampleDodge.pdf",sep = ""),height = 4, width = 6)
ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  facet_grid(.~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()

len <- dim(dist_df)[1]
dist_df$colonization <- factor(unlist(str_split(dist_df$samp,"-"))[seq(2,len*4,4)], levels=c("PBS","HA107"))
dist_df$timepoint <- unlist(str_split(dist_df$samp,"-"))[seq(4,len*4,4)]
dist_df$stimulation <- factor(unlist(str_split(dist_df$samp,"-"))[seq(3,len*4,4)], levels=c("PBS","LPS"))

pdf(paste(outdir,dato,"_Distribution_PersampleStacked_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=colonization, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black", width=0.8)+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~stimulation)+
  theme(axis.text.x = element_text(angle=90))
dev.off()
pdf(paste(outdir,dato,"_Distribution_PersampleDodge_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=colonization, y=percent, fill=stimulation))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()
pdf(paste(outdir,dato,"_Distribution_PersampleDodge_StimxColxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=stimulation, y=percent, fill=colonization))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()


saveRDS(dist_df, paste(dato,"PBSHA107PBALPS_8wk_clean_distdf.rds",sep="_"))

#### mono trajectory gene exp ####
visu <- readRDS("2510_BM-HA107PBS-LPSPBS-8wk_MonocyteTraj_visudf.rds")

head(visu)

## Monocyte separation - feature plots on UMAP
pdf(paste(outdir,dato,"_tPC2xTmean_res.0.2.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=ATAC_snn_res.0.2))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_color_manual(values = res2col)+
  labs(x="Av. Pseudotime", y="tPC2", colour="res.0.2")+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_IDlabs.pdf",sep = ""),height = 3, width = 5)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=ID_labs))+
  geom_point_rast(size=1)+
  theme_classic()+
  labs(x="Av. Pseudotime", y="tPC2", colour="ID_labs")+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Ly6c2_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t6@data["Ly6c2",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Ly6c2")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Ly6c1_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t6@data["Ly6c1",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Ly6c1")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Ccr2_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t2@data["Ccr2",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Ccr2")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Cd26LSELL_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t2@data["Sell",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Sell")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Cd34_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t6@data["Cd34",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Cd34")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Cx3cr1_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t2@data["Cx3cr1",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Cx3cr1")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Kit_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t6@data["Kit",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Kit")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_Spn_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@assays$imputed_t2@data["Spn",rownames(visu)]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Spn")+
  theme()
dev.off()

#### Monocyte trajectory - PBMC and BM signatures ####
# from https://www.sciencedirect.com/science/article/pii/S1074761317301838#mmc1, Fig 1C
# Cl IV
Ly6chiblood <- c("Fcrls","Irf7","Trem2","Ccr2","Igals3","Ifi30")
# Cl VI
Matmono <- c("Apoe","Cd36","Csf1r","Cx3cr1","Ciita","Itgax","H2-aa","Fcgr4","Nr4a1","Cd74","H2-ab1")

obj <- AddModuleScore(obj, features = list(Ly6chiblood), ctrl = 100, name = "Ly6hiblood", assay = "imputed_t6")
obj <- AddModuleScore(obj, features = list(Matmono), ctrl = 100, name = "Matmonoblood", assay = "imputed_t2")

pdf(paste(outdir,dato,"_tPC2xTmean_Ly6chiBlood.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@meta.data[rownames(visu),"Ly6hiblood1"]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Blood\nLy6c high")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_tPC2xTmean_MatMonoBlood.pdf",sep = ""),height = 3, width = 4)
ggplot(visu, aes(x=T_mean, y=tPC2, colour=obj@meta.data[rownames(visu),"Matmonoblood1"]))+
  geom_point_rast(size=1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="Av. Pseudotime", y="tPC2", colour="Blood\nMat. mono")+
  theme()
dev.off()


obj_sub <- subset(obj, idents = c("Q1","Q2","Q3"))
pdf(paste(outdir,dato,"_VlnQ1-Q3_Mildneretal_Ly6chiblood.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj_sub, features = "Ly6hiblood1",
        group.by = "ID_labs_ext",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_Mildneretal_Matmonoblood.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj_sub, features = "Matmonoblood1",
        group.by = "ID_labs_ext",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp6Cd34.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj_sub, features = "Kit",
        group.by = "ID_labs_ext",assay="imputed_t2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp6Cd34.pdf",sep = ""),height = 4, width = 4)
VlnPlot(obj_sub, features = "Ccr2",split.by = "orig.ident",
        group.by = "ID_labs_ext",assay="imputed_t2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp2Cox8a.pdf",sep = ""),height = 4, width = 6)
VlnPlot(obj_sub, features = "Cox8a",split.by = "orig.ident",
        group.by = "ID_labs_ext",assay="imputed_t2",
        pt.size = 0)
  xlab("")
dev.off()


#### Vln and coverageplot for DAR outputs ####
obj@meta.data$orig.ident <- factor(obj$orig.ident, levels =
                                     c("BM-PBS-PBS-21d","BM-PBS-LPS-21d","BM-HA107-PBS-21d","BM-HA107-LPS-21d",
                                       "BM-PBS-PBS-8wk","BM-PBS-LPS-8wk","BM-HA107-PBS-8wk","BM-HA107-LPS-8wk" ))
obj@meta.data$ID_labs_ext <- factor(obj$ID_labs_ext, levels =
                                     c("Q1","Q2","Q3","Ly6c lo monocytes","Ly6c hi monocytes","LSK", "Neutrophils","Dendritic cells", "NK cells"))
Idents(obj) <- 'ID_labs_ext'

VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Hecw2", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")

VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Rhob", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Itgb1", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Rhoa", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Rac1", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Cdc42", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")

VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Fosb", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Fos", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Jun", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Jund", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")

VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Dsel", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Cspg5", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Hexb", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")

VlnPlot(subset(obj, idents = c("Q1","Q2","Q3")), features = "Ifngr2", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")

#### Cox8a coverage plot ####
# Define the region to highlight
# chr16-5473379-5474330
# chr3-105803366-105804286
# chr19-7217216-7218118
highlight_1 <- GRanges(seqnames = "chr19", ranges = IRanges(start =7217216, end = 7218118))
highlight_2 <- GRanges(seqnames = "chr16", ranges = IRanges(start =5473379, end = 5474330))
highlight_3 <- GRanges(seqnames = "chr3", ranges = IRanges(start =105803366, end = 105804286))

pdf(paste0(outdir,dato,"_",project,"_CovPlotCox8aChr19_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
  region = c("Cox8a"),
  group.by = "orig.ident",
  extend.upstream = 100,
  extend.downstream = 1000,
  ncol = 1,
  region.highlight = highlight_1 #c(highlight_1,highlight_2,highlight_3)
)
dev.off()
pdf(paste0(outdir,dato,"_",project,"_CovPlotRap1aChr3_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region = "Rap1a",
             group.by = "orig.ident",
             extend.upstream = 1000,
             extend.downstream = 3000,
             ncol = 1,
             region.highlight = highlight_3)
dev.off()

## CS / DS metabolism 
CDDS_met <- c( "Hexb", "Cspg5", "Dsel")

obj_sub@meta.data$samp <- factor(obj_sub$orig.ident,levels=rev(levels(obj_sub$orig.ident)))


highlight_3 <- GRanges(seqnames = "chr13", ranges = IRanges(start =97190128, end = 97191027))
pdf(paste0(outdir,dato,"_",project,"_CovPlotHexbChr13_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region = "Hexb",
             group.by = "samp", 
             extend.upstream = 1000,
             extend.downstream = 3000,
             ncol = 1,
             region.highlight = highlight_3)
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp2Hexb.pdf",sep = ""),height = 4, width = 6)
VlnPlot(obj_sub, features = "Hexb", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
dev.off()

highlight_4 <- GRanges(seqnames = "chr9", ranges = IRanges(start =110280869, end = 110281768))
pdf(paste0(outdir,dato,"_",project,"_CovPlotCspg5Chr9_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region=c("chr9-110280869-110281768"),
             #region = "Cspg5",
             group.by = "samp",
             extend.upstream = 20000,
             extend.downstream = 1000,
             region.highlight = highlight_4,
             ncol = 1)
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp2Cspg5.pdf",sep = ""),height = 4, width = 6)
VlnPlot(obj_sub, features = "Cspg5", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
dev.off()


highlight_5 <- GRanges(seqnames = "chr1", ranges = IRanges(start =112457144, end = 112458038))
pdf(paste0(outdir,dato,"_",project,"_CovPlotDselChr1_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region=c("chr1-112457144-112458038"),
             #region = "Dsel",
             group.by = "samp",
             extend.upstream = 1000,
             extend.downstream = 1000,
             region.highlight = highlight_5,
             ncol = 1)
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp2Dsel.pdf",sep = ""),height = 4, width = 6)
VlnPlot(obj_sub, features = "Dsel", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
dev.off()

### Cxcr4
#chr1-128614707-128615647
highlight_6 <- GRanges(seqnames = "chr1", ranges = IRanges(start =128614707, end = 128615647))
pdf(paste0(outdir,dato,"_",project,"_CovPlotCxcr4Chr1_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region=c("chr1-128614707-128615647"),
             region = "Cxcr4",
             group.by = "samp",
             extend.upstream = 2000,
             extend.downstream = 25000,
             region.highlight = highlight_6,
             ncol = 1)+geom_hline(yintercept = 2000)
dev.off()
pdf(paste(outdir,dato,"_VlnQ1-Q3_imp2DCxcr4.pdf",sep = ""),height = 4, width = 6)
VlnPlot(obj_sub, features = "Cxcr4", assay = "imputed_t2", pt.size=0, split.by = "orig.ident")
dev.off()


highlight_6 <- GRanges(seqnames = "chr2", ranges = IRanges(start =102898871, end = 102899803))
pdf(paste0(outdir,dato,"_",project,"_CovPlotChr2Cd44_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region = "Cd44",
             group.by = "samp",
             extend.upstream = 2000,
             extend.downstream = 2000,
             region.highlight = highlight_6,
             ncol = 1)+geom_hline(yintercept = 2000)
dev.off()

## mito metabolism
highlight_6 <- GRanges(seqnames = "chr19", ranges = IRanges(start =43674671, end = 43675551))
pdf(paste0(outdir,dato,"_",project,"_CovPlotChr1_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             region = c("chr2-102898871-102899803"),
             #region = "Slc25a28",
             group.by = "samp",
             extend.upstream = 2000,
             extend.downstream = 2000,
             region.highlight = highlight_6,
             ncol = 1)+geom_hline(yintercept = 2000)
dev.off()

### inflam. response
highlight_7 <- GRanges(seqnames = "chr16", ranges = IRanges(start =91534221, end = 91535120))
pdf(paste0(outdir,dato,"_",project,"_CovPlotChr16Ifngr2_Q1-Q3.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region = c("chr16-91534221-91535120"),
             region = "Ifngr2",
             group.by = "samp",
             extend.upstream = 15000,
             extend.downstream = 2000,
             region.highlight = highlight_7,
             ncol = 1)+geom_hline(yintercept = 2000)
dev.off()


highlight_7 <- GRanges(seqnames = "chr6", ranges = IRanges(start =124477886, end = 124478765))
pdf(paste0(outdir,dato,"_",project,"_CovPlotChr6Clstn3_Q1-Q3.pdf"),height = 6, width = 4 )
p1 <- CoveragePlot(obj_sub,
             #region = c("chr6-124477886-124478765"),
             region = "Clstn3",
             group.by = "samp",
             extend.upstream = 1000,
             extend.downstream =20000,
             region.highlight = highlight_7,
             ncol = 1)
dev.off()



#### Custom violin plots ####
obj_sub$colonization <- factor(obj_sub$colonization, levels=c("PBS","HA107"))
col_cols <- c("PBS"="#F8766D","HA107"="#00BFC4")
stim_cols <- c("PBS"="white","LPS"="lightgrey")
# Step 1: Calculate mean or median per group
gene <- "Kmte"
df <- cbind(obj_sub@meta.data, val = obj_sub@assays$imputed_t2@data[gene,])

summary_stats <- df %>%
  group_by(orig.ident,colonization,stimulation, ID_labs_ext) %>%
  summarise(
    mean_val = mean(val),
    median_val = median(val)
  )

# Step 2: Create violin plot with line connecting means (or medians)
pdf(paste0(outdir,dato,"_VlnQ1-Q3_imp2",gene,".pdf"),height = 4, width = 7)
ggplot(df, aes(x = orig.ident, y = val, fill = stimulation, color = colonization)) +
  geom_violin(trim = FALSE) +
  geom_point(data = summary_stats, aes(y = mean_val), color = "black", size = 2) + # points at mean
  geom_line(data = summary_stats, aes(y = mean_val, group = colonization), color = "black", size = 1) + # line connecting means
  scale_color_manual(values = col_cols)+scale_fill_manual(values = stim_cols)+
  theme_classic() +# geom_hline(yintercept = c(0.57,0.8))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_grid(~ID_labs_ext)+
  labs(y = "Imputed Gene Expression", x = "Monocyte sub cluster", title = gene)
dev.off()

#### Coverage and imp. vlns based on lists ####
# dir to save to
outdir2 <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/commonDAR_cov&vln/"
# violin colors
obj_sub$colonization <- factor(obj_sub$colonization, levels=c("PBS","HA107"))
col_cols <- c("PBS"="#F8766D","HA107"="#00BFC4")
stim_cols <- c("PBS"="white","LPS"="lightgrey")
# GRanges object structure
GRobj <- granges(obj_sub)


# PBS-PBS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSPBS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSPBS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "PBSPBS"

# HA107-PBS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107PBS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107PBS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "HA107PBS"

# PBS-LPS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSLPS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSLPS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "PBSLPS"

# HA107-LPS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107LPS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107LPS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "HA107LPS"


## Plotting each group by exceutong above 6 lines (check everything is correct) and run loop
not_there <- c()
for (gene in gens){
  print(gene)
  # Covergae plots
  highlight_regs <- GRobj[0]
  for (i in seq(1,length(hits))[hits %in% gene]){
    print(regs[i])
    hreg <- unlist(str_split(regs[i],"-"))
    hregs <- GRanges(seqnames = hreg[1], ranges = IRanges(start=as.numeric(hreg[2]), end=as.numeric(hreg[3])))
    highlight_regs <- c(highlight_regs, hregs)
    #pdf(paste0(outdir2,dato,"_",org,"_",gene,"_reg",regs[i],"_CovPlot.pdf"),height = 6, width = 4 )
    p1 <- CoveragePlot(obj_sub,
                       region = c(regs[i]),
                       #region = gene,
                       group.by = "samp",
                       extend.upstream = 2000,
                       extend.downstream =2000,
                       region.highlight = hregs,
                       ncol = 1)
    #print(p1)
    #dev.off()
  }
  pdf(paste0(outdir2,dato,"_",org,"_",gene,"_CovPlot.pdf"),height = 6, width = 4 )
  p1 <- CoveragePlot(obj_sub,
                   region = gene,
                   group.by = "samp",
                   extend.upstream = 8000,
                   extend.downstream =8000,
                   region.highlight = highlight_regs,
                   ncol = 1)
  print(p1)
  dev.off()

  ## Violins
  # first check gene exists in imp2 data
  if (gene %in% rownames(obj_sub@assays$imputed_t2@data)){
  # Step 1: Calculate mean or median per group
  df <- cbind(obj_sub@meta.data, val = obj_sub@assays$imputed_t2@data[gene,])
  
  summary_stats <- df %>%
    group_by(orig.ident,colonization,stimulation, ID_labs_ext) %>%
    summarise(
      mean_val = mean(val),
      median_val = median(val)
    )
  
  # Step 2: Create violin plot with line connecting means (or medians)
  pdf(paste0(outdir2,dato,"_",org,"_",gene,"_VlnQ1-Q3_imp2.pdf"),height = 4, width = 7)
  p2 <- ggplot(df, aes(x = orig.ident, y = val, fill = stimulation, color = colonization)) +
    geom_violin(trim = FALSE) +
    geom_point(data = summary_stats, aes(y = mean_val), color = "black", size = 2) + # points at mean
    geom_line(data = summary_stats, aes(y = mean_val, group = colonization), color = "black", size = 1) + # line connecting means
    scale_color_manual(values = col_cols)+scale_fill_manual(values = stim_cols)+
    theme_classic()+ 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    facet_grid(~ID_labs_ext)+
    labs(y = "Imputed Gene Expression", x = "Monocyte sub cluster", title = gene)
  print(p2)
  dev.off()
  } else {not_there <- c(not_there,gene)}
  }
print("Sorry,following genes were not in the imp data and couldn't be plotted with violin plots:")
print(not_there)


#### Extras ####
gene <- "Atxn1" #Etv4, Etv3, Egfr
df <- cbind(obj_sub@meta.data, val = obj_sub@assays$imputed_t2@data[gene,])

summary_stats <- df %>%
  group_by(orig.ident,colonization,stimulation, ID_labs_ext) %>%
  summarise(
    mean_val = mean(val),
    median_val = median(val)
  )

# Step 2: Create violin plot with line connecting means (or medians)
#pdf(paste0(outdir2,dato,"_",org,"_",gene,"_VlnQ1-Q3_imp2.pdf"),height = 4, width = 7)
p2 <- ggplot(df, aes(x = orig.ident, y = val, fill = stimulation, color = colonization)) +
  geom_violin(trim = FALSE) +
  geom_point(data = summary_stats, aes(y = mean_val), color = "black", size = 2) + # points at mean
  geom_line(data = summary_stats, aes(y = mean_val, group = colonization), color = "black", size = 1) + # line connecting means
  scale_color_manual(values = col_cols)+scale_fill_manual(values = stim_cols)+
  theme_classic()+ 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_grid(~ID_labs_ext)+
  labs(y = "Imputed Gene Expression", x = "Monocyte sub cluster", title = gene)
print(p2)
#dev.off()



# if you want the highlights included, define gene and reset to correct origin of DAR then run run 656-673
pdf(paste0(outdir2,dato,"_PBSPBS_Glipr2_v2_CovPlot.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region = c("chr16-91534221-91535120"),
             region = "Glipr2",
             group.by = "samp",
             extend.upstream = -14000,
             extend.downstream = 2000,
             ncol = 1)
dev.off()

pdf(paste0(outdir2,dato,"_PBSPBS_Jund_reg_v2_CovPlot.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region = c("chr16-91534221-91535120"),
             region = "Jund",
             group.by = "samp",
             extend.upstream = 20000,
             extend.downstream = 20000,
             region.highlight = highlight_regs,
             ncol = 1)
dev.off()

pdf(paste0(outdir2,dato,"_PBSPBS_Sept9_reg_v2_CovPlot.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region = c("chr16-91534221-91535120"),
             region = "Sept9",
             group.by = "samp",
             extend.upstream = -115000,
             extend.downstream = -34000,
             region.highlight = highlight_regs,
             ncol = 1)
dev.off()

pdf(paste0(outdir2,dato,"_PBSPBS_Cic_reg_v2_CovPlot.pdf"),height = 6, width = 4 )
CoveragePlot(obj_sub,
             #region = c("chr16-91534221-91535120"),
             region = "Cic",
             group.by = "samp",
             extend.upstream = 1000,
             extend.downstream = -10000,
             region.highlight = highlight_regs,
             ncol = 1)
dev.off()

