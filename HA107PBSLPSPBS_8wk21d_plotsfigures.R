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
outdir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/BM-PBSvsHA107-PBSvsLPS-21dvs8wk/output/UMAPv2clean/"
## sample combination
project <- "BM-PBSvsHA107-PBSvsLPS-21dvs8wk"

#### ---- read in data ---- ####
combined <- readRDS("25_07_16_PBSHA107PBALPS_8wk21d_clean_v2.rds")
comb_df <- cbind(combined@meta.data, combined@reductions$umap@cell.embeddings)
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
pdf(paste(outdir,dato,"_UMAP_Ly6c2_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=combined@assays$RNA@data["Ly6c2",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ly6c2")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_Ly6c1_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=combined@assays$RNA@data["Ly6c1",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ly6c1")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_Ccr2_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=combined@assays$RNA@data["Ccr2",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Ccr2")+
  theme()
dev.off()

## LSK seperation - featuresplots on UMAP
pdf(paste(outdir,dato,"_UMAP_Cd34_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=combined@assays$RNA@data["Cd34",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Cd34")+
  theme()
dev.off()
pdf(paste(outdir,dato,"_UMAP_cKit_geneact.pdf",sep = ""),height = 3, width = 4)
ggplot(comb_df, aes(x=umap_1, y=umap_2, colour=combined@assays$RNA@data["Kit",]))+
  geom_point_rast(size=0.1)+
  theme_classic()+
  scale_colour_gradientn(colors = mycols)+
  labs(x="UMAP1", y="UMAP2", colour="Kit")+
  theme()
dev.off()


## Modulescore vln plots based on the NicheNetData
# monocytes
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_Monocytes.pdf",sep = ""),height = 3, width = 4)
VlnPlot(combined, features = "Monocytes1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_Monocytes.pdf",sep = ""),height = 4, width = 4)
VlnPlot(combined, features = "Monocytes1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# neutrophils
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_Neutrophils.pdf",sep = ""),height = 3, width = 4)
VlnPlot(combined, features = "Neutrophils1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_Neutrophils.pdf",sep = ""),height = 4, width = 4)
VlnPlot(combined, features = "Neutrophils1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# DCs
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_DendriticCells.pdf",sep = ""),height = 3, width = 4)
VlnPlot(combined, features = "Dendritic.cells1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_DendriticCells.pdf",sep = ""),height = 4, width = 4)
VlnPlot(combined, features = "Dendritic.cells1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()

# NK cells
pdf(paste(outdir,dato,"_Vln_res.0.2_Baccinetal_NKcells.pdf",sep = ""),height = 3, width = 4)
VlnPlot(combined, features = "NK.cells1",
        group.by = "ATAC_snn_res.0.2",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()
pdf(paste(outdir,dato,"_Vln_IDlabs_Baccinetal_NKcells.pdf",sep = ""),height = 4, width = 4)
VlnPlot(combined, features = "NK.cells1",
        group.by = "ID_labs",
        pt.size = 0)+NoLegend()+
  xlab("")
dev.off()



#### ---- cell type distribution ---- ####
dist_df <- perc_function_samp("ID_labs", Cells(combined), combined,"orig.ident")
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

saveRDS(dist_df, paste(dato,"PBSHA107PBALPS_8wk21d_clean_v2_distdf.rds",sep="_"))

#### ---- DAR analysis outputs ---- ####
