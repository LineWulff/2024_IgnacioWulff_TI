#' R script for isolating monocytes and prep data for tSPACE
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
## sample combination
project <- "BM-PBSvsHA107-PBSvsLPS-21dvs8wk"



#### ---- Read in data and subset it ---- ####
combined <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/25_07_16_PBSHA107PBALPS_8wk21d_clean_v2.rds")
Idents(combined) <- 'ID_labs'
DefaultAssay(combined) <- 'ATAC'

sheets <- c("cluster_Ly6c lo monocytes","cluster_Neutrophils","cluster_Dendritic cells","cluster_LSK","cluster_Ly6c hi monocytes","cluster_NK cells")
incl_subs <- c("Ly6c lo monocytes","Ly6c hi monocytes")

# subset cells
combined <- subset(combined, idents = incl_subs)

#### ---- Prep for tSpace - test  ---- ####
table(combined@meta.data$orig.ident)/42051*100
table(combined@meta.data$ID_labs)/42051*100
## For test purposes make downsampled version of dataset
down_comb <- sample(Cells(combined), size = 1000)
down_comb <- subset(combined, cells = down_comb)
table(down_comb@meta.data$orig.ident)/10
table(down_comb@meta.data$ID_labs)/10

saveRDS(down_comb, file = "250818_DS_monocytes.rds")

#### ---- Prep for tSpace - whole object  ---- ####

combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="timepoint")
DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="colonization")
DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="stimulation")

DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1)


comb_low <- RunTFIDF(comb_low)
comb_low <- FindTopFeatures(comb_low, min.cutoff = "q95", assay = "ATAC")
length(VariableFeatures(comb_low)) #208510

combined[["ATAC"]]
combined[["RNA"]]
comb_low <- combined
comb_low@active.assay <- "ATAC"
comb_low[["RNA"]] <- NULL #for smaller version
saveRDS(comb_low, file = "250822_NoRNA_monocytes.rds")

#### ---- Open and adjust tspace objects ---- ####
subs1 <- readRDS(file="/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/tspace_output/250822_NoAssays_monocytes_tspacefile.rds")

head(subs1$ts_file)
visu <- subs1$ts_file
visu <- cbind(visu, combined@meta.data[rownames(visu),])
#visu <- visu[!visu$F2F5 %in% c("cDC2","cDC3","amb","HLA low"),]
#visu <- cbind(visu, DC_subs@reductions$umap@cell.embeddings[rownames(visu),])
#visu <- cbind(visu, idents=Idents(DC_subs)[rownames(visu)])


#distance from monocytes
#mono.trajectories <- subs1$tspace_matrix[,which(colnames(subs1$tspace_matrix) %in% paste0('T_', visu[which(visu$integrated_snn_res.2.8=='24'), 'Index']))]
#colnames(mono.trajectories) <- c("T_1","T_2","T_3","T_4","T_5")
#visu <- cbind(visu, dist_mono=mono.trajectories)
#visu$T_mean <-rowMeans(mono.trajectories)

#### Adding new trajectory based clustering ####
visu_obj <- subset(combined, cells = rownames(visu))
visu_obj@reductions$lsi@cell.embeddings <- as.matrix(visu[,startsWith(colnames(visu),"LSI")])

visu_obj <- RunUMAP(object = visu_obj, reduction = 'lsi', dims = 2:30)

DimPlot(visu_obj, group.by = "ID_labs")

ggplot(visu, aes(x=umap1, y=umap2, color=ID_labs))+
  geom_point()
