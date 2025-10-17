#' R script for inspecting Sahers analysis object
#' Author: Line Wulff
#' Date (created): 25-07-02
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
library(tidyverse)
library(tidyr)
#new test
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

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
project <- "BM-HA107PBSLPSPBS-8wk"


#### ---- read in data ---- ####
# 4 samples of 8wk mice
seu_obj <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/21D-integrated_multiome.rds")

#### ---- IDs ---- ####
DimPlot(seu_obj,
        label = T)

FeaturePlot(
  object = seu_obj,
  features = c("Ly6c2"),
  pt.size = 0.1
)

# Cluster IDs (mine):
# 0 - monocytes - Ly6c hi
# 1 - monocytes - Ly6c hi
# 2 - neutrophils 
# 3 - monocytes - Ly6c hi
# 4 - LSK
# 5 - monocytes - Ly6c hi
# 6 - monocytes - Ly6c lo
# 7 - cDC/cDC progenitors
# 8 - Eosinophils?
# 9 - NK cells

## Collapse and label
Ly6chi_mono <- c(0,1,3,5)

seu_obj@meta.data$ID_labs <- NA
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% Ly6chi_mono,]$ID_labs <- "Ly6c hi monocytes"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(6),]$ID_labs <- "Ly6c lo monocytes"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(2),]$ID_labs <- "Neutrophils"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(4),]$ID_labs <- "LSK"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(7),]$ID_labs <- "cDC"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(8),]$ID_labs <- "Eosinophils"
seu_obj@meta.data[seu_obj@meta.data$integrated_snn_res.0.4 %in% c(9),]$ID_labs <- "NK cells"

DimPlot(seu_obj,
        group.by = "ID_labs",
        split.by = "dataset",
        ncol = 2)

## cluster distributon between samples
clus_dist <- perc_function("ID_labs", Cells(seu_obj), seu_obj)
clus_dist
clus_dist <- perc_function_samp("ID_labs", Cells(seu_obj), seu_obj,"dataset")
clus_dist

ggplot(clus_dist, aes(x=cluster, y=percent, fill=samp))+
  geom_bar(stat = "identity", position = "dodge", colour="black")+
  theme_classic()+
  ylab("% of sample")+
  theme(axis.text.x = element_text(angle=90))

ggplot(clus_dist, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat = "identity", colour="black")+
  theme_classic()+
  ylab("% of sample")+
  theme(axis.text.x = element_text(angle=90))

#### Save meta data from this project for future use in annotation
head(seu_obj@meta.data)

# what is wrong woth pbs_pbs data names?
seu_meta[seu_meta$dataset=="TPBS_PBS",] #single names
seu_meta[seu_meta$dataset=="PBS_LPS",] #double names
seu_meta[seu_meta$dataset=="HA107_LPS",] #double names
seu_meta[seu_meta$dataset=="HA107_PBS",] #double names

#change the rownames first
seu_meta <- seu_obj@meta.data
seu_meta$barcode <- NA
seu_meta[seu_meta$dataset=="TPBS_PBS",]$barcode <- unlist(str_split(rownames(seu_meta[seu_meta$dataset=="TPBS_PBS",]),"_"))[seq(2,5971*2,2)]
seu_meta[!seu_meta$dataset=="TPBS_PBS",]$barcode <- unlist(str_split(rownames(seu_meta[!seu_meta$dataset=="TPBS_PBS",]),"_"))[seq(3,(33938-5971)*3,3)]

rownames(seu_meta) <- paste(seu_meta$dataset, seu_meta$barcode, sep = "-")

write.csv(seu_meta, file = paste(projdir,"/Outputs/Clustering/",dato,"8wkcombined_metadata_SaherAnalysis.csv", sep = ""))
