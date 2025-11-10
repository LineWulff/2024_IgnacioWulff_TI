#' R script for add RNAseq data to established ATACseq object
#' Author: Line Wulff
#' Date (created): 25-10-09

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
library(ChIPseeker)


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
project <- "BM-PBSvsHA107-PBSvsLPS-8wk"

#### --- Read in data ---- ####
combined <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/25_10_06_PBSHA107PBALPS_8wk_clean.rds")
#comb_full <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/25_08_11_PBSHA107PBALPS_8wk21d_clean_v2.rds")
Idents(combined) <- 'ID_labs'
DefaultAssay(combined) <- 'ATAC'

#### Try with object Saher worked with ####
RNAobj <- readRDS("2510_HA107PBSLPSPBS-8wk_RNAdata.rds")
RNAobj
unique(RNAobj@meta.data$orig.ident)
rownames(RNAobj@meta.data)
# correct orig.idents and. names
RNAobj@meta.data$stimulation <- unlist(str_split(RNAobj@meta.data$orig.ident,"-"))[seq(3,length(RNAobj@meta.data$orig.ident)*4,4)]
RNAobj@meta.data$colonization <- unlist(str_split(RNAobj@meta.data$orig.ident,"-"))[seq(2,length(RNAobj@meta.data$orig.ident)*4,4)]
unique(RNAobj@meta.data$colonization)
unique(RNAobj@meta.data$orig.ident)
unique(combined@meta.data$orig.ident) %in% unique(RNAobj@meta.data$orig.ident)
# now barcodes
Cells(RNAobj)
RNAobj@meta.data$barcode <- str_sub(Cells(RNAobj),start = -18,end = -1)
# also barcodes for combined
combined@meta.data$barcode <- str_sub(Cells(combined),start = -18,end = -1)
# change names to match sample_barcode
#RNAobj <- RenameCells(RNAobj, new.names = paste(RNAobj@meta.data$orig.ident,RNAobj@meta.data$barcode,sep = "_")) 

## Check cells from combined are in RNAobj
Cells(combined)[!Cells(combined) %in% Cells(RNAobj)]
combined@meta.data$barcode[!combined@meta.data$barcode %in% RNAobj@meta.data$barcode]

# check how cells are distributed
notinRNA <- Cells(combined)[!Cells(combined) %in% Cells(RNAobj)]
notincomb <- Cells(RNAobj)[!Cells(RNAobj) %in% Cells(combined)]
table(combined@meta.data[notinRNA,]$orig.ident)/length(notinRNA)*100 #only around 500 cells total, distribution not important
table(RNAobj@meta.data[notincomb,]$orig.ident)/length(notincomb)*100 #almost 6000 "cells" likely removed contaminants/dublets/debris, almost evenly distributed among samples

table(combined@meta.data[notinRNA,]$ID_labs) #only around 500 cells total, distribution not important
combined$notinRNA <- "both"
combined@meta.data[notinRNA,]$notinRNA <- "not in RNA"
DimPlot(combined, group.by = "ID_labs",split.by = "notinRNA")
# looks evenly distributed, all good to continue without these

# subset RNA assay and save separate objects with unification of cells
cellsshared <- Cells(combined)[Cells(combined) %in% Cells(RNAobj)]
RNAobj <- subset(RNAobj, cells = cellsshared)
rnaint_assay <- RNAobj@assays$integrated
RNAobj <- subset(combined, cells = cellsshared)
RNAobj[["integrated"]] <- rnaint_assay

FeaturePlot(RNAobj, features = "rna_Ly6c2")+scale_colour_gradientn(colors = mycols)

saveRDS(RNAobj,"/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/2510_PBSHA107PBALPS_8wk_wRNA.rds")


