#' R script for preprocessing samples with scATACseq data
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

#### ---- variables used throughout script ---- ####
rm(list = ls())
projdir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj_data_dir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/"
## sample combination
project <- "BM-PBSvsHA107-PBS-8wk"

#### ---- Read in samples for integration ---- ####
BMPBSPBS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-PBS-PBS-8wk/BM-PBS-PBS-8wk.rds")
BMHAPBS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-HA107-PBS-8wk/BM-HA107-PBS-8wk.rds")

## Identify common peaks
rownames(BMPBSPBS8wk@assays$peaks)[rownames(BMPBSPBS8wk@assays$peaks) %in% rownames(BMHAPBS8wk@assays$peaks)] # only 178?

# try anyway
# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = BMPBSPBS8wk,
  y = BMHAPBS8wk,
  add.cell.ids = c("BM-PBS-PBS-8wk", "BM-HA107-PBS-8wk")
)

length(rownames(BMPBSPBS8wk@assays$peaks))+length(rownames(BMHAPBS8wk@assays$peaks)) # 254136
length(rownames(combined@assays$peaks)) # 144016

length(rownames(BMPBSPBS8wk@assays$RNA))+length(rownames(BMHAPBS8wk@assays$RNA)) # 43616
length(rownames(combined@assays$RNA)) # 21808

## Variable features, Dim reduction
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

