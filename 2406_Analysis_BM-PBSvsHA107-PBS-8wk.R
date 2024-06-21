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
library(GenomicRanges) 

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

## Create common peak sets and subset the object assays based on this set
# read in peak sets
peaks.HA <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-HA107-PBS-8wk/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.PBS <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-PBS-PBS-8wk/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.HA <- makeGRangesFromDataFrame(peaks.HA)
gr.PBS <- makeGRangesFromDataFrame(peaks.PBS)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.HA, gr.PBS))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# make samples from existing data with these peaks
frags.PBS <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-PBS-PBS-8wk/fragments.tsv.gz",
  cells = Cells(BMPBSPBS8wk)
)
frags.HA <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-HA107-PBS-8wk/fragments.tsv.gz",
  cells = Cells(BMHAPBS8wk)
)

PBS.counts <- FeatureMatrix(
  fragments = frags.PBS,
  features = combined.peaks,
  cells = Cells(BMPBSPBS8wk)
)

HA.counts <- FeatureMatrix(
  fragments = frags.HA,
  features = combined.peaks,
  cells = Cells(BMHAPBS8wk)
)

PBS_assay <- CreateChromatinAssay(PBS.counts, fragments = frags.PBS)
BMPBSPBS8wk <- CreateSeuratObject(PBS_assay, assay = "ATAC", meta.data=BMPBSPBS8wk@meta.data)

HA_assay <- CreateChromatinAssay(HA.counts, fragments = frags.HA)
BMHAPBS8wk <- CreateSeuratObject(HA_assay, assay = "ATAC", meta.data=BMHAPBS8wk@meta.data)



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
DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1)

