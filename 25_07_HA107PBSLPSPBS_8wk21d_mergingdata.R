#' R script for combining, clustering and annotating scATAC data - 8 BM samples
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

#### ---- Read in samples for integration ---- ####
# 8 wk samples
PBSPBS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-8wk/BM-PBS-PBS-8wk.rds")
HA107PBS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-8wk/BM-HA107-PBS-8wk.rds")
PBSLPS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-8wk/BM-PBS-LPS-8wk.rds")
HA107LPS8wk <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-8wk/BM-HA107-LPS-8wk.rds")
# 21d samples
PBSPBS21d <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-21d/BM-PBS-PBS-21d.rds")
HA107PBS21d <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-21d/BM-HA107-PBS-21d.rds")
PBSLPS21d <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-21d/BM-PBS-LPS-21d.rds")
HA107LPS21d <- readRDS("/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-21d/BM-HA107-LPS-21d.rds")

# only using cell names and meta data. Will be utilised to implement proper thresholds from each sample.
PBSPBS8wk.meta <- PBSPBS8wk@meta.data
HA107PBS8wk.meta <- HA107PBS8wk@meta.data
PBSLPS8wk.meta <- PBSLPS8wk@meta.data
HA107LPS8wk.meta <- HA107LPS8wk@meta.data
PBSPBS21d.meta <- PBSPBS21d@meta.data
HA107PBS21d.meta <- HA107PBS21d@meta.data
PBSLPS21d.meta <- PBSLPS21d@meta.data
HA107LPS21d.meta <- HA107LPS21d@meta.data

PBSPBS8wk <- Cells(PBSPBS8wk)
HA107PBS8wk <- Cells(HA107PBS8wk)
PBSLPS8wk <- Cells(PBSLPS8wk)
HA107LPS8wk <- Cells(HA107LPS8wk)
PBSPBS21d <- Cells(PBSPBS21d)
HA107PBS21d <- Cells(HA107PBS21d)
PBSLPS21d <- Cells(PBSLPS21d)
HA107LPS21d <- Cells(HA107LPS21d)

## Create common peak sets and subset the object assays based on this set
# read in peak sets
# 8 WK SAMPLES
peaks.PBSPBS8wk <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-8wk/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.HA107PBS8wk <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-8wk/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.PBSLPS8wk <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-8wk/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.HA107LPS8wk <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-8wk/peaks.bed",
  col.names = c("chr", "start", "end"))
# 21D SAMPLES
peaks.PBSPBS21d <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-21d/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.HA107PBS21d <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-21d/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.PBSLPS21d <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-21d/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.HA107LPS21d <- read.table(
  file = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-21d/peaks.bed",
  col.names = c("chr", "start", "end"))

# convert to genomic ranges
gr.PBSPBS8wk <- makeGRangesFromDataFrame(peaks.PBSPBS8wk)
gr.HA107PBS8wk <- makeGRangesFromDataFrame(peaks.HA107PBS8wk)
gr.PBSLPS8wk <- makeGRangesFromDataFrame(peaks.PBSLPS8wk)
gr.HA107LPS8wk <- makeGRangesFromDataFrame(peaks.HA107LPS8wk)
gr.PBSPBS21d <- makeGRangesFromDataFrame(peaks.PBSPBS21d)
gr.HA107PBS21d <- makeGRangesFromDataFrame(peaks.HA107PBS21d)
gr.PBSLPS21d <- makeGRangesFromDataFrame(peaks.PBSLPS21d)
gr.HA107LPS21d <- makeGRangesFromDataFrame(peaks.HA107LPS21d)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- GenomicRanges::reduce(x = c(gr.PBSPBS8wk,
                               gr.HA107PBS8wk,
                               gr.PBSLPS8wk,
                               gr.HA107LPS8wk,
                               gr.PBSPBS21d,
                               gr.HA107PBS21d,
                               gr.PBSLPS21d,
                               gr.HA107LPS21d))

# clean up
rm(peaks.PBSPBS8wk,peaks.PBSLPS8wk,peaks.HA107PBS8wk,peaks.HA107LPS8wk,peaks.PBSPBS21d,peaks.PBSLPS21d,peaks.HA107PBS21d,peaks.HA107LPS21d)
rm(gr.PBSLPS21d,gr.PBSLPS8wk,gr.PBSPBS21d,gr.PBSPBS8wk,gr.HA107LPS21d,gr.HA107LPS8wk,gr.HA107PBS21d,gr.HA107PBS8wk)

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# make samples from existing data with these peaks
frags.PBSPBS8wk <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-8wk/fragments.tsv.gz",
  cells = PBSPBS8wk)
frags.HA107PBS8wk <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-8wk/fragments.tsv.gz",
  cells = HA107PBS8wk)
frags.PBSLPS8wk <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-8wk/fragments.tsv.gz",
  cells = PBSLPS8wk)
frags.HA107LPS8wk <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-8wk/fragments.tsv.gz",
  cells = HA107LPS8wk)
frags.PBSPBS21d <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-PBS-21d/fragments.tsv.gz",
  cells = PBSPBS21d)
frags.HA107PBS21d <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-PBS-21d/fragments.tsv.gz",
  cells = HA107PBS21d)
frags.PBSLPS21d <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-PBS-LPS-21d/fragments.tsv.gz",
  cells = PBSLPS21d)
frags.HA107LPS21d <- CreateFragmentObject(
  path = "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/ATAC/BM-HA107-LPS-21d/fragments.tsv.gz",
  cells = HA107LPS21d)

counts.PBSPBS8wk <- FeatureMatrix(
  fragments = frags.PBSPBS8wk,
  features = combined.peaks,
  cells = PBSPBS8wk)
counts.HA107PBS8wk <- FeatureMatrix(
  fragments = frags.HA107PBS8wk,
  features = combined.peaks,
  cells = HA107PBS8wk)
counts.PBSLPS8wk <- FeatureMatrix(
  fragments = frags.PBSLPS8wk,
  features = combined.peaks,
  cells = PBSLPS8wk)
counts.HA107LPS8wk <- FeatureMatrix(
  fragments = frags.HA107LPS8wk,
  features = combined.peaks,
  cells = HA107LPS8wk)
counts.PBSPBS21d <- FeatureMatrix(
  fragments = frags.PBSPBS21d,
  features = combined.peaks,
  cells = PBSPBS21d)
counts.HA107PBS21d <- FeatureMatrix(
  fragments = frags.HA107PBS21d,
  features = combined.peaks,
  cells = HA107PBS21d)
counts.PBSLPS21d <- FeatureMatrix(
  fragments = frags.PBSLPS21d,
  features = combined.peaks,
  cells = PBSLPS21d)
counts.HA107LPS21d <- FeatureMatrix(
  fragments = frags.HA107LPS21d,
  features = combined.peaks,
  cells = HA107LPS21d)

PBSPBS8wk_assay <- CreateChromatinAssay(counts.PBSPBS8wk, fragments = frags.PBSPBS8wk)
PBSPBS8wk <- CreateSeuratObject(PBSPBS8wk_assay, assay = "ATAC", meta.data=PBSPBS8wk.meta)
HA107PBS8wk_assay <- CreateChromatinAssay(counts.HA107PBS8wk, fragments = frags.HA107PBS8wk)
HA107PBS8wk <- CreateSeuratObject(HA107PBS8wk_assay, assay = "ATAC", meta.data=HA107PBS8wk.meta)
PBSLPS8wk_assay <- CreateChromatinAssay(counts.PBSLPS8wk, fragments = frags.PBSLPS8wk)
PBSLPS8wk <- CreateSeuratObject(PBSLPS8wk_assay, assay = "ATAC", meta.data=PBSLPS8wk.meta)
HA107LPS8wk_assay <- CreateChromatinAssay(counts.HA107LPS8wk, fragments = frags.HA107LPS8wk)
HA107LPS8wk <- CreateSeuratObject(HA107LPS8wk_assay, assay = "ATAC", meta.data=HA107LPS8wk.meta)
PBSPBS21d_assay <- CreateChromatinAssay(counts.PBSPBS21d, fragments = frags.PBSPBS21d)
PBSPBS21d <- CreateSeuratObject(PBSPBS21d_assay, assay = "ATAC", meta.data=PBSPBS21d.meta)
HA107PBS21d_assay <- CreateChromatinAssay(counts.HA107PBS21d, fragments = frags.HA107PBS21d)
HA107PBS21d <- CreateSeuratObject(HA107PBS21d_assay, assay = "ATAC", meta.data=HA107PBS21d.meta)
PBSLPS21d_assay <- CreateChromatinAssay(counts.PBSLPS21d, fragments = frags.PBSLPS21d)
PBSLPS21d <- CreateSeuratObject(PBSLPS21d_assay, assay = "ATAC", meta.data=PBSLPS21d.meta)
HA107LPS21d_assay <- CreateChromatinAssay(counts.HA107LPS21d, fragments = frags.HA107LPS21d)
HA107LPS21d <- CreateSeuratObject(HA107LPS21d_assay, assay = "ATAC", meta.data=HA107LPS21d.meta)


# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = PBSPBS8wk,
  y = list(HA107PBS8wk,
           PBSLPS8wk,
           HA107LPS8wk,
           PBSPBS21d,
           HA107PBS21d,
           PBSLPS21d,
           HA107LPS21d),
  add.cell.ids = c("BM-PBS-PBS-8wk", "BM-HA107-PBS-8wk","BM-PBS-LPS-8wk","BM-HA107-LPS-8wk",
                   "BM-PBS-PBS-21d", "BM-HA107-PBS-21d","BM-PBS-LPS-21d","BM-HA107-LPS-21d"))

# check unique orig.idents exist
unique(combined@meta.data$orig.ident)
combined[["ATAC"]]
length(rownames(combined@assays$ATAC)) # 144016

# clean up
rm(PBSPBS8wk,HA107PBS8wk,PBSLPS8wk,HA107LPS8wk,PBSPBS21d,HA107PBS21d,PBSLPS21d,HA107LPS21d)
rm(counts.PBSPBS8wk,counts.PBSPBS21d,counts.PBSLPS8wk,counts.PBSLPS21d,counts.HA107PBS8wk,counts.HA107PBS21d,counts.HA107LPS8wk,counts.HA107LPS21d)
rm(frags.PBSPBS8wk,frags.PBSPBS21d,frags.PBSLPS8wk,frags.PBSLPS21d,frags.HA107PBS8wk,frags.HA107PBS21d,frags.HA107LPS8wk,frags.HA107LPS21d)
rm(PBSPBS8wk.meta,PBSPBS21d.meta,PBSLPS8wk.meta,PBSLPS21d.meta,HA107PBS8wk.meta,HA107PBS21d.meta,HA107LPS8wk.meta,HA107LPS21d.meta)
rm(PBSPBS8wk_assay,PBSPBS21d_assay,PBSLPS8wk_assay,PBSLPS21d_assay,HA107PBS8wk_assay,HA107PBS21d_assay,HA107LPS8wk_assay,HA107LPS21d_assay)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(combined) <- annotations

## Variable features, Dim reduction
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1)

CoveragePlot(
  object = combined,
  group.by = 'orig.ident',
  region = "chr14-99700000-99760000"
)

VlnPlot(object = combined, group.by = 'orig.ident', features = c("blacklist_fraction","nucleosome_signal"))


combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
res <- seq(0,1,0.1)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = res)

DimPlot(object = combined, label = TRUE,
        group.by = "ATAC_snn_res.0.2",
        split.by = "orig.ident", ncol = 4)
DimPlot(object = combined, label = TRUE,
        group.by = "ATAC_snn_res.0.2",
        split.by = "timepoint")
Idents(combined) <- 'ATAC_snn_res.0.2'

gene.activities <- GeneActivity(combined)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

#### Save object before continuing
saveRDS(combined, paste(dato,"PBSHA107PBALPS_8wk21d.rds"))
rm(combined.peaks)

rownames(combined@assays$RNA@data)[startsWith(rownames(combined@assays$RNA@data), "Ly6")]
#test
VlnPlot(combined, features = "Ly6c2", assay = "RNA",pt.size=0)


# Distribution
dist_df <- perc_function_samp("ATAC_snn_res.0.2", Cells(combined), combined,"orig.ident")
dist_df$samp <- factor(dist_df$samp, 
                          levels = c("BM-PBS-PBS-21d","BM-PBS-PBS-8wk","BM-HA107-PBS-21d","BM-HA107-PBS-8wk","BM-HA107-LPS-21d","BM-HA107-LPS-8wk","BM-PBS-LPS-21d","BM-PBS-LPS-8wk"))
ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black")+
  theme(axis.text.x = element_text(angle=90))

ggplot(dist_df, aes(x=cluster, y=percent, fill=samp))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme(axis.text.x = element_text(angle=90))


#### ---- Transfer labels from Sahers analysis ---- ####
old_meta <- read.csv(paste(projdir,"/Outputs/Clustering/25_07_118wkcombined_metadata_SaherAnalysis.csv",sep = ""), header = T, row.names = 1)
head(old_meta)
tail(old_meta)
head(combined@meta.data)
old_meta[old_meta$dataset=="TPBS_PBS",]$dataset <- "PBS_PBS"
old_meta$timepoint <- "8wk"
old_meta$colonization <- unlist(str_split(old_meta$dataset, "_"))[seq(1,length(rownames(old_meta))*2,2)]
old_meta$stimulation <- unlist(str_split(old_meta$dataset, "_"))[seq(2,length(rownames(old_meta))*2,2)]
unique(old_meta$stimulation)
unique(old_meta$colonization)
old_meta$orig.ident <- paste("BM",old_meta$colonization,old_meta$stimulation,"8wk", sep = "-")
unique(old_meta$orig.ident)
rownames(old_meta) <- paste(old_meta$orig.ident, old_meta$barcode, sep = "_")
length(rownames(old_meta))
rownames(old_meta)[rownames(old_meta) %in% rownames(combined@meta.data)]
rownames(old_meta)[1]
rownames(combined@meta.data)[1]

old_meta$barcode <- str_sub(old_meta$barcode, start = 1, end = -3)
old_meta$barcode <- paste(old_meta$barcode, "-1", sep = "")
oldone <- "ACGACAAAGGCGGGTA-1"
oldupt <- "TGCTGTTTCCGCCCAT-1"

combined@meta.data$barcode <- unlist(str_split(Cells(combined),"_"))[seq(2,length(Cells(combined))*2,2)]
combined@meta.data$barcode <- str_sub(combined@meta.data$barcode, start = 1, end = -3)
old_meta$barcode[old_meta$barcode %in% combined@meta.data$barcode]
## Somehow not working

# change back to working with peaks instead of gene activities
DefaultAssay(combined) <- 'ATAC'


## If running FindConservedMarkers specify assay, as will otherwise assume RNA assay
## nCount_peaks vs peak_region_fragments + essentially the same, nCount calc. by Seurat, 
## highly correlated though
## cor(combined$peak_region_fragments, combined$nCount_peaks) = 0.999848
DA_peaks_res.0.2_FC1 <- FindAllMarkers(
  object = combined,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 1)
  