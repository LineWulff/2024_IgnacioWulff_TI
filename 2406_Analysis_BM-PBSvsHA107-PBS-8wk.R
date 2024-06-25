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
library(colorRamp2)

#### ---- variables used throughout script ---- ####
rm(list = ls())
projdir <- getwd()
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

length(rownames(combined@assays$ATAC)) # 144016

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
        group.by = "ATAC_snn_res.0.4",
        split.by = "orig.ident")
Idents(combined) <- 'ATAC_snn_res.0.4'

gene.activities <- GeneActivity(combined)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)


rownames(combined@assays$RNA@data)[startsWith(rownames(combined@assays$RNA@data), "Ly6")]
#test
VlnPlot(combined, features = "Ly6g", assay = "RNA", split.by = "orig.ident")


# change back to working with peaks instead of gene activities
DefaultAssay(combined) <- 'ATAC'


## If running FindConservedMarkers specify assay, as will otherwise assume RNA assay
## nCount_peaks vs peak_region_fragments + essentially the same, nCount calc. by Seurat, 
## highly correlated though
## cor(combined$peak_region_fragments, combined$nCount_peaks) = 0.999848
DA_peaks_res.0.4_FC1 <- FindAllMarkers(
  object = combined,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 1
)

head(DA_peaks_res.0.4_FC1)

## add gene name to the regions to infer activity
open_regs <- rownames(DA_peaks_res.0.4_FC1)
closest_genes_FC1 <- ClosestFeature(combined, regions = open_regs)
head(closest_genes_FC1[,c("gene_name","query_region")])
DA_peaks_res.0.4_FC1$gene_name <- NA
DA_peaks_res.0.4_FC1[rownames(DA_peaks_res.0.4_FC1) %in% closest_genes_FC1$query_region,]$gene_name <- closest_genes_FC1$gene_name

head(DA_peaks_res.0.4_FC1)
write.csv(DA_peaks_res.0.4_FC1, file = paste("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/Clustering/",dato, project,"_res.0.4_PeaksPrClus_FC1.csv",sep = ""))


#### Cluster IDs #### 
rownames(combined@assays$RNA@data)[startsWith(rownames(combined@assays$RNA@data),"Cd34")]

# HSPC markers - see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969381/
FeaturePlot(combined, features = "rna_Cd34")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Gata2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Procr")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Kit")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Fgd5")+scale_colour_gradientn(colors = mycols)


# monocyte markers
FeaturePlot(combined, features = "rna_Ly6c2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly6c1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Lyz2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ccr2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Csf2ra")+scale_colour_gradientn(colors = mycols)

# DC progenitor markers
FeaturePlot(combined, features = "rna_Flt3")+scale_colour_gradientn(colors = mycols)
