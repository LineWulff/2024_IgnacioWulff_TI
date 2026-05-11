#' R script for finding TFs in genes associated to DARs
#' Author: Line Wulff
#' Date (created): 26-05-11
rm(list=ls())

#### ---- Initiate libraries ---- ####
library(readxl)
library(ggplot2)
library(ggrastr)
library(dplyr)
library(stringr)
library(tidyr)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(GenomicRanges)
library(scales)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10) # DNA sequence - BSgenome for coordinates/sequence
library(ChIPseeker)
library(patchwork)
library(ggseqlogo)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(viridis)

#### Read in seurat data and DAR outputs ####
# Seurat
obj <- readRDS("260212_BM-HA107-PBS-LPS-PBS-8wk_impRNAseq_wlabs.rds")
Idents(obj) <- "ID_labs_ext"
obj_sub <- subset(obj, idents = c("Q1","Q2","Q3"))


# region associated gene annotations
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),TxDb = edb)

#### --- common DARs to use ---- ####
# make df of all for best testing
DARs_tot <- data.frame()
# PBS-PBS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSPBS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSPBS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj_sub, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "PBS-PBS"

ass_genanno <- ClosestFeature(obj_sub, regions = regs, annotation = peakAnno.edb@anno)
ass_genanno <- ass_genanno$annotation
DARs_sub <- data.frame(row.names = regs, genes = hits, orig.ident = rep(org,length(hits)), annotation = ass_genanno)
DARs_tot <- rbind(DARs_tot, DARs_sub)

# HA107-PBS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107PBS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107PBS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj_sub, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "HA107-PBS"

ass_genanno <- ClosestFeature(obj_sub, regions = regs, annotation = peakAnno.edb@anno)
ass_genanno <- ass_genanno$annotation
DARs_sub <- data.frame(row.names = regs, genes = hits, orig.ident = rep(org,length(hits)), annotation = ass_genanno)
DARs_tot <- rbind(DARs_tot, DARs_sub)

# PBS-LPS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSLPS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsPBSLPS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj_sub, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "PBS-LPS"

ass_genanno <- ClosestFeature(obj_sub, regions = regs, annotation = peakAnno.edb@anno)
ass_genanno <- ass_genanno$annotation
DARs_sub <- data.frame(row.names = regs, genes = hits, orig.ident = rep(org,length(hits)), annotation = ass_genanno)
DARs_tot <- rbind(DARs_tot, DARs_sub)

# HA107-LPS
regs <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107LPS_Q1-Q3_intersectedregions.txt")
gens <- readLines("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/26_03_30BM_8wk_IDlabsHA107LPS_Q1-Q3_intersectedgenes.txt")
hits <- ClosestFeature(obj_sub, regions = regs); hits <- hits$gene_name
hits[!hits %in% gens] #should be none 
length(regs)==length(hits) # should be TRUE
org <- "HA107-LPS"

ass_genanno <- ClosestFeature(obj_sub, regions = regs, annotation = peakAnno.edb@anno)
ass_genanno <- ass_genanno$annotation
DARs_sub <- data.frame(row.names = regs, genes = hits, orig.ident = rep(org,length(hits)), annotation = ass_genanno)
DARs_tot <- rbind(DARs_tot, DARs_sub)

# Check total DAR df
nrow(DARs_tot) # 141
unique(DARs_tot$orig.ident)
length(unique(DARs_tot$genes)) # 133

DARs_tot$colonization <- unlist(str_split(DARs_tot$orig.ident,"-"))[seq(1,length(DARs_tot$orig.ident)*2,2)]
DARs_tot$stimulation <- unlist(str_split(DARs_tot$orig.ident,"-"))[seq(2,length(DARs_tot$orig.ident)*2,2)]

#### ---- Read in mouse TF DB ---- ####
tf_df <- read.delim("/Users/linewulff/Documents/work/SetUp/Databases/Mus_musculus_TF.txt")
head(tf_df)
tf_df <- unique(tf_df$Symbol)
tf_df

head(DARs_tot)

DARs_tf <- sort(tf_df[tf_df %in% DARs_tot$genes])

#### ---- available for motif analysis standard ---- ####
DARs_tf

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Get available seqlevels from your BSgenome
available_seqlevels <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
regions <- obj_sub@assays$ATAC@ranges

# Filter regions
regions_filtered <- regions[seqnames(regions) %in% available_seqlevels]
# Then use regions_filtered in AddMotif

# add motif information
obj_sub <- AddMotifs(
  object = obj_sub,
  regions = regions_filtered,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

obj_motifs <- unlist(obj_sub@assays$ATAC@motifs@motif.names)
DARs_tf[str_to_upper(DARs_tf) %in% obj_motifs] # 13 there
DARs_tf[!str_to_upper(DARs_tf) %in% obj_motifs] # 6 not there

## Of these only 2 where available q. weighted matrices on  https://cisbp.ccbr.utoronto.ca/
## Cic and Fosb analysis was included in script 2604_commonDAR_motifenrich.R
