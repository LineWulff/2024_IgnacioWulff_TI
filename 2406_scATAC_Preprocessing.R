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
# Update sample between each
samp <- "LP-PBS-PBS-8wk"

#### ---- Creating object from downloaded data ---- ####
## # extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
## # change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# GCRm38 should be same as mm10, mm10 is the UCSC format - but same in terms of sequences, changing name, since ref gen. was mm10
genome(annotations) <- "mm10"

## Sample object creation
samp_dir <- paste(proj_data_dir, samp, sep = "")
counts <- Read10X_h5(filename = paste(samp_dir,"filtered_peak_bc_matrix.h5",sep = "/"))
motifs <- Read10X_h5(filename = paste(samp_dir,"filtered_tf_bc_matrix.h5",sep = "/"))

metadata <- read.csv(
  file = paste(samp_dir,"singlecell.csv",sep="/"),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste(samp_dir,'fragments.tsv.gz',sep="/"),
  min.cells = 10,
  min.features = 200
)

sampleobj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  motifs = motifs,
  meta.data = metadata
)


sampleobj$blacklist_fraction <- FractionCountsInRegion(
  object = sampleobj, 
  assay = 'peaks',
  regions = blacklist_mm10
)

## Add meta data and gene annotations
Annotation(sampleobj) <- annotations
sampleobj@meta.data$orig.ident <- samp
sampleobj@meta.data$tissue <- unlist(str_split(sampleobj@meta.data$orig.ident,"-"))[seq(1,nrow(sampleobj@meta.data)*4,4)]
sampleobj@meta.data$colonization <- unlist(str_split(sampleobj@meta.data$orig.ident,"-"))[seq(2,nrow(sampleobj@meta.data)*4,4)]
sampleobj@meta.data$stimulation <- unlist(str_split(sampleobj@meta.data$orig.ident,"-"))[seq(3,nrow(sampleobj@meta.data)*4,4)]
sampleobj@meta.data$timepoint <- unlist(str_split(sampleobj@meta.data$orig.ident,"-"))[seq(4,nrow(sampleobj@meta.data)*4,4)]

## after last sample remove the unnecesseary objects to save env. space
rm(chrom_assay, counts, metadata, annotations)

## Looking at the files at a glance
sampleobj[['peaks']]
granges(sampleobj)
head(sampleobj@assays$peaks@data)
sampleobj@assays$peaks@annotation
head(sampleobj@meta.data)


#### ---- Sample QC ---- ####
## Set directory to individual sample QC folder to save all info on QC and threshold subsetting
setwd(paste0(proj_data_dir,samp,"/QC"))

## Five measures:
## Nucleosome banding pattern - calculate, saved as nucleosome_signal
## Transcription starting site - calculate, saved TSS.Enrichment
## Total number of fragments in peaks - precalc. from CR,
## Fraction of fragments in peaks - precalc. from CR, 
## Ratio reads in generic blacklist regions - precalc from CR, 

# compute nucleosome signal score per cell
# ratio of reads per cell from  mononucleosome (147-294 bp) to nucleosome free (<147 bp)
sampleobj <- NucleosomeSignal(object = sampleobj)

# compute TSS enrichment score per cell
sampleobj <- TSSEnrichment(object = sampleobj, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
sampleobj$pct_reads_in_peaks <- sampleobj$peak_region_fragments / sampleobj$passed_filters * 100

pdf(paste(dato,samp,"scATAC_QC_VlnPlots.pdf", sep="_"),height = 8, width = 8)
VlnPlot(
  object = sampleobj,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 3
)
dev.off()

pdf(paste(dato,samp,"scATAC_QC_DensityPlot_CountPeaksvsTSSenrich.pdf", sep="_"),height = 4, width = 5)
DensityScatter(sampleobj, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf(paste(dato,samp,"scATAC_QC_Histo_CountPeaks.pdf", sep="_"),height = 4, width = 5)
ggplot(sampleobj@meta.data, aes(x = nCount_peaks))+
  geom_histogram(fill="grey",bins = 100)+
  geom_vline(xintercept = c(3000,26000), colour = "red")
dev.off()

pdf(paste(dato,samp,"scATAC_QC_DotPlot_CountPeakvsNucleosomeSign_ColPctReadsPeaks.pdf", sep="_"),height = 4, width = 5.5)
ggplot(sampleobj@meta.data, aes(x = nCount_peaks, y = nucleosome_signal, colour = pct_reads_in_peaks))+
  geom_point_rast()+scale_color_viridis_c()+
  geom_vline(xintercept = c(3000,26000))
dev.off()

pdf(paste(dato,samp,"scATAC_QC_DotPlot_CountPeakvsNucleosomeSign_ColBlacklist.pdf", sep="_"),height = 4, width = 5.5)
ggplot(sampleobj@meta.data, aes(x = nCount_peaks, y = nucleosome_signal, colour = blacklist_fraction))+
  geom_point_rast()+scale_color_viridis_c()+
  geom_vline(xintercept = c(3000,26000))
dev.off()

## Nucleosome banding patterns
# grouping cells based on their mononuc/nfr ration, here 2:1 - 2x mononucleosome bound to nfr 
sampleobj$nucleosome_group <- ifelse(sampleobj$nucleosome_signal > 2, 'NS > 2', 'NS < 2') 
# plotted with fragment histo 
length(sampleobj$nucleosome_group[sampleobj$nucleosome_group=='NS > 2'])

pdf(paste(dato,samp,"scATAC_QC_DotPlot_FragmentDist_SplitByNucleosomesign.pdf", sep="_"),height = 4, width = 8)
FragmentHistogram(sampleobj, group.by = 'nucleosome_group', region = "chr1-1-20000000")+geom_vline(xintercept = 147)
dev.off()

## TSS enrichment
sampleobj$high.tss <- ifelse(sampleobj$TSS.enrichment > 3, 'High', 'Low')

pdf(paste(dato,samp,"scATAC_QC_DotPlot_TSSenrich_SplitHighvsLow.pdf", sep="_"),height = 4, width = 8)
TSSPlot(sampleobj, group.by = 'high.tss') + NoLegend() + geom_hline(yintercept = 3)
dev.off()
TSSPlot(sampleobj, group.by = 'nucleosome_group') + NoLegend()


#### ---- Subset and save thresholds ---- ####
nCount_low = 3000
nCount_high = 26000
perc_readspeaks = 15
nuc_sign = 2
blacklist_th = 0.05
TSS.enrich = 3

npre <- length(Cells(sampleobj))

sampleobj <- subset(
  x = sampleobj,
  subset = nCount_peaks > nCount_low &
    nCount_peaks < nCount_high &
    pct_reads_in_peaks > perc_readspeaks &
    blacklist_fraction < blacklist_th &
    nucleosome_signal < nuc_sign &
    TSS.enrichment > TSS.enrich
)

npost <- length(Cells(sampleobj))
remperc <- (npre-npost)/npost*100
remperc

sampleobj

ThresFile <- file(paste(dato,"SubsetThresholds",samp,".txt",sep="_"))
writeLines(c(paste("Pre subsetting there were",npre,"cells."),
             paste("Removing ~",remperc,"% of cells during single cell QC."),
             paste("Post subsetting there are:",npost,"cells"),
             "Thresholds were set to:",
             paste("Peaks per cell:",nCount_low,"to",nCount_high),
             paste("Percentage of reads in peaks: >",perc_readspeaks),
             paste("Nucleosome signal (Mononuc/NFR ratio based): <",nuc_sign),
             paste("Mean TSS enrichment was set to: >",TSS.enrich),
             paste("Blacklist fraction: <",blacklist_th,"%")), 
           ThresFile)
close(ThresFile)

#### ---- Adding gene activity matrix based on open chromatin regions ---- ####
gene.activities <- GeneActivity(sampleobj)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
sampleobj[['RNA']] <- CreateAssayObject(counts = gene.activities)
sampleobj <- NormalizeData(
  object = sampleobj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sampleobj$nCount_RNA)
)

#### ---- save object ---- ####
## to sample folder
saveRDS(sampleobj, file = paste0(proj_data_dir,samp,"/",samp,".rds"))
