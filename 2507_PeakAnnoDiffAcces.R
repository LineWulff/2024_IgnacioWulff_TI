#' R script for annotating regions of ATAC comp results
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


##### ----- Read in data ------ ####
acc_dat <- read.csv("/Volumes/Promise RAID/Saher/For Aline/BM 8 Weeks/Differential Accessibility Analysis with Genes - ATAC/DA cluster by cluster/DA_peaks_cluster2_HA107_LPS_vs_PBS_LPS_LPS.csv")
seu_obj <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/21D-integrated_multiome.rds")

#### ----- Peak annotation ------- ####
# extract gene annotations from EnsDb
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"

annotations <- GetGRangesFromEnsDb(ensdb = edb)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(seu_obj[["ATAC"]]) <- annotations

#edb <- TxDb.Mmusculus.UCSC.mm10.knownGene
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"

peakAnno.edb <- annotatePeak(seu_obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
peakAnno.edb
plotAnnoPie(peakAnno.edb)
plotAnnoBar(peakAnno.edb)

annoStat.edb <- as.data.frame(peakAnno.edb@annoStat)
colnames(annoStat.edb)[2] <- "combined"

#### for the seperate samples
## calulate each and add statistic to common dataframe for plotting
ha107_lps <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/HA107_LPS_pbmc_multiome.rds")
pbs_lps   <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM PBS LPS 8 8wk/outs/PBS_LPS_pbmc_multiome.rds")
ha107_pbs <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-PBS-8wk/outs/HA107_PBS_pbmc_multiome.rds")
pbs_pbs  <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-PBS-PBS-8wk/outs/TPBS_PBS_pbmc_multiome.rds")

## ha107_lps
peakAnno.edb <- annotatePeak(ha107_lps@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
peakAnno.edb
plotAnnoPie(peakAnno.edb)
annoStat.edb <- cbind(annoStat.edb, ha107_lps = peakAnno.edb@annoStat[,2])

## ha107_pbs
peakAnno.edb <- annotatePeak(ha107_pbs@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
peakAnno.edb
plotAnnoPie(peakAnno.edb)
annoStat.edb <- cbind(annoStat.edb, ha107_pbs = peakAnno.edb@annoStat[,2])

## pbs_lps
peakAnno.edb <- annotatePeak(pbs_lps@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
peakAnno.edb
plotAnnoPie(peakAnno.edb)
annoStat.edb <- cbind(annoStat.edb, pbs_lps = peakAnno.edb@annoStat[,2])

## pbs_pbs
peakAnno.edb <- annotatePeak(pbs_pbs@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
peakAnno.edb
plotAnnoPie(peakAnno.edb)
annoStat.edb <- cbind(annoStat.edb, pbs_pbs = peakAnno.edb@annoStat[,2])

## transform dataframe and plot
annoStat.edb <- gather(annoStat.edb, key='data_orig', value='freq', combined:pbs_pbs)
ggplot(annoStat.edb, aes(x=Feature, y=freq,fill=data_orig))+
  geom_bar(stat = "identity", colour="black", position = 'dodge')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))


