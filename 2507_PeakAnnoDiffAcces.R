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
seu_obj <- readRDS("/Volumes/Promise RAID/Saher/FINAL MULTIOME COUNTS/BM-HA107-LPS-8wk/outs/21D-integrated_multiome.rds")

#### ----- Peak annotation ------- ####
# extract gene annotations from EnsDb
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"

annotations <- GetGRangesFromEnsDb(ensdb = edb)

# change to UCSC style since the data was mapped to hg19
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
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

distanceToNearest(seu_obj, subject = Annotation(seu_obj))
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

#### ---- Annotation in Diff. Acc. Tests ---- ####
acc_dat <- read.csv("/Volumes/Promise RAID/Saher/For Aline/BM 8 Weeks/Differential Accessibility Analysis with Genes - ATAC/DA cluster by cluster/DA_peaks_cluster2_HA107_LPS_vs_PBS_LPS_LPS.csv")
head(acc_dat)

## why are we doing this?
## show bar plot with dist. of region types as standard in signac
acc_stat_df <- as.data.frame(summary(as.factor(acc_dat$type)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  theme_classic()

## same as above but split into up vs down reg
acc_dat$regulation <- 'down'
acc_dat[acc_dat$avg_log2FC>0,]$regulation <- "up"
acc_dat$type_reg <- paste(acc_dat$type,acc_dat$regulation,sep = '_')
acc_stat_df <- as.data.frame(summary(as.factor(acc_dat$type_reg)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
acc_stat_df$regulation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(2,18,2)]
acc_stat_df$annotation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(1,18,2)]
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  facet_wrap(.~regulation)+
  theme_classic()

#### Now add the new annotations
peakAnno.edb <- annotatePeak(seu_obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)

rownames(acc_dat) <- acc_dat$peak
open_up <- rownames(acc_dat[acc_dat$avg_log2FC>0, ])
close_open_up <- ClosestFeature(seu_obj@assays$ATAC, regions = open_up, annotation = peakAnno.edb@anno)
head(close_open_up)

open_down <- rownames(acc_dat[acc_dat$avg_log2FC<0, ])
close_open_down <- ClosestFeature(seu_obj@assays$ATAC, regions = open_down, annotation = peakAnno.edb@anno)
head(close_open_down)

## Collapse intron after 1 and same for exon
# introns
close_open_up$annotation_upd <- close_open_up$annotation
for (ann in close_open_up[startsWith(close_open_up$annotation, "Intron"),]$annotation){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    close_open_up[close_open_up$annotation==ann,]$annotation_upd <- "1st Intron"
  }
  else {
    close_open_up[close_open_up$annotation==ann,]$annotation_upd <- "Other Intron"
  }
}
# exons
for (ann in close_open_up[startsWith(close_open_up$annotation, "Exon"),]$annotation){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    close_open_up[close_open_up$annotation==ann,]$annotation_upd <- "1st Exon"
  }
  else {
    close_open_up[close_open_up$annotation==ann,]$annotation_upd <- "Other Exon"
  }
}
close_open <- close_open_up$annotation_upd
close_open <- cbind(close_open, rep("up", times = length(close_open)))
colnames(close_open) <- c("annotation","regulation")

# now for the downreg/less accessible regions
# introns
close_open_down$annotation_upd <- close_open_down$annotation
for (ann in close_open_down[startsWith(close_open_down$annotation, "Intron"),]$annotation){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    close_open_down[close_open_down$annotation==ann,]$annotation_upd <- "1st Intron"
  }
  else {
    close_open_down[close_open_down$annotation==ann,]$annotation_upd <- "Other Intron"
  }
}
# exons
for (ann in close_open_down[startsWith(close_open_down$annotation, "Exon"),]$annotation){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    close_open_down[close_open_down$annotation==ann,]$annotation_upd <- "1st Exon"
  }
  else {
    close_open_down[close_open_down$annotation==ann,]$annotation_upd <- "Other Exon"
  }
}

close_open2 <- close_open_down$annotation_upd
close_open2 <- cbind(close_open2, rep("down", times = length(close_open2)))
colnames(close_open2) <- c("annotation","regulation")
close_open <- rbind(close_open, close_open2)
close_open <- as.data.frame(close_open)
head(close_open)
tail(close_open)

### Plot the updated annotations in bar plots for the same DAR
acc_stat_df <- as.data.frame(summary(as.factor(close_open$annotation)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

## same as above but split into up vs down reg
acc_dat <- close_open
acc_dat$type_reg <- paste(acc_dat$annotation,acc_dat$regulation,sep = '_')
acc_stat_df <- as.data.frame(summary(as.factor(acc_dat$type_reg)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
acc_stat_df$regulation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(2,32,2)]
acc_stat_df$annotation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(1,32,2)]
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  facet_wrap(.~regulation)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

