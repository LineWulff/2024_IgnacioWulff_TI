#' R script for DAR between conditions in seperate clusters
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
project <- "BM-PBSvsHA107-PBSvsLPS-21dvs8wk"

#### --- Read in data ---- ####
combined <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/25_07_16_PBSHA107PBALPS_8wk21d_clean_v2.rds")
Idents(combined) <- 'ID_labs'
DefaultAssay(combined) <- 'ATAC'


#### ---- Running DAR between conditions --- ####
# 6 loops
# 3 per timepoint
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(combined@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions"
  
#### ---- PBS-PBS-8wk vs HA107-PBS-8wk ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-8wk', 
    ident.2 = 'BM-HA107-PBS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-PBS-8wk"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-HA107-PBS-8wk"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-PBS-8wk'="#F8766D",'BM-HA107-PBS-8wk'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSvsHA107.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                         DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond.xlsx", sep = ""),
           rowNames = T)

#### ---- PBS-LPS-8wk vs HA107-LPS-8wk ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-LPS-8wk', 
    ident.2 = 'BM-HA107-LPS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-LPS-8wk"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-HA107-LPS-8wk"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-LPS-8wk'="#F8766D",'BM-HA107-LPS-8wk'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSLPSvsHA107LPS.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)

#### ---- PBS-PBS-8wk vs PBS-LPS-8wk ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-8wk', 
    ident.2 = 'BM-PBS-LPS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-PBS-8wk"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-LPS-8wk"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-PBS-8wk'="#F8766D",'BM-PBS-LPS-8wk'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSPBSvsPBSLPS.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)

#### ---- PBS-PBS-21d vs HA107-PBS-21d ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-21d', 
    ident.2 = 'BM-HA107-PBS-21d',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-PBS-21d"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-HA107-PBS-21d"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-PBS-21d'="#F8766D",'BM-HA107-PBS-21d'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSvsHA107_21d.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSHA107PBSPBS21d.xlsx", sep = ""),
           rowNames = T)

#### ---- PBS-LPS-21d vs HA107-LPS-21d ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-LPS-21d', 
    ident.2 = 'BM-HA107-LPS-21d',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-LPS-21d"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-HA107-LPS-21d"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-LPS-21d'="#F8766D",'BM-HA107-LPS-21d'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSvsHA107LPS_21d.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSHA107LPSLPS21d.xlsx", sep = ""),
           rowNames = T)

#### ---- PBS-PBS-21d vs PBS-LPS-21d ---- ####
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"

for (clus in unique(combined@meta.data$ID_labs)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-21d', 
    ident.2 = 'BM-PBS-LPS-21d',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  DA_peaks_IDlabs$gene_name <-NA
  DA_peaks_IDlabs[rownames(DA_peaks_IDlabs) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  # beacuse not all have significant values, should include option of none
  DA_peaks_IDlabs$sign <- "not sign."
  if (any(DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-PBS-21d"}
  if (any(DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05)){
    DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]$sign <- "BM-PBS-LPS-21d"}
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_IDlabs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-PBS-21d'="#F8766D",'BM-PBS-LPS-21d'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_cluster",clus,"PBSvsLPS_21d.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSPBSPBSLPS21d.xlsx", sep = ""),
           rowNames = T)


