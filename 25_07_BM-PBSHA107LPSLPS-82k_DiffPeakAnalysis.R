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

sheets <- c("cluster_Ly6c lo monocytes","cluster_Neutrophils","cluster_Dendritic cells","cluster_LSK","cluster_Ly6c hi monocytes","cluster_NK cells")

#### ---- Running DAR between conditions --- ####
# 6 loops
# 3 per timepoint
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(combined@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
# color control for annotations
annotations_gen <- peakAnno.edb@anno$annotation
for (ann in annotations_gen[startsWith(annotations_gen, "Intron")]){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    annotations_gen[annotations_gen==ann]<- "1st Intron"}
  else {
    annotations_gen[annotations_gen==ann] <- "Other Intron"
  }}
# exons
for (ann in annotations_gen[startsWith(annotations_gen, "Exon")]){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    annotations_gen[annotations_gen==ann] <- "1st Exon"}
  else {
    annotations_gen[annotations_gen==ann] <- "Other Exon"
  }}
length(annotations_gen) # should be 208k
# make a color chart
ann_col_val <- hue_pal()(length(unique(annotations_gen)))
names(ann_col_val) <- unique(annotations_gen)
ann_col_val
show_col(ann_col_val)

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

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_PBSvsPBS_8wk/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx", sep = ""),
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
### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_LPSvsLPS_8wk/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx",
         colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

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

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/PBSvsPBS_PBSvsLPS_8wk/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

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
### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_PBSvsPBS_21d/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSHA107PBSPBS21d.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107PBSvsPBSPBS21d.xlsx", sep = ""),
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

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_LPSvsLPS_21d/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSHA107LPSLPS21d.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107PBSvsPBLPSLPS21d.xlsx", sep = ""),
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

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/PBSvsPBS_PBSvsLPS_8wk/DiffPeaksConditions25_07_16_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }
  
  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_IDlabs_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)


#### ---- test for ditribution plots ---- ####
## Read in and correct csv file


for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(combined@assays$ATAC, cells = rownames(combined@meta.data[combined@meta.data$ID_labs==clus,])),
                                  regions = open_regs, annotation = peakAnno.edb@anno)
  ## Collapse annotations for exon and introns
  closest_genes$annotation_upd <- closest_genes$annotation
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Intron"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Intron"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Intron"
    }
  }
  # exons
  for (ann in closest_genes[startsWith(closest_genes$annotation, "Exon"),]$annotation){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "1st Exon"
    }
    else {
      closest_genes[closest_genes$annotation==ann,]$annotation_upd <- "Other Exon"
    }
  }

  
  DA_peaks_conditions[[clus]]$annotation <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$annotation <- closest_genes$annotation_upd
  # print head to check output
  print(head(DA_peaks_conditions[[clus]]))}


clus <- "cluster_Ly6c hi monocytes"
clus <- "cluster_LSK"
clus <- "cluster_Ly6c lo monocytes"

## Bar and volcano plots highlighting annotations
acc_stat_df <- as.data.frame(summary(as.factor(DA_peaks_conditions[[clus]]$annotation)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = ann_col_val)+
  theme(axis.text.x = element_text(angle = 90))

## same as above but split into up vs down reg
acc_dat <- DA_peaks_conditions[[clus]]
acc_dat$type_reg <- paste(acc_dat$annotation,acc_dat$sign,sep = '_')
acc_stat_df <- as.data.frame(summary(as.factor(acc_dat$type_reg)))
colnames(acc_stat_df) <- 'amount'
acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
acc_stat_df$sign <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(2,length(acc_stat_df$amount)*2,2)]
acc_stat_df$annotation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(1,length(acc_stat_df$amount)*2,2)]
ggplot(acc_stat_df, aes(x=annotation, y=amount, fill=annotation))+
  geom_bar(stat="identity")+
  facet_wrap(.~sign)+
  scale_fill_manual(values = ann_col_val)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

## freq of DAR annotations
DAsum <- rowsum(acc_stat_df$amount, group=acc_stat_df$sign)
acc_stat_df$totDEG <- 0; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[1],]$totDEG <- DAsum[1]; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[2],]$totDEG <- DAsum[2];
acc_stat_df$freq <- 0; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[1],]$freq <- DAsum[1]; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[2],]$freq <- DAsum[2];
acc_stat_df$freq <- acc_stat_df$amount/acc_stat_df$freq
acc_stat_df$annotation <- factor(acc_stat_df$annotation ,
                                 levels = c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                                            "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)"))
pdf(paste(paste0(outdir,dato,"_HA107PBS_PBSPBS_21d_annotatedDARs_dist_Ly6clomono.pdf")),height = 2, width = 6)
ggplot(acc_stat_df, aes(x=sign, y=freq, fill=annotation))+
  geom_bar(stat="identity", colour = "black")+
  scale_y_continuous(labels = scales::percent)+
  geom_text(data=acc_stat_df[acc_stat_df$annotation=="Promoter (<=1kb)",], aes(x=sign, label=totDEG, y=1.1, fill=NULL))+
  scale_fill_manual(values = ann_col_val)+
  labs(y="",x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  coord_flip()
dev.off()

## volcano plot
ggplot(DA_peaks_conditions[[clus]], aes(x = avg_log2FC, y = -log10(p_val_adj), colour = annotation))+
  geom_point_rast()+
  geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
  geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
  geom_vline(xintercept = c(0))+ #0
  scale_color_manual(values = ann_col_val)+
  theme_classic()+
  ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
  guides(colour=guide_legend(title="OR Annotation"))

gene_int <- c("Btg2", "Neat1", "Fos",   "Tex14", "Zfp36", "Fos")
  #c("Fos","Cxcr4","Ifnar1","Zfp36","Sirpa","Litaf" ) #"E2f2","Ebf2","Spi1"


ggplot(DA_peaks_conditions[[clus]], aes(x = avg_log2FC, y = -log10(p_val_adj), colour = annotation))+
  geom_point_rast()+
  geom_label_repel(data = DA_peaks_conditions[[clus]][DA_peaks_conditions[[clus]]$gene_name %in% gene_int,], 
                   aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene_name, fill = annotation), 
                   color = "black", nudge_y = 75, show.legend=F, max.overlaps = 15)+
  geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
  geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
  geom_vline(xintercept = c(0))+ #0
  scale_color_manual(values = ann_col_val)+scale_fill_manual(values = ann_col_val)+
  theme_classic()+
  ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
  guides(colour=guide_legend(title="Significance"))

#### ---- test for venn diagrams ---- ####
Ly6himono_df <- list()
Ly6clomono_df <- list()
LSK_df <- list()

## HA107-PBS vs PBS-PBS 8wk
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_PBSvsPBS_8wk/DiffPeaksConditions25_07_22_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107PBSvsPBSPBS8wk.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}
Ly6himono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c hi monocytes"]]
Ly6clomono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c lo monocytes"]]
LSK_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_LSK"]]

## HA107-PBS vs PBS-PBS 21d
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_PBSvsPBS_21d/DiffPeaksConditions25_07_22_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107PBSvsPBSPBS21d.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}
Ly6himono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c hi monocytes"]]
Ly6clomono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c lo monocytes"]]
LSK_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_LSK"]]

## HA107-LPS vs PBS-LPS 21d
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_LPSvsLPS_21d/DiffPeaksConditions25_07_22_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107PBSvsPBLPSLPS21d.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}
Ly6himono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c hi monocytes"]]
Ly6clomono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c lo monocytes"]]
LSK_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_LSK"]]

## HA107-LPS vs PBS-LPS 8wk
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets[1:5]))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/DiffPeaksConditions/HA107vsPBS_LPSvsLPS_8wk/DiffPeaksConditions25_07_22_LSKMonoNeu_IDlabs_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx",
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}
Ly6himono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c hi monocytes"]]
Ly6clomono_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_Ly6c lo monocytes"]]
LSK_df[[paste(unique(DA_peaks_conditions[[1]]$sign)[1],unique(DA_peaks_conditions[[1]]$sign)[2], sep = "_vs_")]] <- DA_peaks_conditions[["cluster_LSK"]]


### test - both up and down 
## more accessible in HA107, less in PBS
# HA107-PBS vs PBS-PBS, shared 8wk and 21d
# 1
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-HA"),]$gene_name,
                      Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-HA"),]$gene_name))

# HA107-PBS vs PBS-PBS, shared 8wk and 21d and shared with HA107-LPS vs PBS-LPS 8wk
# 1
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-HA"),]$gene_name,
                      Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-HA"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-HA"),]$gene_name))
# HA107-PBS vs PBS-PBS and HA107-LPS vs PBS-LPS shared 8wk in HA107 samples
# 12
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-HA"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-HA"),]$gene_name))

# HA107-PBS vs PBS-PBS 21d shared with HA107-LPS vs PBS-LPS 8wk
# 12
Reduce(intersect,list(Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-HA"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-HA"),]$gene_name))


## less accessible in HA107, more in PBS
# HA107-PBS vs PBS-PBS, shared 8wk and 21d
# 0
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-PB"),]$gene_name,
                      Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-PB"),]$gene_name))

# HA107-PBS vs PBS-PBS, shared 8wk and 21d and shared with HA107-LPS vs PBS-LPS 8wk
# 0
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-PB"),]$gene_name,
                      Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-PB"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-PB"),]$gene_name))

# HA107-PBS vs PBS-PBS and HA107-LPS vs PBS-LPS shared 8wk in HA107 samples
# 16
Reduce(intersect,list(Ly6himono_df[[1]][startsWith(Ly6himono_df[[1]]$sign,"BM-PB"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-PB"),]$gene_name))

# HA107-PBS vs PBS-PBS 21d shared with HA107-LPS vs PBS-LPS 8wk
# 2
Reduce(intersect,list(Ly6himono_df[[2]][startsWith(Ly6himono_df[[2]]$sign,"BM-PB"),]$gene_name,
                      Ly6himono_df[[4]][startsWith(Ly6himono_df[[4]]$sign,"BM-PB"),]$gene_name))

#### ---- Plotting genomic ranges with peaks --- ####
## OBS!! CONNECT TO RAID SERVER, otherwise error will prompt because of missing path to index of fragments
#regions_highlight <- subsetByOverlaps(StringToGRanges(int), LookupGeneCoords(combined, "Klf13"))

CoveragePlot(
  object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ID_labs=="Ly6c hi monocytes" &
                                                                  combined@meta.data$timepoint=='8wk',])),
  region = "Fos",
  #region.highlight = regions_highlight,
  extend.upstream = 3000,
  extend.downstream = 1000
)

frags <- Fragments(combined)
frags[[1]]@path
