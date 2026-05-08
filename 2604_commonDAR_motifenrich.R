#' R script for DAR motif to gene x gene corr.
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on hhttps://stuartlab.org/signac/articles/pbmc_vignette.html
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
source('2602_motifexpcorfunction.R')

#### Variables to use throughout ####
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/"
outdir2 <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/commonDAR_cov&vln/"

dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj <- "BM_8wk_IDlabs"; project <- proj
clus <- "Q2"
file_path <- paste0(outdir,"26_01_23_MonoTraj_Q1-Q3_FilteredOverview_Q2.xlsx")

# colours 
col_fun <- colorRamp2(
  c(-0.75,0, 0.75),
  rev(c('#d7191c','white','darkblue')))

#### Read in seurat data and DAR outputs ####
# Seurat
obj <- readRDS("260212_BM-HA107-PBS-LPS-PBS-8wk_impRNAseq_wlabs.rds")
Idents(obj) <- "ID_labs_ext"
obj_sub <- subset(obj, idents = c("Q1","Q2","Q3"))

# violin colors
obj_sub$colonization <- factor(obj_sub$colonization, levels=c("PBS","HA107"))
col_cols <- c("PBS"="#F8766D","HA107"="#00BFC4")
stim_cols <- c("PBS"="white","LPS"="lightgrey")
# GRanges object structure
GRobj <- granges(obj_sub)

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
  pfm = pfm_custom
)

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

#### ---- Enrichment analysis loop ---- ####
motif_mat <- obj_sub@assays$ATAC@motifs@data
tail(colnames(motif_mat))

for (org in unique(DARs_tot$orig.ident)){
  print(org)
  col <- unique(DARs_tot[DARs_tot$orig.ident==org,]$colonization)
  obj_sub1 <- subset(obj_sub, cells = rownames(obj_sub@meta.data[obj_sub@meta.data$colonization==col,]))
  enriched_motifs <- FindMotifs(
    object = obj_sub1,
    features = rownames(DARs_tot[DARs_tot$orig.ident==org,]))  # your foreground regions
  print("Pre and post padjust threshold:")
  print(dim(enriched_motifs)[1])
  enriched_motifs <- enriched_motifs[enriched_motifs$p.adjust<0.05,]
  print(dim(enriched_motifs)[1])
  
  if (dim(enriched_motifs)[1]>0){
  dar_motifs <- motif_mat[rownames(DARs_tot[DARs_tot$orig.ident==org,]),]
  # Rows = DARs, Columns = motifs, Values = motif score / presence
  
  # now subset the dar relevant motifs by the significant motifs
  dar_motifs <- dar_motifs[,enriched_motifs$motif]
  hits <- DARs_tot[DARs_tot$orig.ident==org,]$genes
  # stopifnot(
  #   nrow(dar_motifs) == length(ass_genenames),
  #   nrow(dar_motifs) == length(ass_genanno)
  # )
  
  
  # for visualization purposes
  #rownames(dar_motifs) <- DARs_tot[DARs_tot$orig.ident==org,]$genes
  
  #plottable data - binary
  mat <-  apply(as.matrix(dar_motifs), 2, as.numeric)
  rownames(mat) <- hits
  colnames(mat) <- enriched_motifs$motif.name
  
  out_mat <- motifexpcor(obj_sub, mat, enriched_motifs$motif.name, hits)
  
  ## Same but without NAs for each enroched motif
  # for (i in seq(1,dim(out_mat)[1])){
  #   mat_i <- out_mat[i,]
  #   mat_i <- na.omit(mat_i)
  #   hi <- length(mat_i)
  #   nami <- rownames(out_mat)[i]
  #   Ht_i <- Heatmap(mat_i,
  #                   rect_gp = gpar(col = "black", lwd = 0.5),
  #                   col = col_fun,
  #                   show_column_names = FALSE,
  #                   show_row_dend = FALSE,
  #                   column_title = nami)
  #   pdf(paste0(outdir2,dato,"_",org,"_",nami,"MotifEnrichment_CorExp.pdf"),height=hi, width = 2)
  #   print(Ht_i)
  #   dev.off()
  # }
   if ("Cic" %in% enriched_motifs$motif.name){
     mat_i <- out_mat[,"Cic"]
     mat_i <- na.omit(mat_i)
     hi <- length(mat_i)
     Ht_i <- Heatmap(mat_i,
                     rect_gp = gpar(col = "black", lwd = 0.5),
                     col = col_fun,
                     show_column_names = FALSE,
                     show_row_dend = FALSE,
                     column_title = paste0("Cic -",org))
     pdf(paste0(outdir2,"withCic/",dato,"_",org,"_Cic_MotifEnrichment_CorExp.pdf"),height=hi, width = 2)
     print(Ht_i)
     dev.off()
   }
  
  }
  else {
    print("No enriched motifs, so nothing to plot.")
    print(org)
  }}

#### Cic missing
# Check if Cic exists in your pfm object
library(TFBSTools)
names_list <- name(pfm)
grep("Cic", names_list, ignore.case = TRUE, value = TRUE)
# Also try related names
grep("capicua", names_list, ignore.case = TRUE, value = TRUE)

# not there so adding it manually with:
# Cic core binding motif - from SetUp/Databases/CisBP_2026_04_30_7_52_pm/, # Downloaded Cic motif from https://cisbp.ccbr.utoronto.ca/
Cic_mat <- t(matrix(as.numeric(unlist(str_split("0.279066852762909	0.181426232213218	0.181426232213218	0.358080682810656
0.289551319577206	0.189084572544428	0.212314722658093	0.309049385220273
0.000105099534514162	0.000105099534514162	0.000105099534514162	0.999684701396457
0.000105099534514162	0.000105099534514162	0.999684701396457	0.000105099534514162
0.000105099534514162	0.999684701396457	0.000105099534514162	0.000105099534514162
0.000105099534514162	0.000105099534514162	0.000105099534514162	0.999684701396457
0.000105099534514162	0.000105099534514162	0.999684701396457	0.000105099534514162
0.999684701396457	0.000105099534514162	0.000105099534514162	0.000105099534514162
0.0686788673212966	0.613130189946389	0.162332011531521	0.155858931200794
0.241038162691977	0.147189533347557	0.147189533347557	0.464582770612909",pattern="\\s+"))), nrow=10, byrow =T, dimnames = list(NULL,c("A","C","G","T"))))

cic_pfm <- PFMatrix(
  ID = "CIC_custom",
  name = "Cic",
  matrixClass = "HMG",
  strand = "+",
  bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
  tags = list(species = "Mus musculus"),
  profileMatrix = Cic_mat*100)
cic_pfm <- PFMatrixList(cic_pfm)
names(cic_pfm) <- "CIC_custom"
pfm_custom <- c(pfm, cic_pfm)
tail(names(pfm_custom))