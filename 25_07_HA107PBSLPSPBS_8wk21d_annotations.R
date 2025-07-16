#' R script for annotating scATAC data - 8 BM samples
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

#### --- Read in data ---- ####
combined <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/25_07_14 PBSHA107PBALPS_8wk21d.rds")
Idents(combined) <- 'ATAC_snn_res.0.2'

#### ---- Cluster IDs using singular genes ---- #### 
rownames(combined@assays$RNA@data)[startsWith(rownames(combined@assays$RNA@data),"Csf")]

# HSPC markers - see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969381/
FeaturePlot(combined, features = "rna_Cd34")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cd48")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Gata2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Procr")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Kit")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Fgd5")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Hoxb5")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Epc1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly6a")+scale_colour_gradientn(colors = mycols) #Sca-1
# more - Tgfbr3, Slc22a3, Gata2, and Kit - from https://www.cell.com/immunity/fulltext/S1074-7613(24)00260-7
FeaturePlot(combined, features = "rna_Slc22a3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Tgfbr3")+scale_colour_gradientn(colors = mycols)

# lymphoid myeloid primed - Cd34, Kit, Flt3, Meis1, Dach1, Ikzf2, Mecom, and Hlf - from https://www.cell.com/immunity/fulltext/S1074-7613(24)00260-7
FeaturePlot(combined, features = "rna_Cd34")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Kit")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Flt3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Meis1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Dach1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ikzf2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Hlf")+scale_colour_gradientn(colors = mycols)


# monocyte markers - see see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969381/
FeaturePlot(combined, features = "rna_Ly6c2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly6c1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Lyz2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ccr2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Csf2ra")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Elane")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cebpe")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Csf1r")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly86")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cx3cr1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Fcer1g")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Itgam")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_H2-Ab1")+scale_colour_gradientn(colors = mycols)
# non-classical monocytes - Sirpa, Csf1r, Adgre1, Tnfrsf1b, from https://www.cell.com/immunity/fulltext/S1074-7613(24)00260-7
FeaturePlot(combined, features = "rna_Sirpa")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Csf1r")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Adgre1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Tnfrsf1b")+scale_colour_gradientn(colors = mycols)

# DC progenitor markers - MDP/CDP (Calr, Mpo, and Ctsg) - from https://www.cell.com/immunity/fulltext/S1074-7613(24)00260-7
FeaturePlot(combined, features = "rna_Flt3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Irf8")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Batf3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Irf4")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Calr")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Mpo")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ctsg")+scale_colour_gradientn(colors = mycols)
# cDC1 (same paper) - (Bcl2, Tcf3, and Nfil3, Id2, Batf3, Zbtb46, and Cd226)
FeaturePlot(combined, features = "rna_Zbtb46")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cd226")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Batf3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Id2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Nfil3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Tcf3")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Bcl2")+scale_colour_gradientn(colors = mycols)
#cDC2 (same paper) - Itgam, Cx3cr1, Sirpa, Csf1r, Ciita, and H2-Ab1 
FeaturePlot(combined, features = "rna_Itgam")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cx3cr1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Sirpa")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Csf1r")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ciita")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_H2-Ab1")+scale_colour_gradientn(colors = mycols)


## Contaminants?
# Granulocytes: Anxa1, Zmpste24, Ly75, Ncam1, and Syne1
FeaturePlot(combined, features = "rna_Anxa1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Zmpste24")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly75")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ncam1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Syne1")+scale_colour_gradientn(colors = mycols)

# NK cells - https://www.nature.com/articles/s41467-019-11947-7 markers
FeaturePlot(combined, features = "rna_Klrd1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Nkg7")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Gnly")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ncam1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Cd7")+scale_colour_gradientn(colors = mycols)

# Neutrophils
FeaturePlot(combined, features = "rna_S100a8")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly6g")+scale_colour_gradientn(colors = mycols)

#
FeaturePlot(combined, features = "rna_Siglecf")+scale_colour_gradientn(colors = mycols)

#### --- Doublets and debris ---- ####
VlnPlot(combined,
        features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
        pt.size = 0,
        ncol = 3)

FeaturePlot(combined, features = "nCount_peaks")+scale_colour_gradientn(colors = mycols)

DimPlot(combined, label = T, group.by = "ATAC_snn_res.0.2")

combined@meta.data$pos_doub <- NA
combined@meta.data[combined@meta.data$nCount_peaks>90000,]$pos_doub <- "doublet"
dim(combined@meta.data)[1]
dim(combined@meta.data[!is.na(combined@meta.data$pos_doub),])[1]
dim(combined@meta.data[!is.na(combined@meta.data$pos_doub),])[1]/dim(combined@meta.data)[1]*100
DimPlot(combined, label = T, group.by = "pos_doub")

combined@meta.data$pos_deb <- NA
combined@meta.data[combined@meta.data$nCount_peaks<4500,]$pos_deb <- "debris"
dim(combined@meta.data)[1]
dim(combined@meta.data[!is.na(combined@meta.data$pos_deb),])[1]
dim(combined@meta.data[!is.na(combined@meta.data$pos_deb),])[1]/dim(combined@meta.data)[1]*100
DimPlot(combined, label = T, group.by = "pos_deb")


#### --- removing cells 1 ---- ####
## clus 17 at re.0.9 - debris
rem_cells <- c(rownames(combined@meta.data[!is.na(combined@meta.data$pos_doub),]),
               rownames(combined@meta.data[combined@meta.data$ATAC_snn_res.0.9==17,]))
length(Cells(combined))
length(rem_cells)
length(rem_cells)/length(Cells(combined))*100
combined <- subset(combined, cells = Cells(combined)[!Cells(combined) %in% rem_cells])


FeaturePlot(combined, features = "nCount_peaks")+scale_colour_gradientn(colors = mycols)


#### ---- Use external data set to identify cells instead ---- ####
# data set from https://www.nature.com/articles/s41556-019-0439-6 (preprocessed Seurat object)
# droplet-based scRNAseq12 of cells from total mouse BM
load(paste(RAID_dir,"10x_ext_datasets/RNAMagnetDataBundle/NicheData10x.rda",sep = "/"))
NicheData10x <- UpdateSeuratObject(NicheData10x)
# top DEGs of idnets from same data set
ND_DEGs <- read.csv(paste(RAID_dir,"10x_ext_datasets/RNAMagnetDataBundle/NicheData10x_topDEGs.csv", sep="/"))
head(ND_DEGs)

# inspect
NicheData10x
DimPlot(NicheData10x, label = T)#+NoLegend()
head(NicheData10x@meta.data)
NicheData10x@meta.data$ID <- Idents(NicheData10x)
# Add also upper level IDs as in Fig 1
mesenchymal <- c("Myofibroblasts","Smooth muscle","Fibro/Chondro p.","Stromal fibro.","Arteriolar fibro.",
                 "Osteo-CAR","Chondrocytes","Endosteal fibro.","Osteoblasts","Ng2+ MSCs","Adipo-CAR")
immune <- c("NK cells","B cell","Dendritic cells","small pre-B.","Neutrophils","T cells","Monocytes","pro-B",
            "large pre-B.")
HSPC <- c("Mk prog.","Erythroblasts","Eo/Baso prog.","Ery prog.","LMPPs","Gran/Mono prog.","Mono prog.",
          "Neutro prog.","Ery/Mk prog.")
EC <- c("Arteriolar ECs","Sinusoidal ECs")
Neuronal <- c("Schwann cells")

unique(NicheData10x@meta.data$ID)[!unique(NicheData10x@meta.data$ID) %in% c(mesenchymal,immune,HSPC,Neuronal,EC)]
# all covered, now add the upper level ID
upIDs <- c(rep("mesenchymal",length(mesenchymal)),
           rep("immune",length(immune)),
           rep("HSPC",length(HSPC)),
           rep("neuronal",length(Neuronal)),
           rep("EC",length(EC)))
names(upIDs) <- c(mesenchymal,immune,HSPC,Neuronal,EC)


NicheData10x@meta.data$uplev_ID <- NA
for (clus in unique(NicheData10x@meta.data$ID)){
  NicheData10x@meta.data[NicheData10x@meta.data$ID==clus,]$uplev_ID <- upIDs[clus]
}
# check worked
DimPlot(NicheData10x, group.by = "uplev_ID", label = T)+NoLegend()

# Run function from clustercorrelation script
# first get variable feature gene_names from combined object
var_gens <- ClosestFeature(combined, VariableFeatures(combined))
var_gens <- unique(var_gens$gene_name)

uplev_corr <- cluster_corr(combined,"ATAC_snn_res.0.2","RNA",var1 = var_gens ,NicheData10x,"uplev_ID","RNA")
heatmap(as.matrix(uplev_corr), scale = "none", trace="none", density ="none",
        dendrogram = "column",
        col = c(rep("#313695",15),mycols),
        RowSideColors = unlist(as.list(hue_pal()(nrow(uplev_corr)))), ylab = deparse(substitute(combined)),# labRow = "",
        ColSideColors = unlist(as.list(hue_pal()(ncol(uplev_corr)))), xlab = deparse(substitute(NicheData10x)),# labCol = "",
        margins = c(11,5))

ID_corr <- cluster_corr(combined,"ATAC_snn_res.0.2","RNA",var1 = var_gens ,NicheData10x,"ID","RNA")
heatmap(as.matrix(ID_corr), scale = "none", trace="none", density ="none",
        dendrogram = "column",
        col = c(rep("#313695",15),mycols),
        RowSideColors = unlist(as.list(hue_pal()(nrow(ID_corr)))), ylab = deparse(substitute(combined)),# labRow = "",
        ColSideColors = unlist(as.list(hue_pal()(ncol(ID_corr)))), xlab = deparse(substitute(NicheData10x)),# labCol = "",
        margins = c(8,5))

ID_sub_corr <- cluster_corr(combined,"ATAC_snn_res.0.2","RNA",var1 = var_gens ,
                            subset(NicheData10x, cells = rownames(NicheData10x@meta.data[NicheData10x@meta.data$uplev_ID %in% c("immune","HSPC","EC"),])),
                            "ID","RNA")
heatmap(as.matrix(ID_sub_corr), scale = "none", trace="none", density ="none",
        dendrogram = "column",
        col = c(rep("#313695",15),mycols),
        RowSideColors = unlist(as.list(hue_pal()(nrow(ID_sub_corr)))), ylab = deparse(substitute(combined)),# labRow = "",
        ColSideColors = unlist(as.list(hue_pal()(ncol(ID_sub_corr)))), xlab = deparse(substitute(NicheData10x)),# labCol = "",
        margins = c(8,5))

# max values from above
ID_max <- apply(ID_sub_corr, 1, max)
ID_max_match <- c()
for (i in seq(1,length(ID_max))){
  if (ID_max[i] %in% ID_sub_corr[names(ID_max[i]),]){
    print(names(ID_max[i]))
    print(colnames(ID_sub_corr)[ID_sub_corr[names(ID_max[i]),] %in% ID_max[i]])
    ID_max_match <- c(ID_max_match, colnames(ID_sub_corr)[ID_sub_corr[names(ID_max[i]),] %in% ID_max[i]])}
  
}

## Do top DEGs as module scores
head(ND_DEGs)
DefaultAssay(combined) <- 'RNA'

# Add modulescores for each DEG set 
for (col in colnames(ND_DEGs)){
  col_genes <- ND_DEGs[,col]
  col_genes <- col_genes[col_genes != ""]
  if (length(col_genes)>100) {len_col_genes = 100} else{len_col_genes = length(col_genes)}
  print(col)
  print(paste("existing genes: ", length(col_genes[col_genes %in% rownames(combined@assays$RNA@data)]), " of ",length(col_genes), sep = ""))
  combined <- AddModuleScore(combined, features = list(col_genes),
                             assay = "RNA", ctrl = len_col_genes,
                             name = col)
  
}
# print and save the featureplots and vlnplots
for (col in colnames(ND_DEGs)){
  feat_plot <- FeaturePlot(combined, features = paste(col,1, sep = ""))+
    scale_colour_gradientn(colors = mycols)
  vln_plot <- VlnPlot(combined, features = paste(col,1, sep = ""), pt.size = 0)
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "AllCells_res0.2_UMAPfeatPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(feat_plot)
  dev.off()
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "AllCells_res0.2_VlnPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(vln_plot)
  dev.off()
}

#### ---- Removing cells 2 ---- ####
## CLuster 7 non immune cells
## Cluster 9 - non cells
rem_cells <- rownames(combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% c(7,9),])
length(Cells(combined))
length(rem_cells)
length(rem_cells)/length(Cells(combined))*100
combined <- subset(combined, cells = Cells(combined)[!Cells(combined) %in% rem_cells])

# save res.0.2 in meta for future reference
combined@meta.data$v1_res.0.2 <- combined@meta.data$ATAC_snn_res.0.2 

## recalculate clustering and umap
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
res <- seq(0,1,0.1)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = res)

DimPlot(combined, group.by = "v1_res.0.2", label = T)
DimPlot(combined, group.by = "ATAC_snn_res.0.1", label = T)
DimPlot(combined, group.by = "ATAC_snn_res.0.2", label = T)

## Now we can annotate monocytes in more detail:
FeaturePlot(combined, features = "rna_Ly6c2")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ly6c1")+scale_colour_gradientn(colors = mycols)
FeaturePlot(combined, features = "rna_Ccr2")+scale_colour_gradientn(colors = mycols)

## at res.0.2
# 0 - Ly6c hi monocytes
# 1,5,2 - Ly6c lo monocytes

#### ---- Adding cluster IDs ---- ####
Ly6chimono <- 0
Ly6clomono <- c(1,2,5)
DCs <- 6
NKcells <- 7
LSK <- 4
Neutro <- 3

combined@meta.data$ID_labs <- NA
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% Ly6clomono,]$ID_labs <- "Ly6c lo monocytes"
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% Ly6chimono,]$ID_labs <- "Ly6c hi monocytes"
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% DCs,]$ID_labs <- "Dendritic cells"
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% NKcells,]$ID_labs <- "NK cells"
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% LSK,]$ID_labs <- "LSK"
combined@meta.data[combined@meta.data$ATAC_snn_res.0.2 %in% Neutro,]$ID_labs <- "Neutrophils"
unique(combined@meta.data$ID_labs)
combined@meta.data$ID_labs <- factor(combined@meta.data$ID_labs, levels = c("Ly6c lo monocytes","Ly6c hi monocytes",
                                                                            "LSK","Neutrophils","Dendritic cells",
                                                                            "NK cells"))
#combined@meta.data[is.na(combined@meta.data$ID_labs),]
DimPlot(combined, group.by = "ID_labs")

#### ---- save object for further analysis from this point ---- ####
saveRDS(combined, paste(dato,"PBSHA107PBALPS_8wk21d_clean_v2.rds",sep="_"))
