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
library(scales)
library(matrixStats)
library(openxlsx)


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
rownames(combined@assays$RNA@data)[startsWith(rownames(combined@assays$RNA@data),"Csfr")]

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


VlnPlot(combined,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0,
  ncol = 3
)
FeaturePlot(combined, features = "nCount_peaks")+scale_colour_gradientn(colors = mycols)
# is cl 5 just debris??? - Yes

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

uplev_corr <- cluster_corr(combined,"ATAC_snn_res.0.4","RNA",var1 = var_gens ,NicheData10x,"uplev_ID","RNA")
heatmap(as.matrix(uplev_corr), scale = "none", trace="none", density ="none",
          dendrogram = "column",
          col = c(rep("#313695",15),mycols),
          RowSideColors = unlist(as.list(hue_pal()(nrow(uplev_corr)))), ylab = deparse(substitute(combined)),# labRow = "",
          ColSideColors = unlist(as.list(hue_pal()(ncol(uplev_corr)))), xlab = deparse(substitute(NicheData10x)),# labCol = "",
          margins = c(11,5))

ID_corr <- cluster_corr(combined,"ATAC_snn_res.0.4","RNA",var1 = var_gens ,NicheData10x,"ID","RNA")
heatmap(as.matrix(ID_corr), scale = "none", trace="none", density ="none",
        dendrogram = "column",
        col = c(rep("#313695",15),mycols),
        RowSideColors = unlist(as.list(hue_pal()(nrow(ID_corr)))), ylab = deparse(substitute(combined)),# labRow = "",
        ColSideColors = unlist(as.list(hue_pal()(ncol(ID_corr)))), xlab = deparse(substitute(NicheData10x)),# labCol = "",
        margins = c(8,5))

ID_sub_corr <- cluster_corr(combined,"ATAC_snn_res.0.4","RNA",var1 = var_gens ,
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
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "AllCells_res0.4_UMAPfeatPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(feat_plot)
  dev.off()
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "AllCells_res0.4_VlnPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(vln_plot)
  dev.off()
}

#### --- Distribution plots ---- ####
## Run function script first
res0.4_dist <- perc_function_samp("ATAC_snn_res.0.4", colnames(combined), combined,"colonization")
head(res0.4_dist)
tapply(combined@meta.data$ATAC_snn_res.0.4, combined@meta.data$orig.ident, summary)

res0.4_dist$samp <- factor(res0.4_dist$samp, levels = c("PBS","HA107"))

ggplot(res0.4_dist, aes(x=samp, y=percent, fill=cluster))+geom_bar(position = "stack", stat = "identity", colour = "black")
ggplot(res0.4_dist, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(position = "dodge", stat = "identity", colour = "black")+
  facet_grid(.~cluster)+
  theme_classic()+
  xlab("colonization")+
  theme(axis.text.x = element_text(angle = 90))





print(hej)
  


#### ---- save object version and metadata---- ####
# Clustering res.0.4
# 0 - monocytes
# 1 - monocytes
# 2 - monocytes
# 3 – neutrophil (progenitors?) **
# 4 - monocytes
# 5 - Debris *
# 6 – LSK sorted cells
# 7 - Dendritic cell / progenitors
# 8 - NK cells *
# 9 - Eosinophils *

## Add IDs and save metadata
res0.4_IDs <- c("0"="0_monocytes", "1"="1_monocytes", "2"="2_monocytes", "3"="neutrophils", "4"="4_monocytes",
                "5"="debris", "6"="LSK", "7"="DC prog.", "8"="NK cells", "9"="eosinophils")
combined@meta.data$res0.4_IDs <- NA
for (clus in as.character(unique(combined@meta.data$ATAC_snn_res.0.4))){
  combined@meta.data[as.character(combined@meta.data$ATAC_snn_res.0.4) == clus,]$res0.4_IDs <- res0.4_IDs[[clus]]
}
DimPlot(combined, group.by = "res0.4_IDs", label = T)

# save metadata
saveRDS(combined@meta.data, paste(dato,"_SeuObj_",project,"_MetaDataOnly.rds",sep = ""))
saveRDS(combined, paste(dato,"_SeuObj_",project,".rds",sep = ""))

#### ---- Reclustering based on res.0.4 ---- ####
## clusterID
# Clustering res.0.4
# 0 - monocytes
# 1 - monocytes
# 2 - monocytes
# 3 – neutrophil (progenitors?) **
# 4 - monocytes
# 5 - Debris *
# 6 – LSK sorted cells
# 7 - Dendritic cell / progenitors
# 8 - NK cells *
# 9 - Eosinophils *

# save metadata

save_cells <- rownames(combined@meta.data[combined@meta.data$ATAC_snn_res.0.4 %in% c(0,1,2,3,4,6,7),])
length(save_cells) #10325 cells in object after removing contaminants
table(combined@meta.data[save_cells,]$orig.ident)
length(Cells(combined)) - length(save_cells) # removing 915 cells/debris
table(combined@meta.data$orig.ident)-table(combined@meta.data[save_cells,]$orig.ident) # 491 from HA107, 424 from PBS

### Run merge from line 96 again 
combined <- merge(
  x = BMPBSPBS8wk,
  y = BMHAPBS8wk,
  add.cell.ids = c("BM-PBS-PBS-8wk", "BM-HA107-PBS-8wk")
)

combined <- subset(combined, cells = save_cells)
length(Cells(combined)) # check fots with line 465

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

## Add IDs from metadata of former round
r1_meta <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/24_07_22_SeuObj_BM-PBSvsHA107-PBS-8wk_MetaDataOnly.rds")
combined@meta.data$res0.4_IDs <- r1_meta[save_cells,]$res0.4_IDs

DimPlot(combined, group.by = 'res0.4_IDs', pt.size = 0.1)

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
        group.by = "ATAC_snn_res.0.1",
        split.by = "orig.ident")
Idents(combined) <- 'ATAC_snn_res.0.1'

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
VlnPlot(combined, features = "Ccr2", assay = "RNA", split.by = "orig.ident")


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

## monocyte markers
FeaturePlot(combined, features = "rna_Fcgr3")+scale_colour_gradientn(colors = mycols) # CD16


#### ---- Distribution plot ----- ####
## Run function script first
res0.1_dist <- perc_function_samp("ATAC_snn_res.0.1", colnames(combined), combined,"colonization")
head(res0.1_dist)
tapply(combined@meta.data$ATAC_snn_res.0.1, combined@meta.data$orig.ident, summary)

res0.1_dist$samp <- factor(res0.1_dist$samp, levels = c("PBS","HA107"))

ggplot(res0.1_dist, aes(x=samp, y=percent, fill=cluster))+geom_bar(position = "stack", stat = "identity", colour = "black")
ggplot(res0.1_dist, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(position = "dodge", stat = "identity", colour = "black")+
  facet_grid(.~cluster)+
  theme_classic()+
  xlab("colonization")+
  theme(axis.text.x = element_text(angle = 90))

#### ---- compare to literature ---- ####
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
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "AllCells_res0.4_UMAPfeatPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(feat_plot)
  dev.off()
  pdf(paste(RAID_dir,"/", project, "/output/", dato, "LSKMonoNeu_res0.1_VlnPlot_BaccinData_",col , ".pdf", sep = ""), height = 5, width = 5)
  print(vln_plot)
  dev.off()
}

#### ---- DE between clusters and conditions ---- ####
# change back to working with peaks instead of gene activities
DefaultAssay(combined) <- 'ATAC'


## If running FindConservedMarkers specify assay, as will otherwise assume RNA assay
## nCount_peaks vs peak_region_fragments + essentially the same, nCount calc. by Seurat, 
## highly correlated though
## cor(combined$peak_region_fragments, combined$nCount_peaks) = 0.999848
DA_peaks_res.0.1_FC1 <- FindAllMarkers(
  object = combined,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 1
)

head(DA_peaks_res.0.1_FC1)

## add gene name to the regions to infer activity
open_regs <- rownames(DA_peaks_res.0.1_FC1)
closest_genes_FC1 <- ClosestFeature(combined, regions = open_regs)
head(closest_genes_FC1[,c("gene_name","query_region")])
DA_peaks_res.0.1_FC1$gene_name <- NA
DA_peaks_res.0.1_FC1[rownames(DA_peaks_res.0.1_FC1) %in% closest_genes_FC1$query_region,]$gene_name <- closest_genes_FC1$gene_name

head(DA_peaks_res.0.1_FC1)
write.csv(DA_peaks_res.0.1_FC1, file = paste("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Outputs/Clustering/",dato, project,"LSKMonoNeu_res.0.1_PeaksPrClus_FC1.csv",sep = ""))

## Comparison between conditions per cluster
DA_peaks_conditions <- list()
Idents(combined) <- "orig.ident"
  
for (clus in unique(combined@meta.data$ATAC_snn_res.0.1)){
  DA_peaks_res.0.1 <- FindMarkers(
    object = subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ATAC_snn_res.0.1==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-8wk', 
    ident.2 = 'BM-HA107-PBS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_res.0.1)
  closest_genes <- ClosestFeature(subset(combined, cells = rownames(combined@meta.data[combined@meta.data$ATAC_snn_res.0.1==clus,])),
                                  regions = open_regs)
  DA_peaks_res.0.1$gene_name <-NA
  DA_peaks_res.0.1[rownames(DA_peaks_res.0.1) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  # Add sign. groups (mainly for plotting), all non sign. will be removed before saving
  DA_peaks_res.0.1$sign <- "not sign."
  DA_peaks_res.0.1[,]$sign <- "BM-PBS-PBS-8wk"
  DA_peaks_res.0.1[,]$sign <- "BM-HA107-PBS-8wk"
  
  # plot as volcano
  volc_plot <- ggplot(DA_peaks_res.0.1, aes(x = , y =, colour = ))+
    geom_point_rast()+
    geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = c('BM-PBS-PBS-8wk'="#F8766D",'BM-HA107-PBS-8wk'="#00BFC4", 'not sign.'='lightgrey'))+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    theme()
  pdf(paste0(RAID_dir,dato,"cluster",clus,".pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  
  DA_peaks_conditions[[clus]] <- DA_peaks_res.0.1 
}
names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in csv format
write.xlsx(DA_peaks_conditions, file = paste(projdir,"/Outputs/Clustering/", dato, "_LSKMonoNeu_res.0.1_PerClusCompvsCond.xlsx", sep = ""),
           rowNames = T)


