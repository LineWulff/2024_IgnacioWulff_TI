#' R script for isolating monocytes and prep data for tSPACE
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on hhttps://stuartlab.org/signac/articles/pbmc_vignette.html

#### ---- Initiate libraries ---- ####
library(ggplot2)
library(plotly)
library(stringr)
library(ggrastr)
library(viridis)
library(scales)
library(Signac)
library(Seurat)
library(colorRamp2)
library(scales)
library(matrixStats)
library(openxlsx)
library(ggrastr)
library(plot3D)
library(rgl)
library(crosstalk)
library(DT)
library(dplyr)
library(future)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(GenomicRanges)
source("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/Distribution_functions.R")


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

### Allowing to use multiple cores when possible
options(future.globals.maxSize = 10 * 1024^3)
plan(multisession, workers = 2)

#### ---- Read in data and subset it ---- ####
combined <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/25_07_16_PBSHA107PBALPS_8wk21d_clean_v2.rds")
Idents(combined) <- 'ID_labs'
DefaultAssay(combined) <- 'ATAC'

sheets <- c("cluster_Ly6c lo monocytes","cluster_Neutrophils","cluster_Dendritic cells","cluster_LSK","cluster_Ly6c hi monocytes","cluster_NK cells")
incl_subs <- c("Ly6c lo monocytes","Ly6c hi monocytes")

# subset cells
combined <- subset(combined, idents = incl_subs)

#### ---- Prep for tSpace - test  ---- ####
table(combined@meta.data$orig.ident)/42051*100
table(combined@meta.data$ID_labs)/42051*100
## For test purposes make downsampled version of dataset
down_comb <- sample(Cells(combined), size = 1000)
down_comb <- subset(combined, cells = down_comb)
table(down_comb@meta.data$orig.ident)/10
table(down_comb@meta.data$ID_labs)/10

saveRDS(down_comb, file = "250818_DS_monocytes.rds")

#### ---- Prep for tSpace - whole object  ---- ####

combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="timepoint")
DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="colonization")
DimPlot(combined, group.by = 'ID_labs', pt.size = 0.1,split.by="stimulation")

DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1)


comb_low <- RunTFIDF(comb_low)
comb_low <- FindTopFeatures(comb_low, min.cutoff = "q99", assay = "ATAC")
length(VariableFeatures(comb_low)) #208510

combined[["ATAC"]]
combined[["RNA"]]
comb_low <- combined
comb_low@active.assay <- "ATAC"
comb_low[["RNA"]] <- NULL #for smaller version
saveRDS(comb_low, file = "250822_NoRNA_monocytes.rds")

combined[["ATAC"]]
combined[["RNA"]]
comb_low <- combined
comb_low@active.assay <- "RNA"
comb_low[["ATAC"]] <- NULL #for smaller version
comb_low <-FindVariableFeatures(comb_low, nfeatures = 1000, selection.method = "vst")
saveRDS(comb_low, file = "250925_NoATAC_monocytes.rds")


#### ---- Open and adjust tspace objects ---- ####
subs1 <- readRDS(file="/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/tspace_output/250925_NoATAC_monocytes_tspacefile.rds")
subs1 <- readRDS(file="/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/tspace_output/251002_NoATAC_monocytes_tspacefile.rds")


subs1$ts_file[1:5,1:20]
visu <- subs1$ts_file
visu <- cbind(visu, combined@meta.data[rownames(visu),])
rownames(subs1$tspace_matrix) <- rownames(visu)
# exclude cells obased on below
ggplot(visu, aes(x=umap1, y=umap2, color=ID_labs))+
  geom_point()+geom_hline(yintercept = -10)+geom_vline(xintercept = -7)


Excl_cells <- rownames(visu[visu$umap1<(-7) & visu$umap2<(-10),])
visu <- visu[!rownames(visu) %in% Excl_cells,]
#visu <- visu[!visu$F2F5 %in% c("cDC2","cDC3","amb","HLA low"),]
#visu <- cbind(visu, DC_subs@reductions$umap@cell.embeddings[rownames(visu),])
#visu <- cbind(visu, idents=Idents(DC_subs)[rownames(visu)])
incl_genes <- colnames(visu)[14:1013]


#distance from monocytes
#mono.trajectories <- subs1$tspace_matrix[,which(colnames(subs1$tspace_matrix) %in% paste0('T_', visu[which(visu$integrated_snn_res.2.8=='24'), 'Index']))]
#colnames(mono.trajectories) <- c("T_1","T_2","T_3","T_4","T_5")
#visu <- cbind(visu, dist_mono=mono.trajectories)
#visu$T_mean <-rowMeans(mono.trajectories)

#### Adding new trajectory based clustering ####
visu_obj <- subset(combined, cells = rownames(visu))
visu_obj@reductions$lsi@cell.embeddings <- as.matrix(visu[,startsWith(colnames(visu),"tPC")])

visu_obj <- RunUMAP(object = visu_obj, reduction = 'lsi', dims = 2:10, n.components = 2)
visu_obj@reductions$umap@cell.embeddings <- as.matrix(visu[,startsWith(colnames(visu),"umap")])
colnames(visu_obj@reductions$umap@cell.embeddings) <- c("umap_1","umap_2")
visu_obj <- FindNeighbors(object = visu_obj, reduction = 'lsi', dims = 2:10)
res <- seq(0.1,1,0.1)
visu_obj <- FindClusters(object = visu_obj, verbose = FALSE, resolution = res)
res <- seq(0.01,0.09,0.01)
visu_obj <- FindClusters(object = visu_obj, verbose = FALSE, resolution = res)

## Add clustering to visu df
tclust <- visu_obj@meta.data[,c("ATAC_snn_res.0.09","ATAC_snn_res.0.08","ATAC_snn_res.0.07","ATAC_snn_res.0.1","ATAC_snn_res.0.2")]
colnames(tclust) <- paste0("tsp_",colnames(tclust))
visu <- cbind(visu,tclust)

DimPlot(visu_obj, dims = c(1,2),group.by = "ID_labs")
DimPlot(visu_obj, dims = c(1,2), group.by = "ATAC_snn_res.0.07")

DimPlot(visu_obj, dims = c(1,2),group.by = "colonization")
DimPlot(visu_obj, dims = c(1,2),group.by = "stimulation", split.by = "colonization")
DimPlot(visu_obj, dims = c(1,2),group.by = "timepoint", split.by = "colonization")
FeaturePlot(visu_obj, features = c("Ly6c2"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
#FeaturePlot(visu_obj, features = c("Ly6c1"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Ccr2"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
# trajectory strating point, based:
# https://rupress.org/jem/article/213/11/2293/42007/CXCR4-identifies-transitional-bone-marrow
FeaturePlot(visu_obj, features = c("Myb"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Cxcr4"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Kit"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Mpo"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Ctsg"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Elane"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, features = c("Cdca7"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))


# potential starting points based on above
ggplot(visu, aes(x=umap1, y=umap2, color=ID_labs))+
  geom_point()+geom_hline(yintercept = 9)+geom_vline(xintercept = c(-7,-2))

ggplot(visu, aes(x=umap1, y=umap2, color=ATAC_snn_res.0.1))+
  geom_point()+geom_hline(yintercept = 9)+geom_vline(xintercept = c(-7,-2))


ggplot(visu, aes(x=tPC1, y=tPC2, color=ID_labs))+
  geom_point()
ggplot(visu, aes(x=tPC3, y=tPC4, color=stimulation))+
  geom_point()

#### 3D visualization ####
visu <- cbind(visu, visu_obj@reductions$umap@cell.embeddings)
visu$colIDlabs <- NA
for (i in 1:length(unique(visu$ID_labs))){
  cl <- sort(unique(visu$ID_labs))[i]
  visu[visu$ID_labs == cl,]$colIDlabs <- hue_pal()(length(unique(visu$ID_labs)))[i]
  print(hue_pal()(length(unique(visu$ID_labs)))[i])
  print(cl)
}


plot3d(x=visu[,"umap_1"],y=visu[,"umap_2"],z=visu[,"umap_3"],
       col=visu$colIDlabs,
       xlab = "UMAP_1",ylab = "UMAP_2",zlab = "UMAP_3")
writeWebGL(dir ="/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/trajectories/tspace/monomac_woprol/",filename=paste(dato,project,traj_sub,"res2.8_3D.html", sep="_"))

#### --- Distance from monocytes ---- ####
#distance from Ly6c hi monocytes at top corner - geom_hline(yintercept = 9)+geom_vline(xintercept = c(-7,-2))
ly6lotop <- rownames(visu[visu$umap2>9 & visu$umap1>(-7) & visu$umap1<(-2) ,]) #& visu$umap1<(-1)
ly6lotop_cells <- rownames(visu)[which(visu$Index %in% as.numeric(str_sub(colnames(subs1$tspace_matrix)[which(colnames(subs1$tspace_matrix) %in% paste0('T_', visu[ly6lotop, 'Index']))],start=3,end=-1)))]
ly6lo.trajectories <- subs1$tspace_matrix[rownames(visu),which(colnames(subs1$tspace_matrix) %in% paste0('T_', visu[ly6lotop, 'Index']))]
colnames(ly6lo.trajectories) <- c("T_1","T_2","T_3")
visu <- cbind(visu, dist_ly6lotop=ly6lo.trajectories)
visu$T_mean <-rowMeans(ly6lo.trajectories)

ggplot(visu, aes(x=umap1,y=umap2,color=T_mean))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=umap1,y=umap2),color="lightgrey",shape=8,size=4,)+
#  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
    scale_color_viridis_c(option = "magma")+
  theme_classic()


ggplot(visu, aes(x=umap1,y=umap2,color=ATAC_snn_res.0.4))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=umap1,y=umap2),color="lightgrey",shape=8,size=4,)+
  theme_classic()

ggplot(visu, aes(x=umap1,y=umap2,color=ID_labs))+
  geom_density_2d()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=umap1,y=umap2),color="black",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  facet_grid(stimulation~colonization)

#### Trying something different ####
ggplot(visu, aes(x=T_mean,y=tPC2,color=T_mean))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="lightgrey",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  scale_color_viridis_c(option = "magma")+
  theme_classic()

ggplot(visu, aes(x=T_mean,y=tPC2,color=ID_labs))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="lightgrey",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  theme_classic()

ggplot(visu, aes(x=T_mean,y=tPC2,color=ID_labs))+
  geom_density_2d()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="black",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  facet_grid(stimulation~colonization)

ggplot(visu, aes(x=T_mean,y=tPC2,color=ID_labs))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="black",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  facet_grid(stimulation~colonization)

ggplot(visu, aes(x=T_mean,y=tPC2,color=timepoint))+
  geom_density_2d()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="black",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  facet_grid(stimulation~colonization)

visu_obj@reductions$traj_emb <- visu_obj@reductions$umap
visu_obj@reductions$traj_emb@cell.embeddings <- as.matrix(visu[,c("T_mean","tPC2")])
visu_obj@reductions$traj_emb@cell.embeddings[,1] <- visu_obj@reductions$traj_emb@cell.embeddings[,1]*10
visu_obj@reductions$traj_emb@cell.embeddings[,2] <- visu_obj@reductions$traj_emb@cell.embeddings[,2]*100
colnames(visu_obj@reductions$traj_emb@cell.embeddings) <- c("umap_1","umap_2")


# trajectory strating point, based:
# https://rupress.org/jem/article/213/11/2293/42007/CXCR4-identifies-transitional-bone-marrow
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Myb"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Cxcr4"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Kit"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Mpo"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Ctsg"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Elane"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Cdca7"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))

FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Ly6c2"), pt.size = 0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))
FeaturePlot(visu_obj, reduction = "traj_emb", features = c("Ccr2"), pt.size =0.5)+scale_color_gradientn(colors=c(mycols,rep('#a50026',10)))

ggplot(visu, aes(x=T_mean,y=tPC2,color=tsp_ATAC_snn_res.0.09))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="black",shape=8,size=4,)
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  #facet_grid(stimulation~colonization)


#### ---- Distribution of new tsp based clustering ---- ####
dist_df <- perc_function_samp("ATAC_snn_res.0.09", Cells(visu_obj), visu_obj,"orig.ident")
dist_df$samp <- factor(dist_df$samp, 
                       levels = c("BM-PBS-PBS-21d","BM-PBS-PBS-8wk","BM-HA107-PBS-21d","BM-HA107-PBS-8wk","BM-HA107-LPS-21d","BM-HA107-LPS-8wk","BM-PBS-LPS-21d","BM-PBS-LPS-8wk"))

pdf(paste(outdir,dato,"_Distribution_PersampleStacked.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black")+
  theme_classic()+
  labs(x="", y="% of sample")+
  theme(axis.text.x = element_text(angle=90))
dev.off()

pdf(paste(outdir,dato,"_Distribution_PersampleDodge.pdf",sep = ""),height = 4, width = 6)
ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  facet_grid(.~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()

len <- dim(dist_df)[1]
dist_df$colonization <- factor(unlist(str_split(dist_df$samp,"-"))[seq(2,len*4,4)], levels=c("PBS","HA107"))
dist_df$timepoint <- unlist(str_split(dist_df$samp,"-"))[seq(4,len*4,4)]
dist_df$stimulation <- factor(unlist(str_split(dist_df$samp,"-"))[seq(3,len*4,4)], levels=c("PBS","LPS"))

pdf(paste(outdir,dato,"_Distribution_PersampleStacked_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=colonization, y=percent, fill=cluster))+
  geom_bar(stat="identity", colour="black", width=0.8)+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~stimulation)+
  theme(axis.text.x = element_text(angle=90))
dev.off()
pdf(paste(outdir,dato,"_Distribution_PersampleDodge_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=colonization, y=percent, fill=stimulation))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()
pdf(paste(outdir,dato,"_Distribution_PersampleDodge_StimxColxTP.pdf",sep = ""),height = 4, width = 4.5)
ggplot(dist_df, aes(x=stimulation, y=percent, fill=colonization))+
  geom_bar(stat="identity", colour="black",position = "dodge")+
  theme_classic()+
  labs(x="", y="% of sample")+
  facet_grid(timepoint~cluster)+
  theme(axis.text.x = element_text(angle=90))
dev.off()

#### groupings supervised ####
visu$t_split <- NA
visu <- visu %>%
  mutate(
    # y-values on each abline at given x
    line1 = -0.01 * T_mean + 0.0065,
    line2 =  0.005 * T_mean - 0.0025,
    
    # Define region based on position relative to the two lines
    t_split = case_when(
      tPC2 > line1 & tPC2 > line2 ~ "T3",  # above both lines
      tPC2 > line1 & tPC2 <= line2 ~ "T3", # between line1 and line2 (upper part)
      tPC2 <= line1 & tPC2 > line2 ~ "T1", # between line2 and line1 (lower part)
      tPC2 <= line1 & tPC2 <= line2 ~ "T2" # below both lines
    )
  ) %>%
  select(-line1, -line2) 

ggplot(visu, aes(x=T_mean,y=tPC2,color=t_split))+
  geom_point()+
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="lightgrey",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  theme_classic()

visu_8 <- readRDS("2510_BM-HA107PBS-LPSPBS-8wk_MonocyteTraj_visudf.rds")
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/Outputs/tspace/"
visu$t_split_8wk <- NA
inboth <- rownames(visu)[rownames(visu) %in% rownames(visu_8)]
visu[inboth,]$t_split_8wk <- visu_8[inboth,]$t_split

ggplot(visu, aes(x=T_mean,y=tPC2))+
  geom_point(data=visu[visu$timepoint=="21d",], color = "lightgrey")+
  geom_point(data=visu[visu$timepoint=="8wk",],aes(color=t_split_8wk))
  geom_point(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="lightgrey",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  theme_classic()

  visu_obj <- subset(visu_obj, cells=rownames(visu))
  visu_obj@meta.data$t_split <- visu$t_split
  visu_obj@meta.data$t_split_v2 <- visu$t_split_v2
  
  visu$stimulation <- factor(visu$stimulation, levels=c("PBS","LPS"))
  visu$colonization <- factor(visu$colonization, levels=c("PBS","HA107"))

### distribution of above  
visu_obj$t_split <- visu$t_split
dist_df <- perc_function_samp("t_split", Cells(visu_obj), visu_obj,"orig.ident")
dist_df$samp <- factor(dist_df$samp, 
                         levels =  c("BM-PBS-PBS-21d","BM-PBS-PBS-8wk","BM-HA107-PBS-21d","BM-HA107-PBS-8wk","BM-HA107-LPS-21d","BM-HA107-LPS-8wk","BM-PBS-LPS-21d","BM-PBS-LPS-8wk"))
  
  pdf(paste(outdir,dato,"_Distribution_PersampleStacked.pdf",sep = ""),height = 4, width = 4.5)
  ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
    geom_bar(stat="identity", colour="black")+
    theme_classic()+
    labs(x="", y="% of sample")+
    theme(axis.text.x = element_text(angle=90))
  dev.off()
  
  pdf(paste(outdir,dato,"_Distribution_PersampleDodge.pdf",sep = ""),height = 4, width = 6)
  ggplot(dist_df, aes(x=samp, y=percent, fill=cluster))+
    geom_bar(stat="identity", colour="black",position = "dodge")+
    theme_classic()+
    facet_grid(.~cluster)+
    theme(axis.text.x = element_text(angle=90))
  dev.off()
  
  len <- dim(dist_df)[1]
  dist_df$colonization <- factor(unlist(str_split(dist_df$samp,"-"))[seq(2,len*4,4)], levels=c("PBS","HA107"))
  dist_df$timepoint <- unlist(str_split(dist_df$samp,"-"))[seq(4,len*4,4)]
  dist_df$stimulation <- factor(unlist(str_split(dist_df$samp,"-"))[seq(3,len*4,4)], levels=c("PBS","LPS"))
  
  pdf(paste(outdir,dato,"_Distribution_PersampleStacked_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
  ggplot(dist_df, aes(x=colonization, y=percent, fill=cluster))+
    geom_bar(stat="identity", colour="black", width=0.8)+
    theme_classic()+
    labs(x="", y="% of sample")+
    facet_grid(timepoint~stimulation)+
    theme(axis.text.x = element_text(angle=90))
  dev.off()
  pdf(paste(outdir,dato,"_Distribution_PersampleDodge_ColxStimxTP.pdf",sep = ""),height = 4, width = 4.5)
  ggplot(dist_df, aes(x=colonization, y=percent, fill=stimulation))+
    geom_bar(stat="identity", colour="black",position = "dodge")+
    theme_classic()+
    labs(x="", y="% of sample")+
    facet_grid(timepoint~cluster)+
    theme(axis.text.x = element_text(angle=90))
  dev.off()
  pdf(paste(outdir,dato,"_Distribution_PersampleDodge_StimxColxTP.pdf",sep = ""),height = 4, width = 4.5)
  ggplot(dist_df, aes(x=stimulation, y=percent, fill=colonization))+
    geom_bar(stat="identity", colour="black",position = "dodge")+
    theme_classic()+
    labs(x="", y="% of sample")+
    facet_grid(timepoint~cluster)+
    theme(axis.text.x = element_text(angle=90))
  dev.off()
  
  
#### Plots of tsplit
  pdf(paste(outdir,dato,"_tSP_TmeanxtPC2_IDlabs_split.pdf",sep = ""),height = 6, width = 7)
  ggplot(visu, aes(x=T_mean,y=tPC2,color=ID_labs))+
    geom_point_rast()+
    geom_point(data=visu[ly6lotop_cells,],color="black",shape=8,size=4,)+
    geom_abline(intercept = 0.0065, slope = -0.01)+
    geom_abline(intercept = - 0.0025, slope = 0.005)+
    #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
    facet_grid(stimulation~colonization)
  dev.off()
  
  pdf(paste(outdir,dato,"_tSP_TmeanxtPC2_IDlabsContour_split.pdf",sep = ""),height = 6, width = 13)
  ggplot(visu, aes(x=T_mean,y=tPC2,color=ID_labs))+
    geom_density_2d()+
    geom_abline(intercept = 0.0065, slope = -0.01)+
    geom_abline(intercept = - 0.0025, slope = 0.005)+
    geom_point(data=visu[ly6lotop_cells,],color="black",shape=8,size=4,)+
    #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
    facet_grid(timepoint~stimulation+colonization)
  dev.off()
  
  pdf(paste(outdir,dato,"_tSP_TmeanxtPC2_SampleContour_split.pdf",sep = ""),height = 6, width = 7)
  ggplot(visu, aes(x=T_mean,y=tPC2,color=orig.ident))+
    geom_density_2d()+
    geom_point(data=visu[ly6lotop_cells,],color="black",shape=8,size=4,)+
    #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
    facet_grid(stimulation~colonization)
  dev.off()
  
  pdf(paste(outdir,dato,"_tSP_TmeanxtPC2_tsplit.pdf",sep = ""),height = 6, width = 7)
  ggplot(visu, aes(x=T_mean,y=tPC2,color=t_split))+
    geom_point_rast()+
    geom_abline(intercept = 0.0065, slope = -0.01)+
    geom_abline(intercept = - 0.0025, slope = 0.005)+
    geom_point(data=visu[ly6lotop_cells,],color="black",shape=8,size=4,)
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  dev.off()
  
  pdf(paste(outdir,dato,"_tSP_TmeanxtPC2_OrigidentContour_split.pdf",sep = ""),height = 6, width = 13)
  ggplot(visu, aes(x=T_mean,y=tPC2))+
    stat_density_2d(aes(fill = ..level..), geom = "polygon",bins=20)+
    scale_fill_gradientn(colours = mycols)+
    geom_abline(intercept = 0.0065, slope = -0.01)+
    geom_abline(intercept = - 0.0025, slope = 0.005)+
    geom_point(data=visu[ly6lotop_cells,],color="black",shape=8,size=4,)+
    #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
    facet_grid(timepoint~stimulation+colonization)+
    theme_minimal()
  dev.off()
  
saveRDS(visu, paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visudf.rds"))
saveRDS(visu_obj, paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visuobj.rds"))

visu <- readRDS(paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visudf.rds"))
visu_obj <- readRDS(paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visuobj.rds"))
#### ---- Diff. Peak analysis on Q1-Q3 ---- ####
DefaultAssay(visu_obj) <- 'ATAC'

edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(visu_obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)
sheets <- c("cluster_Q1","cluster_Q2","cluster_Q3")

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
length(annotations_gen) # should be 172k
# make a color chart
ann_col_val <- hue_pal()(length(unique(annotations_gen)))
names(ann_col_val) <- unique(annotations_gen)
ann_col_val
show_col(ann_col_val)



##### 'BM-PBS-PBS-8wk' vs 'BM-HA107-PBS-8wk'####
DA_peaks_conditions <- list()
Idents(visu_obj) <- "orig.ident"

for (clus in unique(visu_obj@meta.data$t_split)){
  visu_sub <- subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,]))
  DA_peaks_IDlabs <- FindMarkers(
    object = visu_sub,
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-8wk', 
    ident.2 = 'BM-HA107-PBS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_MonoTrajSplit",clus,"PBSvsHA1078wk.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS8wk.xlsx", sep = ""),
           rowNames = T)

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx(paste0(outdir,"25_10_29_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS8wk.xlsx"),
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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

write.xlsx(DA_peaks_conditions, file = paste(outdir,"/", dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS8wk.xlsx", sep = ""),
           rowNames = T)

##### 'BM-PBS-PBS-21d' vs 'BM-HA107-PBS-21d'####
DA_peaks_conditions <- list()
Idents(visu_obj) <- "orig.ident"

for (clus in unique(visu_obj@meta.data$t_split)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-21d', 
    ident.2 = 'BM-HA107-PBS-21d',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_MonoTrajSplit",clus,"PBSvsHA10721d.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS21d.xlsx", sep = ""),
           rowNames = T)

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx(paste0(outdir,"25_10_29_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS21d.xlsx"),
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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

write.xlsx(DA_peaks_conditions, file = paste(outdir,"/", dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107PBSvsPBSPBS21d.xlsx", sep = ""),
           rowNames = T)


#### ---- PBS-LPS-8wk vs HA107-LPS-8wk ---- ####
DA_peaks_conditions <- list()
Idents(visu_obj) <- "orig.ident"

for (clus in unique(visu_obj@meta.data$t_split)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-LPS-8wk', 
    ident.2 = 'BM-HA107-LPS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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
  pdf(paste0(outdir,"/",dato,project,"VolcanoPlotDAR_MonoTrajSplit",clus,"PBSLPSvsHA107LPS.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir,"/", dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)
### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx(paste0(outdir,"/25_10_23_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx"),
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)

DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks/25_10_23_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_HA107LPSvsPBSLPS8wk.xlsx"
                       ,
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

#### ---- PBS-PBS-8wk vs PBS-LPS-8wk ---- ####
DA_peaks_conditions <- list()
Idents(visu_obj) <- "orig.ident"

for (clus in unique(visu_obj@meta.data$t_split)){
  DA_peaks_IDlabs <- FindMarkers(
    object = subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
    only.pos = FALSE,
    ident.1 = 'BM-PBS-PBS-8wk', 
    ident.2 = 'BM-PBS-LPS-8wk',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    logfc.threshold = 0)
  # add ass. genename
  open_regs <- rownames(DA_peaks_IDlabs)
  closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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
  pdf(paste0(outdir,dato,project,"VolcanoPlotDAR_MonoTrajSplit",clus,"PBSPBSvsPBSLPS.pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  # subset to only include significant values, logFC 0.5 and p.adj 0.05
  DA_peaks_IDlabs <- DA_peaks_IDlabs[DA_peaks_IDlabs$avg_log2FC>0.25 & DA_peaks_IDlabs$p_val_adj<0.05 |
                                       DA_peaks_IDlabs$avg_log2FC<(-0.25) & DA_peaks_IDlabs$p_val_adj<0.05,]
  DA_peaks_conditions[[clus]] <- DA_peaks_IDlabs 
}

names(DA_peaks_conditions) <- paste("cluster", names(DA_peaks_conditions), sep = "_")

## save in excel format
write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)

### Read in and add gene names and annotations
DA_peaks_conditions <- list()
for (i in seq(1,length(sheets))){
  print(i)
  DA_clus <- read.xlsx(paste0(outdir,"25_10_23_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx"),
                       colNames = T, sheet = sheets[i], rowNames = T)
  DA_peaks_conditions[[sheets[i]]] <- DA_clus
}

for (clus in names(DA_peaks_conditions)){
  ## Add gene names
  open_regs <- rownames(DA_peaks_conditions[[clus]])
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                  regions = open_regs)
  DA_peaks_conditions[[clus]]$gene_name <-NA
  DA_peaks_conditions[[clus]][rownames(DA_peaks_conditions[[clus]]) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
  
  ## Add annotation
  closest_genes <- ClosestFeature(subset(visu_obj@assays$ATAC, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
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

write.xlsx(DA_peaks_conditions, file = paste(outdir, dato, "_LSKMonoNeu_MonoTrajQ1-3_PerClusCompvsCond_PBSPBSvsPBSLPS8wk.xlsx", sep = ""),
           rowNames = T)

##### DA split in genomic annotations ####
# make a color chart
ann_col_val <- hue_pal()(length(unique(annotations_gen)))
names(ann_col_val) <- c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                        "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)")
ann_col_val
show_col(ann_col_val)

clus <- "cluster_Q1"
clus <- "cluster_Q2"
clus <- "cluster_Q3"

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
pdf(paste(paste0(outdir,"/",dato,"_HA107PBS_PBSPBS_8wk_annotatedDARs_dist_",clus,".pdf")),height = 2, width = 6)
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
