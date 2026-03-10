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
dist_df <- perc_function_samp("t_split", Cells(visu_obj), visu_obj,"orig.ident")
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

tsplit_8wk <- ggplot(visu, aes(x=T_mean,y=tPC2))+
  geom_point_rast(data=visu[visu$timepoint=="21d",], color = "lightgrey")+
  geom_point_rast(data=visu[visu$timepoint=="8wk",],aes(color=t_split_8wk))+
  #geom_point_rast(data=visu[ly6lotop_cells,],aes(x=T_mean,y=tPC2),color="lightgrey",shape=8,size=4,)+
  #  geom_point(data=visu[ly6lotop,],aes(x=umap1,y=umap2),color="black")+
  theme_classic()

pdf(paste(outdir,dato,"_tSPvsavpsud_tsplit8wk.pdf",sep = ""),height = 4, width = 5.5)
print(tsplit_8wk)
dev.off()

IDlabs_plot <- ggplot(visu, aes(x=T_mean,y=tPC2, color=ID_labs))+
  geom_point_rast()+
  theme_classic()

pdf(paste(outdir,dato,"_tSPvsavpsud_IDlabscol.pdf",sep = ""),height = 4, width = 5.5)
print(IDlabs_plot)
dev.off()




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

#### Read in saved data objects ####
visu <- readRDS(paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visudf.rds"))
visu_obj <- readRDS(paste0(outdir,"2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visuobj.rds"))
