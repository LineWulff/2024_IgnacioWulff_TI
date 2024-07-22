# Libraries
library(tidyverse)
library(tidyr)
library(scales)
library(Seurat)

dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### ---- Function ---- ####
## calcualte correlation between two Seurat objets ##
## if one is mouse put that one first
cluster_corr <- function(obj1,res1,assay1,var1=F,obj2,res2,assay2,mouse=F){
  Idents(obj1) <- obj1@meta.data[,res1]
  Idents(obj2) <- obj2@meta.data[,res2]
  obj1_av <- AverageExpression(obj1, return.seurat = TRUE, assays = assay1); obj2_av <- AverageExpression(obj2, return.seurat = TRUE, assays = assay2)
  if (var1[1] == F){
    obj1_var <- VariableFeatures(obj1)} else {obj1_var <- var1[var1 %in% rownames(obj1_av@assays[["RNA"]]$data)]}
  obj2_var <- VariableFeatures(obj2)
  print(length(obj1_var))
  print(length(obj2_var))
  if (mouse==T){
    inc_feats <- obj1_var[toupper(obj1_var) %in% obj2_var]}
  else {inc_feats <- obj1_var[obj1_var %in% obj2_var]}
  print(paste(length(inc_feats),"overlapping variable genes bewteen the two objects"))
  obj1_mat <- as.matrix(obj1_av@assays[[assay1]]$data[inc_feats,])
  rownames(obj1_mat) <- toupper(rownames(obj1_mat))
  if (mouse==T){inc_feats <- toupper(inc_feats)}
  write.csv(inc_feats,file=paste(dato,"corr_",obj1@project.name,res1,obj2@project.name,res2,".csv",sep = "_"), row.names = F)
  obj2_mat <- as.matrix(obj2_av@assays[[assay2]]$data[inc_feats,])
  rownames(obj2_mat) <- toupper(rownames(obj2_mat))
  corr_mat <- as.data.frame(cor(x=obj1_mat, y = obj2_mat, method = "pearson"))
  print(corr_mat)
  # hm <- heatmap.2(as.matrix(corr_mat), scale = "none", col =  bluered(150), trace="none", density ="none",
  #           dendrogram = "column",
  #           RowSideColors = unlist(as.list(hue_pal()(nrow(corr_mat)))), ylab = deparse(substitute(obj2)),# labRow = "",
  #           ColSideColors = unlist(as.list(hue_pal()(ncol(corr_mat)))), xlab = deparse(substitute(obj1)),# labCol = "",
  #           margins = c(5,5))
  # png(paste(dato,"correlationplot",deparse(substitute(obj1)),deparse(substitute(res1)),deparse(substitute(obj2)),deparse(substitute(res2)),".png",sep="_"), height = 1000, width = 1000, res = 150)
  # print(hm)
  # dev.off()
  return(corr_mat)
}


#### --- Examples ---- ###

new_hm <- cluster_corr(SI_graft, "RNA_snn_res.0.2", "RNA",SI_strom,"integrated_snn_res.0.2","integrated")

png(paste(dato,project,))
heatmap.2(as.matrix(new_hm), scale = "none", col =  bluered(150), trace="none", density ="none",
          dendrogram = "column",
          RowSideColors = unlist(as.list(hue_pal()(nrow(new_hm)))), ylab = deparse(substitute(SI_graft)),# labRow = "",
          ColSideColors = unlist(as.list(hue_pal()(ncol(new_hm)))), xlab = deparse(substitute(SI_strom)),# labCol = "",
          margins = c(5,5))
dev.off()




png(paste(dato,"adultLIres0.2graftLIres0.2clusters_pearsoncorrelation_heatmap2_wlabels.png",sep="_"), height = 1000, width = 750, res = 150)
heatmap.2(cor_mat, scale = "none", col =  bluered(150), trace="none", density ="none",
          dendrogram = "column",
          RowSideColors = unlist(as.list(hue_pal()(nrow(cor_mat)))), ylab = "LI graft",# labRow = "",
          ColSideColors = LI_cols, xlab = "LI adult",# labCol = "",
          margins = c(12,4) )
dev.off()
png(paste(dato,"adultLIres0.2graftLIres0.2clusters_pearsoncorrelation_heatmap2.png",sep="_"), height = 800, width = 850, res = 150)
heatmap.2(as.matrix(cor_mat[,1:6]), scale = "none", col =  bluered(150), trace="none", density ="none",
          dendrogram = "column",
          RowSideColors = unlist(as.list(hue_pal()(nrow(cor_mat)))), ylab = "LI graft", labRow = "",
          ColSideColors = LI_cols, xlab = "LI adult", labCol = "",
          margins = c(3,3) )
dev.off()