##### Function to create dataframe for visualization (like barplots) of distribution within a Seurat single cell object
# meta data column (discrete) contribution to specific cells (could be all cell in an object)
# can be plotted as in 22_03_F5_monomac_update.R - not currently on github

# one discrete column with no splits, e.g. a clustering
perc_function <- function(meta_col, cell_vec, Seu_obj){
  Seu_obj_sub <- subset(Seu_obj, cells = cell_vec)
  res_df <- data.frame(cluster=unique(Seu_obj_sub@meta.data[,meta_col]))
  res_df <- cbind(res_df, percent=NA)
for (clus in res_df$cluster){
  print(clus)
  res_df[res_df$cluster==clus,"percent"] <- nrow(Seu_obj_sub@meta.data[Seu_obj_sub@meta.data[,meta_col]==clus,])/length(cell_vec)*100
  }
  return(res_df)
  }

# one discrete column + a split, for instance on per sample base to see clustering
# equally distribution not dominated by single sample of many cells
perc_function_samp <- function(meta_col, cell_vec, Seu_obj,splitgroup){
  Seu_obj_sub <- subset(Seu_obj, cells = cell_vec)
  res_df <- data.frame(cluster=rep(unique(Seu_obj_sub@meta.data[,meta_col]),length(unique(Seu_obj_sub@meta.data[,splitgroup]))))
  res_df$samp <- NA
  for (i in seq(1,length(unique(Seu_obj_sub@meta.data[,splitgroup])))){
    gro <- unique(Seu_obj_sub@meta.data[,splitgroup])[i]
    print(gro)
    start = i*length(unique(Seu_obj_sub@meta.data[,meta_col]))-length(unique(Seu_obj_sub@meta.data[,meta_col]))+1
    stop = i*length(unique(Seu_obj_sub@meta.data[,meta_col]))
    res_df[start:stop,]$samp <- gro
  }
  #res_df <- cbind(res_df, group = c(rep(unique(Seu_obj_sub@meta.data[,splitgroup]),length(unique(Seu_obj_sub@meta.data[,meta_col])))))
  res_df$percent <- NA
  print(res_df)
  for (clus in res_df$cluster){
    print(clus)
    for (gro in unique(res_df$samp)){
    res_df[res_df$cluster==clus & res_df$samp==gro,"percent"] <- nrow(Seu_obj_sub@meta.data[Seu_obj_sub@meta.data[,meta_col]==clus & Seu_obj_sub@meta.data[,splitgroup]==gro,])/nrow(Seu_obj_sub@meta.data[Seu_obj_sub@meta.data[,splitgroup]==gro,])*100
  }}
  return(res_df)
}
