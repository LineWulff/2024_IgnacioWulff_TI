# dato and libraries read in
#### ---- Initiate libraries ---- ####
library(ggplot2)
library(plotly)
library(stringr)
library(ggrastr)
library(scales)
library(Signac)
library(Seurat)
library(scales)
library(matrixStats)
library(openxlsx)
library(ggrastr)
library(dplyr)
library(future)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(GenomicRanges)
source("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/CollapseAnno.R")

# variables used 
dir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

# Initial check of location and start time for script
print(paste0("Seurat script started at: ",Sys.time()))
print(paste0("I am here: ",dir,""))

control_df <- read.csv(paste(dir,"/control_df.csv",sep=""), header = T)
control_df <- control_df[,c("parameter","value")]
print(control_df)

# Read in sheets of comparisons and control sheet
comp_sheet <- read.csv("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/comp_sheet.csv", header = T)

#create output folder
if (file.exists("output")){print("Directory already exist.")} else {
  dir.create(file.path(dir, "output"))}

#variables from control csv file
obj_name <- as.character(control_df[control_df$parameter=="Seurat_obj",]$value)
assay <- as.character(control_df[control_df$parameter=="assay",]$value)
var <- as.character(control_df[control_df$parameter=="var",]$value)
cond <- as.character(control_df[control_df$parameter=="cond",]$value)
proj <- as.character(control_df[control_df$parameter=="proj",]$value)

#### read in Seurat data ####
obj <- readRDS(paste(dir, "/", obj_name, ".rds", sep=""))

# Test that all varaibles exust in the object you are working with (spelling etc.)
head(comp_sheet)
for (i in rownames(comp_sheet)){
  if (comp_sheet[i,"var1"] %in% unique(obj@meta.data[,cond]) & comp_sheet[i,"var2"] %in% unique(obj@meta.data[,cond])){}
  else{stop(paste0("There is an error in comparison: ",i,"."))}}
print("All condition comparisons exist in object.")

# Set assay and identity and latent variables to include
DefaultAssay(obj) <- assay
Idents(obj) <- cond
latent.vars.incl <- if (assay=="ATAC"){latent.vars.incl='nCount_peaks'} else if(assay=="RNA"){latent.vars.incl='nCount_RNA'} else{latent.vars.incl=NULL}

#If assay is ATAC, peak annotation should be initiated here
if (assay=="ATAC"){
  edb <- EnsDb.Mmusculus.v79
  seqlevelsStyle(edb) <- "UCSC"
  peakAnno.edb <- annotatePeak(obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),TxDb = edb)
}

# Run loop of FindMarkers
for (clus in unique(obj@meta.data[,var])){
  DA_peaks_clus <- list()
  for (i in rownames(comp_sheet)){
    DA_peaks_conditions <- list()
    cond1 <- comp_sheet[i,"var1"]
    cond2 <- comp_sheet[i,"var2"]
    # check that the subsetted data includes enough cells to run 3 =< per groupxcond
    if (dim(obj@meta.data[obj@meta.data[,var]==clus & obj@meta.data[,cond]==cond1,])[1]>3 & 
        dim(obj@meta.data[obj@meta.data[,var]==clus & obj@meta.data[,cond]==cond2,])[1]>3){
      
      # subset data
      obj_sub <- subset(obj, cells = rownames(obj@meta.data[obj@meta.data[,var]==clus & obj@meta.data[,cond] %in% c(cond1,cond2),]))
      
      # launch parallelization
      future({
      # Run FindMArkers on the subset
      DA_peaks_conditions <- FindMarkers(object = obj_sub,
                                         only.pos = FALSE,
                                         ident.1 = cond1, 
                                         ident.2 = cond2,
                                         test.use = 'LR',
                                         latent.vars = latent.vars.incl,
                                         logfc.threshold = 0)
      # siginificance thresholds
      DA_peaks_conditions$sign <- "not sign."
      if (any(DA_peaks_conditions$avg_log2FC>0.25 & DA_peaks_conditions$p_val_adj<0.05)){
        DA_peaks_conditions[DA_peaks_conditions$avg_log2FC>0.25 & DA_peaks_conditions$p_val_adj<0.05,]$sign <- cond1}
      if (any(DA_peaks_conditions$avg_log2FC<(-0.25) & DA_peaks_conditions$p_val_adj<0.05)){
        DA_peaks_conditions[DA_peaks_conditions$avg_log2FC<(-0.25) & DA_peaks_conditions$p_val_adj<0.05,]$sign <- cond2}
      
      # volcano plot
      volc_plot <- ggplot(DA_peaks_conditions, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ))+
        geom_point_rast()+
        geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
        geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
        geom_vline(xintercept = c(0))+ #0
        scale_color_manual(values = c(cond1="#F8766D",cond2="#00BFC4", 'not sign.'='lightgrey'))+
        theme_classic()+
        ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
        guides(colour=guide_legend(title="Significance"))
      pdf(paste0(outdir,dato,"_",proj,"_VolcanoPlot_",clus,"_",cond1,"X",cond2,".pdf"),height = 4, width = 5.5)
      print(volc_plot)
      dev.off()
      
      
      ## For ATAC + add gene names and overall annotations FIX - needs to collapse the intron / exon states - write function
      if (assay=="ATAC"){
        # First add genename
        open_regs <- rownames(DA_peaks_conditions)
        closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                        regions = open_regs)
        DA_peaks_conditions$gene_name <-NA
        DA_peaks_conditions[rownames(DA_peaks_conditions) %in% closest_genes$query_region,]$gene_name <- closest_genes$gene_name
        # nor overall annotaion
        closest_genes <- ClosestFeature(subset(visu_obj, cells = rownames(visu_obj@meta.data[visu_obj@meta.data$t_split==clus,])),
                                        regions = open_regs, annotation = peakAnno.edb@anno)
        DA_peaks_conditions$annotation <-NA
        DA_peaks_conditions[rownames(DA_peaks_conditions) %in% closest_genes$query_region,]$annotation <- CollapseAnno(closest_genes$annotation)
        
      }
      
      # Subset data to only include significant data points for the saved data file
      DA_peaks_conditions <- DA_peaks_conditions[DA_peaks_conditions$avg_log2FC>0.25 & DA_peaks_conditions$p_val_adj<0.05 |
                                                   DA_peaks_conditions$avg_log2FC<(-0.25) & DA_peaks_conditions$p_val_adj<0.05,]
      DA_peaks_clus[[paste(cond1,cond2,sep = "X")]] <- DA_peaks_conditions
    }) #end future/parallization
      }
    else {
      # If analysis could not be run due to inadequate number of cells
      DA_peaks_clus[[paste(cond1,cond2,sep = "X")]] <- NULL
      print(paste("There were too few cells in one condition to do the comparison for",clus,"between",cond1,"and",cond2))
    }
  }
  
  # write xlsx for clus
  write.xlsx(DA_peaks_clus, file = paste0(outdir,dato,"_",proj,"_Overview_",clus,".xlsx"),
             rowNames = T)
  
}
# Script finished running at
print(paste0("Seurat script ended at: ",Sys.time()))



