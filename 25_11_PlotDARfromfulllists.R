# Usage: Rscript volcano_from_excel.R path/to/results.xlsx

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

rm(list=ls())

# Create output folder
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSvsHA107-PBSvsLPS-8wk/DiffPeaks_v2/"
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj <- "BM_8wk_IDlabs"; project <- proj
clus <- "Neutrophils"
file_path <- paste0(outdir,"26_01_23_BM-LSKNeu-8wk_TotalOverview_Neutrophils.xlsx")
#object related to analysis
obj <- readRDS(paste0("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/2510_BM-HA107PBS-LPSPBS-8wk_MonocyteTraj_visuobj.rds"))
obj <- readRDS("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/25_10_06_PBSHA107PBALPS_8wk_clean.rds") 

edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)

### color control for annotations ####
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


ann_col_val <- hue_pal()(length(unique(annotations_gen)))
names(ann_col_val) <- c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                        "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)")
ann_col_val
show_col(ann_col_val)

tot_ann <- cbind(peakAnno.edb@annoStat,ID=rep(1,length(peakAnno.edb@annoStat$Feature)))
levels(tot_ann$Feature) <- c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                                                       "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)")

pdf(paste(outdir,dato,project,"_DistTotalAnnotationObject.pdf",sep=""),height = 2, width = 6)
ggplot(tot_ann, aes(x=ID,y=Frequency, fill=Feature))+
  geom_bar(stat="identity", colour = "black")+
  scale_fill_manual(values = ann_col_val)+
  labs(y="",x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  coord_flip()
dev.off()

### Sheet processing start ####
# Get sheet names
sheets <- readxl::excel_sheets(file_path)
cat("Found", length(sheets), "sheets:\n", paste(sheets, collapse = ", "), "\n\n")


for (sh in sheets) {
  cat("Processing sheet:", sh, "\n")
  df <- readxl::read_excel(file_path, sheet = sh)
  if (grepl("_", sh, ignore.case = FALSE, fixed = FALSE)){
    tiss <- unlist(str_split(sh,"-"))[1]
    timep <- unlist(str_split(unlist(str_split(sh,"-"))[2],"_"))[1]
    cond1 <- unlist(str_split(unlist(str_split(sh,"_"))[2],"X"))[1]
    cond2 <- unlist(str_split(unlist(str_split(sh,"_"))[2],"X"))[2]
    cond1 <- paste(tiss,cond1,timep, sep = "-")
    cond2 <- paste(tiss,cond2,timep, sep = "-")
    cat("Condition 1:", cond1, "\n")
    cat("Condition 2:", cond2, "\n")
  }
  else {cond1 <- unlist(str_split(sh,"X"))[1]
  cond2 <- unlist(str_split(sh,"X"))[2]
  cat("Condition 1:", cond1, "\n")
  cat("Condition 2:", cond2, "\n")}
  # Check that needed columns exist
  if (!all(c("avg_log2FC", "p_val_adj","sign") %in% names(df))) {
    warning(paste("Skipping sheet", sh, "- required columns missing"))
    next}
  # Check that needed columns exist
  if (length(unique(df$sign))<2) {
    warning(paste("Skipping sheet", sh, "- No significant values"))
    next}
  
  # Clean data
  df <- df %>%
    mutate(
      avg_log2FC = as.numeric(avg_log2FC),
      p_val_adj = as.numeric(p_val_adj),
      sign = ifelse(is.na(sign), "not sign.", sign),
      annotation = annotation
    )

  # Split by significance so non-significant points are plotted first
  nonsig <- df %>% dplyr::filter(sign == "not sign.")
  sig <- df %>% dplyr::filter(sign != "not sign.")
  
  # Define consistent colors
  cols <- c("lightgrey", "#F8766D", "#00BFC4")
  names(cols) <- c("not sign.",cond1,cond2)


  # Make plot
  volc_plot <- ggplot()+
    geom_point_rast(data = nonsig, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign ),
               color = "lightgrey" )+
    geom_point_rast(data = sig, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sign )) +
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ # sign. threshold
    geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed")+ # sign. threshold
    geom_vline(xintercept = c(0))+ #0
    scale_color_manual(values = cols)+
    theme_classic()+
    ylab("-log10(adj. p-value)")+xlab("avg. log2FC")+
    guides(colour=guide_legend(title="Significance"))

  pdf(paste0(outdir,"/",dato,"_",proj,"_VolcanoPlot_",clus,"_",cond1,"X",cond2,".pdf"),height = 4, width = 5.5)
  print(volc_plot)
  dev.off()
  
  ##### Annotations #####
  ## Bar and volcano plots highlighting annotations
  acc_stat_df <- as.data.frame(summary(as.factor(sig$annotation)))
  colnames(acc_stat_df) <- 'amount'
  acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))

  ## same as above but split into up vs down reg
  acc_dat <- sig
  acc_dat$type_reg <- paste(acc_dat$annotation,acc_dat$sign,sep = '_')
  acc_stat_df <- as.data.frame(summary(as.factor(acc_dat$type_reg)))
  colnames(acc_stat_df) <- 'amount'
  acc_stat_df <- cbind(acc_stat_df, annotation=rownames(acc_stat_df))
  acc_stat_df$sign <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(2,length(acc_stat_df$amount)*2,2)]
  acc_stat_df$annotation <- unlist(str_split(acc_stat_df$annotation, "_"))[seq(1,length(acc_stat_df$amount)*2,2)]

  ## freq of DAR annotations
  DAsum <- rowsum(acc_stat_df$amount, group=acc_stat_df$sign)
  DAsum_df <- data.frame(sign = rownames(DAsum),totDEG = as.numeric(DAsum))
  acc_stat_df <- acc_stat_df %>%
    left_join(DAsum_df, by = "sign") %>%
    mutate(totDEG = tidyr::replace_na(totDEG, 0))
  acc_stat_df$freq <- acc_stat_df$amount/acc_stat_df$totDEG
  acc_stat_df$annotation <- factor(acc_stat_df$annotation ,
                                   levels = c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                                              "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)"))
  anno_plot <- ggplot(acc_stat_df, aes(x=sign, y=freq, fill=annotation))+
    geom_bar(stat="identity", colour = "black")+
    scale_y_continuous(labels = scales::percent)+
    geom_text(data=acc_stat_df[acc_stat_df$annotation=="Promoter (<=1kb)",], aes(x=sign, label=totDEG, y=1.1, fill=NULL))+
    scale_fill_manual(values = ann_col_val)+
    labs(y="",x="")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))+
    coord_flip()
  pdf(paste(paste0(outdir,"/",dato,proj,"_annotatedDARs_dist_",clus,"_",cond1,"X",cond2,".pdf")),height = 2, width = 6)
  print(anno_plot)
  dev.off()

}


#### Cell numbers per group ####
visu_obj <- readRDS(paste0("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/2510_BM-HA107PBS-LPSPBS-8wk_MonocyteTraj_visuobj.rds"))

## first add Q1-Q3 info to monos
obj@meta.data$ID_labs_ext <- as.character(obj$ID_labs)
obj@meta.data[Cells(visu_obj),]$ID_labs_ext <- visu_obj$t_split
df_count <- as.data.frame(table(obj$ID_labs_ext,obj$orig.ident))
colnames(df_count) <- c("CellID_ext","Sample","CellCount")

df_count %>% ggplot(aes(x=Sample,y=CellCount,fill=CellID_ext))+
  geom_bar(stat="identity",colour="black")+
  facet_wrap(.~CellID_ext)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept = c(250,500), linetype = "dashed")
  

#### CovergaePlots ####
CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q1",])),
  group.by = "orig.ident",
  region = "Fos",
  #region.highlight = regions_highlight,
  extend.upstream = 4000,
  extend.downstream = 3000
)

CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q1",])),
  group.by = "orig.ident",
  region = "Ier2",
  #region.highlight = regions_highlight,
  extend.upstream = 3000,
  extend.downstream = 3000
)

CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q2",])),
  group.by = "orig.ident",
  region = "Fosb",
  #region.highlight = regions_highlight,
  extend.upstream = 3000,
  extend.downstream = 3000
)

CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q1",])),
  group.by = "orig.ident",
  region = "Jund",
  #region.highlight = regions_highlight,
  extend.upstream = 3000,
  extend.downstream = 3000
)

CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q1",])),
  group.by = "orig.ident",
  region = c("Cspg5","Dsel"),
  #region.highlight = regions_highlight,
  extend.upstream = 5000,
  extend.downstream = 5000
)

CoveragePlot(
  object = subset(obj, cells = rownames(obj@meta.data[obj@meta.data$ID_labs_ext=="Q1",])),
  group.by = "orig.ident",
  region = c("Cr2"),
  #region.highlight = regions_highlight,
  extend.upstream = 5000,
  extend.downstream = 5000
)
