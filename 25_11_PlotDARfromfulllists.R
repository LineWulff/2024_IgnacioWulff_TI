# Usage: Rscript volcano_from_excel.R path/to/results.xlsx

library(readxl)
library(ggplot2)
library(ggrastr)
library(dplyr)
library(stringr)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(GenomicRanges)
library(scales)
library(Signac)
library(Seurat)


file_path <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/Outputs/tspace_DAR/25_11_11_MonoTrajQ1-Q3comp_TotalOverview_T1.xlsx"

# Create output folder
outdir <- "/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/Outputs/tspace_DAR"
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj <- "MonoTrajQ1-Q3comp"
clus <- "T1"
#object related to analysis
obj <- readRDS(paste0("/Users/linewulff/Documents/work/projects/2024_IgnacioWulff_TI/BM-PBSHA107_PBSLPS_21d8wk/Outputs/tspace/","2510_BM-HA107PBS-LPSPBS-21d8wk_MonocyteTraj_visuobj.rds"))


edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(obj@assays$ATAC@ranges, tssRegion=c(-3000, 3000),
                             TxDb = edb)

# Get sheet names
sheets <- readxl::excel_sheets(file_path)
cat("Found", length(sheets), "sheets:\n", paste(sheets, collapse = ", "), "\n\n")


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


ann_col_val <- hue_pal()(length(unique(annotations_gen)))
names(ann_col_val) <- c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR",
                        "1st Intron","Other Intron","1st Exon","Other Exon","3' UTR","Downstream (<=300bp)")
ann_col_val
show_col(ann_col_val)


for (sh in sheets) {
  cat("Processing sheet:", sh, "\n")
  df <- readxl::read_excel(file_path, sheet = sh)
  cond1 <- unlist(str_split(sh,"X"))[1]
  cond2 <- unlist(str_split(sh,"X"))[2]
  cat("Condition 1:", cond1, "\n")
  cat("Condition 2:", cond2, "\n")
  
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
  nonsig <- df %>% filter(sign == "not sign.")
  sig <- df %>% filter(sign != "not sign.")
  
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
  acc_stat_df$totDEG <- 0; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[1],]$totDEG <- DAsum[1]; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[2],]$totDEG <- DAsum[2];
  acc_stat_df$freq <- 0; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[1],]$freq <- DAsum[1]; acc_stat_df[acc_stat_df$sign==rownames(DAsum)[2],]$freq <- DAsum[2];
  acc_stat_df$freq <- acc_stat_df$amount/acc_stat_df$freq
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




