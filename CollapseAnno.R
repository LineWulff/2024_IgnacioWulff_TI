## Collapse annotations for exon and introns
CollapseAnno <- function(anno_vec){
  anno_vec_upd <- anno_vec
  for (ann in anno_vec[startsWith(anno_vec, "Intron")]){
    if ( unlist(str_split(ann," "))[4] == "1" ){
      anno_vec_upd[anno_vec==ann] <- "1st Intron"
  }
    else {
      anno_vec_upd[anno_vec==ann] <- "Other Intron"
  }}
# exons
  for (ann in anno_vec[startsWith(anno_vec, "Exon")]){
  if ( unlist(str_split(ann," "))[4] == "1" ){
    anno_vec_upd[anno_vec==ann] <- "1st Exon"
  }
  else {
    anno_vec_upd[anno_vec==ann] <- "Other Exon"
  }}
  return(anno_vec_upd)
} #end function
