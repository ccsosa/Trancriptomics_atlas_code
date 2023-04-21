# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rtracklayer")
#call rtracklayer
require("rtracklayer")
#define directory to call files
#DO for CNE filtered motif
in_dir <- "/scratch/bis_klpoe/chsos/data/MOTIFS/COMBINED"
#list files
x_C <- list.files(in_dir,pattern = ".bed")
#load files and combined in one
x_C_len <- lapply(1:length(x_C),function(i){
  #i <-1
  message(paste0("i: ",i, "/ ", round((i/length(x_C)*100),2)," %"),"")
  xx <- rtracklayer::import(paste0(in_dir,"/",x_C[[i]]),format = "BED")
  x <-data.frame(MOTIF=sub(".bed","",x_C[[i]]),TARGET=unique(xx$name))
  return(x)
})
#joining in one
x_C_len <- do.call(rbind,x_C_len)

#saving feature files
write.table(x_C_len,"/scratch/bis_klpoe/chsos/data/MOTIFS/feature_comb.tsv",na = "",row.names = F,sep = "\t",quote = F)
#write.table(x_C_len,"/scratch/bis_klpoe/chsos/data/MOTIFS/feature_comb_NH.tsv",na = "",row.names = F,sep = "\t",quote = F,col.names = F)
