# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rtracklayer")
require("rtracklayer");library(parallel);require(ggpubr)


in_dir <- "/scratch/bis_klpoe/chsos/data/MOTIFS"
CB_dir <-paste0(in_dir,"/CLUSTER_BUSTER")
FIM_dir <- paste0(in_dir,"/FIMO")
out_dir <- "/scratch/bis_klpoe/chsos/data/MOTIFS/COMBINED"

x_COMB1 <- list.files(CB_dir,pattern = ".bed")
x_COMB2 <- list.files(FIM_dir,pattern = ".bed")

x_COMB <- intersect(x_COMB1,x_COMB2)
x_COMB_len <- lapply(1:length(x_COMB),function(i){
  
  message(paste0("i: ",i, "/ ", round((i/length(x_COMB)*100),2)," %"),"")
  
  if(!file.exists(paste0(out_dir,"/",sub(".bed.gz",".bed",x_COMB[[i]])))){
    message("processing")
    
  x1 <- rtracklayer::import(paste0(CB_dir,"/",sub(".bed.gz",".bed",x_COMB[[i]])),format = "BED")
  
  if(length(x1)>0){
  x1$method <- "CLUSTER_BUSTER"
  }
  x2 <- rtracklayer::import(paste0(FIM_dir,"/",sub(".bed.gz",".bed",x_COMB[[i]])),format = "BED")
  
  if(length(x2)>0){
    x2$method <- "FIMO"
  }
  
  if(length(x1)>0 & length(x2)>0){
    x <- c(x1,x2) 
  } else if(length(x1)==0 & length(x2)>0){
    x <- x2
  } else if(length(x1)>0 & length(x2)==0){
    x <- x1
  }
  
  x <- x[order(GenomeInfoDb::seqnames(x),IRanges::ranges(x),decreasing = F)] 
  
  rtracklayer::export(x,paste0(out_dir,"/",sub(".bed.gz",".bed",x_COMB[[i]])))
  
  } else {
    message("No coords found")
  }
  
  x <- data.frame(motif=sub(".bed","",x_COMB[[i]]),
                  genes= length(x),
                  unique_genes = length(unique(x$name)),
                  FIMO=sum(x$method=="FIMO"),
                  CB=sum(x$method=="CLUSTER_BUSTER"),
                  #threshold =thr,
                  method="COMBINED")
  return(x)
  })
  #x


x_COMB_len <- do.call(rbind,x_COMB_len)
x_df <- read.table("length2015.tsv")
################################################################################

################################################################################
x_df_comb <- rbind(x_df,x_COMB_len[,-c(4,5)])
####################################
# Combine histogram and density plots
gghistogram(x_df_comb, x = "unique_genes",
            add = "median", rug = TRUE,
            color = "method",
            fill = "method", palette = c("#00AFBB", "#E7B800","red"),
            add_density = TRUE,bins = 200,xlab = "Unique genes")



my_comparisons <- list( c("Cluster buster", "FIMO"),c("Cluster buster", "COMBINED"),c("FIMO", "COMBINED"))
ggviolin(x_df_comb, x = "method", y = "unique_genes", fill = "method",
         palette = c("#00AFBB","red","#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)   

ggviolin(x_df_comb, x = "method", y = "unique_genes", fill = "method",
         palette = c("#00AFBB","red","#E7B800"),
          add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)   
save.image("~/CLUST_COMB.RData")
################################################################################
load("~/CLUST_COMB.RData")

require("rtracklayer");library(parallel);require(ggpubr);require(data.table)
x2 <- fread("/scratch/bis_klpoe/chsos/data/MOTIFS/motifs_filtered_cns.bed",header = F)
colnames(x2) <- c("chrom","start","stop","gene","score","strand","motif")
Un_motif <- unique(x2$motif)
x_df_comb_cns <- lapply(1:length(Un_motif),function(i){
  message(paste(round(i/length(Un_motif),3)*100),"%")
  x_i <- x2[which(x2$motif==Un_motif[[i]]),]
  x <- data.frame(motif=Un_motif[[i]],
                  genes= length(x_i$gene),
                  unique_genes = length(unique(x_i$gene)),
                  FIMO=NA,
                  CB=NA,
  #threshold =thr,
                  method="COMBINED_FILTERED_CNS")  
  return(x)
  
})

x_df_comb_cns <- do.call(rbind,x_df_comb_cns)

x_df_comb_final <- rbind(x_df_comb,x_df_comb_cns[,-c(4,5)])


gghistogram(x_df_comb_final, x = "unique_genes",
            add = "median", rug = TRUE,
            color = "method",
            fill = "method", palette = c("#00AFBB", "#E7B800","green","red"),
            add_density = TRUE,bins = 200,xlab = "Unique genes")


my_comparisons2 <- list( c("Cluster buster", "FIMO"),c("Cluster buster", "COMBINED"),c("FIMO", "COMBINED"),c("COMBINED","COMBINED_FILTERED_CNS"))
ggviolin(x_df_comb_final, x = "method", y = "unique_genes", fill = "method",
         palette = c("#00AFBB","red","#E7B800","green"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)   

ggviolin(x_df_comb_final, x = "method", y = "unique_genes", fill = "method",
         palette = c("#00AFBB","red","#E7B800","green"),
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)

Un_motif <- unique(x2$motif)
feature_file_cns <- lapply(1:length(Un_motif),function(i){
  #i <- 1
  message(paste(round(i/length(Un_motif),3)*100),"%")
  x_i <- x2[which(x2$motif==Un_motif[[i]]),]
  x <- data.frame(MOTIF=Un_motif[[i]],
                  TARGET = unique(x_i$gene))  
  return(x)
  
})





feature_file_cns <- do.call(rbind,feature_file_cns)

write.table(feature_file_cns,"/scratch/bis_klpoe/chsos/data/MOTIFS/feature_comb_CNS_NH.tsv",na = "",row.names = F,sep = "\t",quote = F,col.names = F)

save.image("~/CLUST_COMB_CNS.RData")
load("~/CLUST_COMB_CNS.RData")
require(data.table)
x_C_len <- fread("/scratch/bis_klpoe/chsos/data/MOTIFS/feature_files/feature_comb_NH.tsv",header = F)
#unique genes for combined s
length(unique(x_C_len$V2))
length(unique(feature_file_cns$TARGET))
