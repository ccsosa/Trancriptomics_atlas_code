# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rtracklayer")
require("rtracklayer");library(parallel);require(ggpubr)
################################################################################
CBUSTER_dir <- "/scratch/bis_klpoe/indis/nf_motif_map_calibrate_2K_1K/output/CB_processed_outputs"
FIMO_dir <- "/scratch/bis_klpoe/indis/nf_motif_map_calibrate_2K_1K/output/FIMO_processed_outputs"
out_dir <- "/scratch/bis_klpoe/chsos/data/MOTIFS"
################################################################################
thr <- 0.8
x_C <- list.files(CBUSTER_dir,pattern = ".bed.gz")
numCores <- 3
################################################################################
#cl <- parallel::makeCluster(numCores)


# parallel::clusterExport(cl, varlist=c("CBUSTER_dir","x_C","thr","out_dir"),envir=environment())
# x_C_len <- parallel::parLapplyLB(cl,
                                 #X = seq_len(length(x_C)),
                                 #fun = function (i){
x_C_len <- lapply(1:length(x_C),function(i){
  message(paste0("i: ",i, "/ ", round((i/length(x_C)*100),2)," %"),"")
                                   #i <- 1
                                   if(!file.exists(paste0(out_dir,"/CLUSTER_BUSTER/",sub(".bed.gz",".bed",x_C[[i]])))){
                                      message("processing")
                                     xx <- rtracklayer::import(paste0(CBUSTER_dir,"/",x_C[[i]]),format = "BED")
                                     range_end_start <- end(xx) - start(xx)
                                     range_end_start2 <- range_end_start<1000
                                     xx2 <- xx[range_end_start2]
                                     xx2 <- xx2[order(xx2$score,decreasing = T)]
                                     if(length(xx2)>20000){
                                       
                                       xx2 <- xx2[1:20000]
                                     }
                                     
                                     xx2 <- xx2[order(GenomeInfoDb::seqnames(xx2),IRanges::ranges(xx2),decreasing = F)]
                                     rtracklayer::export(xx2,paste0(out_dir,"/CLUSTER_BUSTER/",sub(".bed.gz",".bed",x_C[[i]])),format = "BED")
                                     #xx2 <- xx2[which(xx2$score>thr)]
                                     # length(xx2)
                                     # length(unique(xx2$name))
                                     
                                   } else {
                                     message("loading processed results")
                                     xx2 <- rtracklayer::import(paste0(out_dir,"/CLUSTER_BUSTER/",sub(".bed.gz",".bed",x_C[[i]])),format = "BED")
                                   }
                                   
                                   x <- data.frame(motif=sub(".bed.gz","",x_C[[i]]),
                                                   genes= length(xx2),
                                                   unique_genes = length(unique(xx2$name)),
                                                   #threshold =thr,
                                                   method="Cluster buster")
                                   return(x)
                                 })

#parallel::stopCluster(cl)

x_C_len <- do.call(rbind,x_C_len)
save.image("~/CLUST20.RData")
#load("~/CLUST.RData")
################################################################################
x_F <- list.files(FIMO_dir,pattern = ".bed")

# #cl <- parallel::makeCluster(numCores)
# parallel::clusterExport(cl, varlist=c("FIMO_dir","x_F","thr"),envir=environment())
#x_F_len <- parallel::parLapplyLB(cl,
#                                  X = seq_len(length(x_F)),
#                                  fun = function (i){
x_F_len <- lapply(1:length(x_F),function(i){
  message(paste0("i: ",i, "/ ", round((i/length(x_F)*100),2)," %"),"")
  #i <- 1
  if(!file.exists(paste0(out_dir,"/FIMO/",sub(".bed.gz",".bed",x_C[[i]])))){
    message("processing")
  xx <- rtracklayer::import(paste0(FIMO_dir,"/",x_F[[i]]),format = "BED")
  range_end_start <- end(xx) - start(xx)
  range_end_start2 <- range_end_start<1000
  xx2 <- xx[range_end_start2]
  if(length(xx2)>0){
  xx2 <- xx2[order(xx2$score,decreasing = T)]
  if(length(xx2)>15000){
    
    xx2 <- xx2[1:15000]
  }
  
  xx2 <- xx2[order(GenomeInfoDb::seqnames(xx2),IRanges::ranges(xx2),decreasing = F)]
  rtracklayer::export(xx2,paste0(out_dir,"/FIMO/",sub(".bed.gz",".bed",x_F[[i]])),format = "BED")
  } else {
    message("No coords found")
  }
  #xx2 <- xx2[which(xx2$score>thr)]
  # length(xx2)
  # length(unique(xx2$name))
  
  } else {
    message("loading processed results")
    xx2 <- rtracklayer::import(paste0(out_dir,"/FIMO/",sub(".bed.gz",".bed",x_C[[i]])),format = "BED")
  }
  x <- data.frame(motif=sub(".bed","",x_F[[i]]),
                  genes= length(xx2),
                  unique_genes = length(unique(xx2$name)),
                  #threshold =thr,
                  method="FIMO")
  return(x)
  #x
})

# parallel::stopCluster(cl)

x_F_len <- do.call(rbind,x_F_len)

#save.image("~/CLUST_FIMO.RData")
################################################################################
x_df <- rbind(x_C_len,x_F_len)
save.image("~/CLUST_FIMO15.RData")
################################################################################
write.table(x_df,"length2015.tsv")
################################################################################
# Combine histogram and density plots
gghistogram(x_df, x = "unique_genes",
            add = "median", rug = TRUE,
            color = "method",
            fill = "method", palette = c("#00AFBB", "#E7B800"),
            add_density = TRUE,bins = 200,xlab = "Unique genes")



my_comparisons <- list( c("Cluster buster", "FIMO") )
save.image("~/CLUST_FIMO15.RData")
ggviolin(x_df, x = "method", y = "unique_genes", fill = "method",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)   
save.image("~/CLUST_FIMO15.RData")

################################################################################
#load 
load("~/CLUST_FIMO15.RData")
