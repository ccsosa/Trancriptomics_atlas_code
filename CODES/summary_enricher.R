require(data.table);require(ggpubr);require(dplyr)
################################################################################
data_dir <- "/scratch/bis_klpoe/chsos/analysis/MOTIF_ENRICHMENT/results"
################################################################################
################################################################################
################################################################################
#LOADING CORE
CORE_FILES_f <- list.files(data_dir,pattern = "CORE")
#CNS_FILES <- CNS_FILES[-c(2,3,6)]

CORE_FILES <- lapply(seq_len(length(CORE_FILES_f)), function(i){
  message(i)
  #i <- 4
  x <- as.data.frame(fread(paste0(data_dir,"/",CORE_FILES_f[[i]]),skip = 6))
  x <- x[order(x$enr_fold,decreasing = T),]
  if(nrow(x)==0){
    x[1,] <- NA
  }
  x$dataset <- sub(".tsv*","",CORE_FILES_f[[i]])
  
  #x$group <- "CNE + MOTIFS INTERSECT"
  return(x)
})

CORE_FILES <- do.call(rbind,CORE_FILES)
unique(CORE_FILES$dataset)

################################################################################
#LOADING PAN
PAN_FILES_f <- list.files(data_dir,pattern = "PAN")
#CNS_FILES <- CNS_FILES[-c(2,3,6)]

PAN_FILES <- lapply(seq_len(length(PAN_FILES_f)), function(i){
  message(i)
  #i <- 4
  x <- as.data.frame(fread(paste0(data_dir,"/",PAN_FILES_f[[i]]),skip = 6))
  x <- x[order(x$enr_fold,decreasing = T),]
  if(nrow(x)==0){
    x[1,] <- NA
  }
  x$dataset <- sub(".tsv*","",PAN_FILES_f[[i]])
  
  #x$group <- "CNE + MOTIFS INTERSECT"
  return(x)
})

PAN_FILES <- do.call(rbind,PAN_FILES)
unique(PAN_FILES$dataset)
`################################################################################
#LOADING GO_2001318
GO_2001318_FILES_f <- list.files(data_dir,pattern = "GO_2001318")

GO_2001318_FILES <- lapply(seq_len(length(GO_2001318_FILES_f)), function(i){
  message(i)
  #i <- 4
  x <- as.data.frame(fread(paste0(data_dir,"/",GO_2001318_FILES_f[[i]]),skip = 6))
  x <- x[order(x$enr_fold,decreasing = T),]
  if(nrow(x)==0){
    x[1,] <- NA
  }
  x$dataset <- sub(".tsv*","",GO_2001318_FILES_f[[i]])
    return(x)
})

GO_2001318_FILES <- do.call(rbind,GO_2001318_FILES)
unique(GO_2001318_FILES$dataset)
################################################################################
#LOADING GO_2001319
GO_2001319_FILES_f <- list.files(data_dir,pattern = "GO_2001319")

GO_2001319_FILES <- lapply(seq_len(length(GO_2001319_FILES_f)), function(i){
  message(i)
  #i <- 4
  x <- as.data.frame(fread(paste0(data_dir,"/",GO_2001319_FILES_f[[i]]),skip = 6))
  x <- x[order(x$enr_fold,decreasing = T),]
  if(nrow(x)==0){
    x[1,] <- NA
  }
  x$dataset <- sub(".tsv*","",GO_2001319_FILES_f[[i]])
  return(x)
})

GO_2001319_FILES <- do.call(rbind,GO_2001319_FILES)
unique(GO_2001319_FILES$dataset)
################################################################################
#join all
final_blend_file <- rbind(CORE_FILES,PAN_FILES,GO_2001318_FILES,GO_2001319_FILES)
#final_blend_file
final_blend_file2 <- final_blend_file[which(!is.na(final_blend_file$`q-val`<0.05)),]
final_blend_file2 <- final_blend_file[which(final_blend_file2$`q-val`<0.05),]
final_blend_file2 <- final_blend_file[which(final_blend_file2$enr_fold>=1.25),]
#writing combined files in a summary file with q value < 0.05 and enc fold >1.25
write.table(final_blend_file2,paste0(data_dir,"/","summary_file.tsv"),sep = "\t",row.names = F)

################################################################################
# JOIN GENES TO DF
#load files to complement data
mot_info_dir <- "/scratch/bis_klpoe/indis/motif_info"
CISBP <- fread(paste0(mot_info_dir,"/",
             "CISBP_TF_Information_all_motifs_plus_tae_available_motifs_with_new_ID_homoeolog_expanded.tsv"))

#Motif_ID, Family_Name, new, homoeolog_expanded
CISBP <- CISBP[,c("Motif_ID","Family_Name","new","homoeolog_expanded")]


#joining complete information
final_blend_file2a <- left_join(final_blend_file2,CISBP,c("ftr_id"="Motif_ID"))
write.table(final_blend_file2a,paste0(data_dir,"/","summary.tsv"),sep = "\t",row.names = F)
################################################################################
# JOIN GENES TO DF
mot_info_dir <- "/scratch/bis_klpoe/indis/motif_info"
JASPAR <- fread(paste0(mot_info_dir,"/",
                      "JASPAR_2020_motif_to_gene_species.tsv"))

#Motif_ID, Family_Name, new, homoeolog_expanded
JASPAR <- JASPAR[,c("ID","Name","Species","Class","Family")]

################################################################################
#filling out the file
file_list <- list()
for(i in 1:nrow(final_blend_file2)){
  #i <- 1
  x_i <- final_blend_file2[i,]
  ############################################################################
  x_i2 <- left_join(x_i,CISBP,c("ftr_id"="Motif_ID"))
  ############################################################################
  x_i2 <- x_i2[which(!is.na(x_i2$new)),]
  x_i2 <- x_i2[which(x_i2$new!=""),]
  ############################################################################
  x_i3 <- left_join(x_i,JASPAR,c("ftr_id"="ID"))
  ############################################################################
  if(length(x_i2)==0){ #is.na(x_i2$new[[1]]) | 
    #x_i$Family_Name <- unique(x_i2$Family_Name)
    x_i$Family_Name <- NA
    x_i$New <- NA
    x_i$New_num <- NA
    x_i$homoeolog_expanded <- NA
    x_i$homoeolog_expanded_num <- NA
    
  } else {
    if(length(x_i2$Family_Name)>0){
      x_i$Family_Name <- unique(x_i2$Family_Name)
    } else {
      x_i$Family_Name <- NA
    }

    x_i$New <- paste(x_i2$new,collapse = "//")
    x_i$New_num <- length(unique(x_i2$new))
    x_i$homoeolog_expanded <- paste(x_i2$homoeolog_expanded,collapse = "//")
    x_i$homoeolog_expanded_num <- length(unique(x_i2$homoeolog_expanded))
  }
    ############################################################################
    if(length(x_i3)==0){ #is.na(x_i2$new[[1]]) | 
      x_i$Family <- NA
      x_i3$Name <- NA
      x_i$Species <- NA
      x_i$Class <- NA
    } else {
      if(length(x_i3$Name)>0){
        x_i$Name <- unique(x_i3$Name)
      } else {
        x_i$Name <- NA
      }  
      x_i$Class <- paste(x_i3$Class,collapse = "//")
      x_i$Species <- paste(x_i3$Species,collapse = "//")
      x_i$Family <- x_i3$Family
  }
  
  x_i$family_final <- NA
  if(!is.na(x_i$Family_Name) & is.na(x_i$Family)){
    x_i$family_final <- x_i$Family_Name
  } else if(is.na(x_i$Family_Name) & !is.na(x_i$Family)) {
    x_i$family_final <- x_i$Family    
  }
  
  file_list[[i]] <- x_i
  rm(x_i,x_i2)
};rm(i)

file_list <- do.call(rbind,file_list)
file_list$datasets_obtained <- NA
file_list2 <- file_list

x_unique <- unique(file_list$dataset)
file_list2 <- file_list2[file_list2$dataset %in% x_unique[c(1:10,11,14,15)],]
x_unique <- unique(file_list2$dataset)



dict <- data.frame(id=x_unique,desc=NA)
dict$desc[[1]] <- "Downregulated core genes (filtered motifs)"
dict$desc[[2]] <- "Downregulated core genes"
dict$desc[[3]] <- "D-U core genes (filtered motifs)"
dict$desc[[4]] <- "Upregulated core genes  (filtered motifs)"
dict$desc[[5]] <- "Upregulated core genes"
dict$desc[[6]] <- "Downregulated pan genes (filtered motifs)"
dict$desc[[7]] <- "D-U pan genes (filtered motifs)"
dict$desc[[8]] <- "D-U pan genes"
dict$desc[[9]] <- "Upregulated pan genes  (filtered motifs)"
dict$desc[[10]] <- "Upregulated pan genes"
dict$desc[[11]] <- "pan candidate genes"
dict$desc[[12]] <- "pan candidate genes (filtered motifs)"
dict$desc[[13]] <- "core candidate genes"
write.table(file_list,paste0(data_dir,"/","summary_file_joined.tsv"),sep = "\t",row.names = F,na = "")
file_list2$family_final[which(file_list2$family_final=="")] <- NA
file_list2$family_final[which(file_list2$family_final=="Myb/SANT domain factors" )] <- "Myb/SANT"

unique(file_list2$family_final)
for(i in 1:nrow(dict)){
  file_list2$datasets_obtained[which(file_list2$dataset==dict$id[i])] <-  dict$desc[i]
  
};rm(i)
write.table(file_list2,paste0(data_dir,"/","summary_file_joined.tsv"),sep = "\t",row.names = F,na = "")

##################################################################
require(UpSetR);require(ComplexHeatmap)
x_unique2 <- unique(file_list2$datasets_obtained)
Mt_list <- list()
for(i in 1:length(x_unique2)){
#  i <- 1
  x_i <- unique(na.omit(file_list2$family_final[which(file_list2$datasets_obtained==x_unique2[[i]])]))
  x_i <- x_i[which(x_i!="")]
  Mt_list[[i]] <- x_i
}
names(Mt_list) <- x_unique2
x_up_list <- list_to_matrix(Mt_list)
x_up_list <- x_up_list[order(rowSums(x_up_list),decreasing = T),]
x_up_list <- x_up_list[c(1:15),]
m2 = make_comb_mat(x_up_list)
#m2 <- t(m2)


x_plot <- UpSet(m2,comb_order = order(comb_size(m2),decreasing = T))

pdf(paste0(data_dir,"/","UPSET_TF.pdf"),height = 20, width = 40)
x_plot
dev.off()

################################################################################