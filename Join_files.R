require(dplyr);require(data.table);require(pheatmap);library(caret)
#Directories
#where are the salmon files
################################################################################
#function to save heatmaps in pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
################################################################################
#get folder names to join them
folders <- c(1,2,3,5,6,7,8,9,10)
folder_name <- folders[[1]]
################################################################################
#listing folders to get files 
data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
vhRR_dir <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input"
#Read salmon outputs to build the atlas
x1 <- read.table(paste0(data_dir,"/","salmon.merged.gene_tpm.tsv"),header = T)
x1$tx <- NULL
x1$gene_name <- NULL
print(paste(nrow(x1)," / ",1))
################################################################################
#join files
for(i in 2:length(folders)){
  message(i)
  folder_name <- folders[[i]]
  data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
  x <- read.table(paste0(data_dir,"/","salmon.merged.gene_tpm.tsv"),header = T)
  x$tx <- NULL
  x1$gene_name <- NULL

  print(identical(x$gene_id,x1$gene_id))
  x$gene_id <- NULL

  x1 <- cbind(x1,x)
  print(paste(nrow(x1)," / ",i))
  rm(x)
};rm(i)

x1$gene_name <- NULL
################################################################################
#writing all samples in one unique atlas
row.names(x1) <- x1$gene_id
################################################################################
#remove problematic samples
colnames(x1)
#removing 8 samples
x1 <- x1[,-c(25:33)]
#write drought transcriptomic atlas
write.table(x1,"/scratch/bis_klpoe/chsos/analysis/atlas_drought_total.tsv",sep = "\t",na = "",row.names = F,col.names =T,quote = F)
write.table(x1,paste0(vhRR_dir,"/atlas_drought_total.tsv"),sep = "\t",na = "",row.names = F,col.names = T,quote = F)
write.table(x1[,-1],paste0(vhRR_dir,"/atlas_drought_total_NC.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)

################################################################################
################################################################################
################################################################################

#subsetting contrasts to see
sets = c("SRP356530/LEAVES_D-LEAVES_C", 
         "SRP098756/CROWN_D-CROWN_C",
         "SRP098756/LEAVES_D-LEAVES_C",
         "SRP098756/ROOTS_D-ROOTS_C",
         ##"SRP098756/ROOTS_D-LEAVES_D",
         #"SRP257474/CK0H-ABA6H",
         #"SRP257474/DT6H-ABA6H",         
         "SRP072216/LEAVES_D-LEAVES_C",
         "SRP072216/LEAVES_TABA-LEAVES_D")
################################################################################
#adding p value to filter
pval <- 0.05
data_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"         
sets_st_cont <- strsplit(sets,"/")
################################################################################
##loading all DEG for all contrasts
x_deg <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_filtered.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg <- do.call(rbind,x_deg)

################################################################################
#getting all genes accross all studies
x_unique_genes <- as.data.frame(matrix(nrow=length(unique(x_deg$V1)),ncol=length(sets)+1))
x_unique_genes$V1 <- unique(x_deg$V1)
colnames(x_unique_genes) <- c("genes",sets)
row.names(x_unique_genes) <- unique(x_deg$V1)

pb <-
  utils::txtProgressBar(min = 0,
                        max = nrow(x_unique_genes),
                        style = 3)

#Presence of a gene in a set of contrasts
for(i in 1:nrow(x_unique_genes)){
  utils::setTxtProgressBar(pb, i)
  x_i <- x_deg[which(x_deg$V1==x_unique_genes$genes[[i]]),]
  x_unique_genes[i,-1] <- sets %in% x_i$group # transforming boolean in 1 and 0s
  x_unique_genes[i,-1] <- (x_unique_genes[i,-1]*1)
};rm(i)

close(pb)
################################################################################
#unique gene to get pan subset genes!
x_unique_genes$total <- rowSums(x_unique_genes[,-1])
#subsetting
unique_genes <- unique(x_unique_genes$genes)
subset_pan <- x1[which(x1$gene_id %in% unique_genes),]
write.table(subset_pan,paste0(vhRR_dir,"/atlas_drought_pan.tsv"),sep = "\t",na = "",row.names = F,col.names = F,quote = F)
################################################################################
################################################################################
#defining core genes with four or more contrasts
#x_unique_genes_dif_contrasts <- x_unique_genes
#get unique genes by contrast at least dos
#x_unique_genes_dif_contrasts <- x_unique_genes_dif_contrasts[which(x_unique_genes_dif_contrasts$total>1),]
#core genes 5 contrasts
#all_contr_genes2 <- x_unique_genes$genes[which(x_unique_genes$total>=5)]
#core genes 4 contrasts
all_contr_genes3 <- x_unique_genes$genes[which(x_unique_genes$total>=4)]
########################################################################################

########################################################################################
#plotting core subset

subset2 <- x1[which(x1$gene_id %in% all_contr_genes3),]; subset2a <- subset2
subset2 <- subset2[rowSums(subset2[,-1]) > 0,]; subset2a <- subset2

row.names(subset2) <- subset2$gene_id; subset2$gene_id <- NULL;subset2$gene_name <- NULL
#write.table(subset1,"/scratch/bis_klpoe/chsos/analysis/atlas_drought_419.tsv",sep = "\t",na = "",row.names = T,col.names = F)
write.table(subset2,paste0(vhRR_dir,"/atlas_drought_261.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)

# #min max
process <- preProcess(subset2, method=c("range"))
norm_scale <- predict(process, subset2)

# z-score
# preproc1 <- preProcess(subset2, method=c("center", "scale"))
# norm_scale <- predict(preproc1, subset2)

Breaks <- seq(floor(min(norm_scale)),ceiling(max(norm_scale)), by = 0.2);Breaks[length(Breaks)] <- max(norm_scale) #by=1
#Breaks <-  round(Breaks,1)
Breaks2 <- as.character(Breaks); Breaks2[length(Breaks2)] <- "Min-Max Normalization (TPM)\n"#"Z-Score Normalization(TPM)\n"


#subset1[subset1==0] <- NA
p <- pheatmap::pheatmap(norm_scale,cluster_rows = T,
                        cluster_cols = F,
                        fontsize = 18,
                        annotation_legend = F,#F,
                        labels_row =NULL,
                        #annotation_col = ann_col,
                        # annotation_colors =  list(genelist=c(DOWNREGULATED="blue",#"#006EC9",
                        #                                      UPREGULATED="red")),#"#FF00FF")),
                        drop_levels=F,
                        # display_numbers = matrix_top_go_final_2,#T,#TRU1E,
                        #number_color = "white",
                        fontsize_row = 0.01,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        fontsize_number=14,#10
                        cellwidth = 40, #48
                        angle_col = 45,
                        cellheight = 13, #16
                        #pvalue
                        legend_breaks = Breaks,
                        legend_labels = Breaks2,
                        #color
                        #color = hcl.colors(16, "Greens"), #Oranges
                        color = hcl.colors(10, "BluYl"),
                        #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
                        legend = T,
                        annotation_names_col = F,
                        na_col = "gray",
                        #number_format= "%.2f",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_4_or_more_contrasts.pdf"),
                  width = 50,
                  height = 52 )
rm(norm_scale,preproc1)
########################################################################################
########################################################################################
########################################################################################
########################################################################################
#Preparing LFC heatmap

x_deg_total <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_fulltable.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg_total <- do.call(rbind,x_deg_total)


########################################################################################
########################################################################################
#419 genes
subset2a <- subset2
subset2a <- x_deg_total[x_deg_total$V1 %in% row.names(subset2a),]
LFC_2a <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subset2a$V1))))
LFC_2a[,1] <- unique(subset2a$V1)
colnames(LFC_2a) <- c("gene",sets)
LFC_2a_p_val <- LFC_2a

for(i in 2:ncol(LFC_2a)){
  
  x_i <- subset2a[which(subset2a$group==sets[i-1]),]
  for(j in 1:nrow(LFC_2a)){
    LFC_2a[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_2a$gene[[j]]),2])
    LFC_2a_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_2a$gene[[j]]),6])
  };rm(j)
};rm(i)

for(i in 2:ncol(LFC_2a)){
  LFC_2a_p_val[,i] <- as.numeric(LFC_2a[,i])
  LFC_2a[,i] <- as.numeric(LFC_2a_p_val[,i])
  LFC_2a_p_val[,i][which(LFC_2a_p_val[,i]>=pval)] <- ""
  LFC_2a_p_val[,i][which(LFC_2a_p_val[,i]<pval)] <- "*"
  LFC_2a_p_val[,i][which(is.na(LFC_2a_p_val[,i]))] <- ""
  LFC_2a[,i][which(is.na(LFC_2a[,i]))] <- 0
};rm(i)




row.names(LFC_2a) <- LFC_2a$gene;
row.names(LFC_2a_p_val) <- LFC_2a$gene;
LFC_2a$gene <- NULL;LFC_2a_p_val$gene <- NULL

#LFC_2a$`SRP257474/CK0H-ABA6H` <- -LFC_2a$`SRP257474/CK0H-ABA6H`
LFC_2a$`SRP072216/LEAVES_TABA-LEAVES_D` <- -LFC_2a$`SRP072216/LEAVES_TABA-LEAVES_D`

#colnames(LFC_2a)[5] <-  "SRP257474/ABA6H-CK0H"
colnames(LFC_2a)[6] <-  "SRP072216/LEAVES_D-LEAVES_TABA"

write.table(LFC_2a,paste0(data_dir,"/LFC_CORE.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)


ann_col <- data.frame(Tissue=c("leaves","crown","leaves","roots","leaves","leaves"))
row.names(ann_col) <- colnames(LFC_2a)#sets
p <- pheatmap::pheatmap(LFC_2a,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Tissue=c(leaves="purple",
                                                          crown="brown",
                                                          roots= "red"#,
                                                          #three_leaf="gray"
                        )),
                        drop_levels=F,
                        fontsize_row = 0.0001,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 40, #48
                        angle_col = 45,
                        cellheight = 4, #16
                        legend_breaks = c(-7.5,-5,-2.5,0,2.5,5,7.5,
                                          max(LFC_2a,na.rm = T)),
                        legend_labels = c("-7.5","-5","-2.5","0","2.5","5","7.5","log2 fold change\n"), 
                        color = hcl.colors(8, "BrBg"),#hcl.colors(6, "BluYl"),
                        legend = T,
                        annotation_names_col = T,
                        na_col = "gray",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_4_or_more_contrasts_LFC.pdf"),
                  width = 20,
                  height = 27)

################################################################################
#LFC for all pan core genes
subsetall <- as.data.frame(subset_pan)
#subsetting only to DEG
subsetall <- as.data.frame(x_deg_total[x_deg_total$V1 %in% row.names(subsetall),])
#fitting a matrix with n unique genes as rows 
LFC_all <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subsetall$V1))))
#get unique genes as first column

LFC_all[,1] <- unique(subsetall$V1)
colnames(LFC_all) <- c("gene",sets)

#copying matrix for filtering steps
LFC_all_p_val <- LFC_all
LFC_all_p_val2 <- LFC_all

#filling out the matrix
for(i in 2:ncol(LFC_all)){
  message(i)
  x_i <- subsetall[which(subsetall$group==sets[i-1]),]
  for(j in 1:nrow(LFC_all)){
    #message(j)
    if(length(as.numeric(x_i[which(x_i$V1==LFC_all$gene[[j]]),2]))>0){
      #filling out with the LFC value
      LFC_all[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_all$gene[[j]]),2])
      #filling out with p values
      LFC_all_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_all$gene[[j]]),6]) 
      #if NA fill out matrices with NAs
      if(is.na(LFC_all_p_val[j,i])){
        LFC_all_p_val2[j,i] <- NA
      } else if(LFC_all_p_val[j,i]>=pval){
        LFC_all_p_val2[j,i] <- NA
      } else if(LFC_all_p_val[j,i]<pval){
        #if p value is < 0,05 cell is 1
        LFC_all_p_val2[j,i] <- 1
      }
    } else {
      LFC_all[j,i] <- NA
      LFC_all_p_val[j,i] <- NA 
      LFC_all_p_val2[j,i] <- NA
    }

  };rm(j)
};rm(i)


row.names(LFC_all) <- LFC_all$gene;row.names(LFC_all_p_val2) <- row.names(LFC_all) 
LFC_all$gene <- NULL;LFC_all_p_val$gene <- NULL;LFC_all_p_val2$gene <- NULL

#invert order to fix the problem
LFC_all$`SRP072216/LEAVES_TABA-LEAVES_D` <- -LFC_all$`SRP072216/LEAVES_TABA-LEAVES_D`

#if a p value is < 0.05 the value is conserved but it is >0.05 is NA. This is to get majority rule

LFC_all_p_val2 <- LFC_all_p_val2*LFC_all
colnames(LFC_all)[6] <-  "SRP072216/LEAVES_D-LEAVES_TABA"
colnames(LFC_all_p_val2)[6] <-  "SRP072216/LEAVES_D-LEAVES_TABA"

#writing the tables
write.table(LFC_all,paste0(data_dir,"/LFC_PAN.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)
write.table(LFC_all_p_val2,paste0(data_dir,"/LFC_PAN_fixed.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)

################################################################################
#get majority rule up and down datasets
#PAN
#obtain LFC values sign to obtain majority rule for DEG with p val < 0,.05
LFC_all2 <- sign(LFC_all_p_val2)#sign(LFC_all[,c(1,2,3,4,6,7,8)])

#adding extra columns
LFC_all2$POS_LEN <- NA
LFC_all2$NEG_LEN <- NA
LFC_all2$RES <- NA
LFC_all2$RES_HC <- NA

#filling out matrix for the majority rule

for(i in 1:nrow(LFC_all2)){
  #counting negative and positive LFC values 
   xpos <- LFC_all2[i,-c(7,8,9,10)][which(LFC_all2[i,-c(7,8,9,10)]==1)]
   xneg <- LFC_all2[i,-c(7,8,9,10)][which(LFC_all2[i,-c(7,8,9,10)]==-1)]
  
   #counting values with LFC positive values 
  if(length(xpos)>0){
    LFC_all2$POS_LEN[[i]] <-    abs(sum(xpos,na.rm = T))
  } else {
    LFC_all2$POS_LEN[[i]] <-    0
  }
   #counting values with LFC negative values 
  if(length(xneg)>0){
    LFC_all2$NEG_LEN[[i]] <-    abs(sum(xneg,na.rm = T))
  } else {
    LFC_all2$NEG_LEN[[i]] <-    0
  }
  
   #approach to get the majority rule
  x_sign <- which.max(c(LFC_all2$POS_LEN[[i]],LFC_all2$NEG_LEN[[i]]))
    
    
    #adding suggested status 
    if(LFC_all2$POS_LEN[[i]]==LFC_all2$NEG_LEN[[i]]){
      LFC_all2$RES[[i]] <- "GENE"
    }else if(x_sign==2){
      LFC_all2$RES[[i]] <- "DOWN"
    } else if(x_sign==1){
      LFC_all2$RES[[i]] <- "UP"
    }
    
  #adding high (HC), medium (MC) and low (LC) confidence values for upregulated DEG
  #Upregulated HC: at least five positive LFC values and no negative LFC
  #Upregulated MC: at least three positive LFC values and no negative LFC 
  #Upregulated LC: at least one positive LFC values and no negative LFC
  
  if(LFC_all2$POS_LEN[[i]]>= 5 & 
     #(ncol(LFC_all2[-c(7:10)]))-1 & 
     LFC_all2$NEG_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "UP_HC"
  } else if((LFC_all2$POS_LEN[[i]]>=3 & LFC_all2$POS_LEN[[i]] <5) &
            LFC_all2$NEG_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "UP_MC"
  } else if((LFC_all2$POS_LEN[[i]]>=1 & LFC_all2$POS_LEN[[i]] <3) &
            LFC_all2$NEG_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "UP_LC"
  }
  #adding high (HC), medium (MC) and low (LC) confidence values for downregulated DEG
  #Downregulated HC: at least five negative LFC values and no positive LFC
  #Downregulated MC: at least three negative LFC values and no positive LFC
  #DownregulatedLC: at least one negative LFC values and no positive LFC
  
  if(LFC_all2$NEG_LEN[[i]]>= 5 & 
     #(ncol(LFC_all2[-c(7:10)]))-1 & 
     LFC_all2$POS_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "DOWN_HC"
  } else if((LFC_all2$NEG_LEN[[i]]>=3 & LFC_all2$NEG_LEN[[i]] <5) &
            LFC_all2$POS_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "DOWN_MC"
  } else if((LFC_all2$NEG_LEN[[i]]>=1 & LFC_all2$NEG_LEN[[i]] <3) &
            LFC_all2$POS_LEN[[i]]==0){
    LFC_all2$RES_HC[[i]] <- "DOWN_LC"
  }
  
  if(is.na(LFC_all2$RES_HC[[i]])){
    LFC_all2$RES_HC[[i]] <- "VARIABLE"
  }  
  
};rm(i)
##printing up,down and tie genes 
tapply(LFC_all2$RES,LFC_all2$RES,length)
#sum(tapply(LFC_all2$RES,LFC_all2$RES,length))
##printing up and down high, medium and low confidence (HC,MC,LC respectively)
tapply(LFC_all2$RES_HC,LFC_all2$RES_HC,length)
#sum(tapply(LFC_all2$RES_HC,LFC_all2$RES_HC,length))
################################################################################
#LOST GENES
 sum(tapply(LFC_all2$RES,LFC_all2$RES,length))-
 sum(tapply(LFC_all2$RES_HC,LFC_all2$RES_HC,length))
################################################################################
#CORE subset genes for majority rule
LFC_all2_core <- LFC_all2
LFC_all2_core <- LFC_all2_core[row.names(subset2),]

#printing up,down and tie genes 
tapply(LFC_all2_core$RES,LFC_all2_core$RES,length)
#printing up and down high, medium and low confidence (HC,MC,LC respectively)
tapply(LFC_all2_core$RES_HC,LFC_all2_core$RES_HC,length)
################################################################################
#saving files for pan and core gene subsets
write.table(LFC_all2,paste0(data_dir,"/LFC_PAN.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)
write.table(LFC_all2_core,paste0(data_dir,"/LFC_CORE.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)

