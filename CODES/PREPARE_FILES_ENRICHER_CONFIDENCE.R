require(data.table)
################################################################################
#analysis dir
an_dir <- "/scratch/bis_klpoe/chsos/analysis"
#DEG dir
DEG_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
#path to save the inputs for  the enrichment analysis
in_dir <- "/scratch/bis_klpoe/chsos/analysis/MOTIF_ENRICHMENT/inputs"
#path to save the results for  the enrichment analysis
res_dir <- "/scratch/bis_klpoe/chsos/analysis/MOTIF_ENRICHMENT/results"
#path to save the bash file to run Enrichr
xx_dir <- "/scratch/bis_klpoe/chsos/analysis/MOTIF_ENRICHMENT"
################################################################################
#read core and pan genes
CORE_GENES <- as.data.frame(fread(paste0(DEG_dir,"/","LFC_CORE.tsv"),header = T))
PAN_GENES <- as.data.frame(fread(paste0(DEG_dir,"/","LFC_PAN.tsv"),header = T))

if(any(colnames(CORE_GENES) %in% "#")){
  row.names(CORE_GENES) <- CORE_GENES[,"#"]
  CORE_GENES[,"#"] <- NULL
}

if(any(colnames(PAN_GENES) %in% "#")){
  row.names(PAN_GENES) <- PAN_GENES[,"#"]
  PAN_GENES[,"#"] <- NULL
}
################################################################################
#printing new numbers for pan
tapply(PAN_GENES$RES_HC,PAN_GENES$RES_HC,length)
#printing new numbers for core
tapply(CORE_GENES$RES_HC,CORE_GENES$RES_HC,length)
################################################################################
###splitting files per categories
##apply for normal groups (silenced)
cats <- unique(PAN_GENES$RES_HC) 
cats <- cats[c(4,5,6,7)]
###############################################################################
#defining  lists to use to create the bash file to run the motif enrichment
core_list1 <- list()
pan_list1 <- list()
core_list2 <- list()
pan_list2 <- list()

################################################################################
#writing inputs for the motif enrichment analysis
for(i in seq_len(length(cats))){
 write.table(row.names(CORE_GENES[which(CORE_GENES$RES_HC==cats[i]),]),
             paste0(in_dir,"/","CORE_",cats[[i]],".tsv"),sep = "\t",na = "",quote = F,col.names = F,row.names = F)
 write.table(row.names(PAN_GENES[which(PAN_GENES$RES_HC==cats[i]),]),
              paste0(in_dir,"/","PAN_",cats[[i]],".tsv"),sep = "\t",na = "",quote = F,col.names = F,row.names = F)
  
  
  #no filtered by CNS
  core_list1[i] <- 
  paste("/group/transreg/Tools/dreec_enrichment/enricherv2.4 /scratch/bis_klpoe/chsos/data/MOTIFS/feature_files/feature_comb_NH.tsv",
  paste0(in_dir,"/","CORE_",cats[[i]],".tsv"), '"-m" "2" "-f" "0.05" "-o"',
  paste0(res_dir,"/","CORE_",cats[[i]],"_N_CNS.tsv"))
  #filtered by CNS
  core_list2[i] <- 
    paste("/group/transreg/Tools/dreec_enrichment/enricherv2.4 /scratch/bis_klpoe/chsos/data/MOTIFS/feature_files/feature_comb_CNS_NH.tsv",
          paste0(in_dir,"/","CORE_",cats[[i]],".tsv"), '"-m" "2" "-f" "0.05" "-o"',
          paste0(res_dir,"/","CORE_",cats[[i]],"_CNS.tsv"))

  #no filtered by CNS
  pan_list1[i] <- 
    paste("/group/transreg/Tools/dreec_enrichment/enricherv2.4 /scratch/bis_klpoe/chsos/data/MOTIFS/feature_files/feature_comb_NH.tsv",
          paste0(in_dir,"/","PAN_",cats[[i]],".tsv"), '"-m" "2" "-f" "0.05" "-o"',
          paste0(res_dir,"/","PAN_",cats[[i]],"_N_CNS.tsv"))
  #filtered by CNS
  pan_list2[i] <- 
    paste("/group/transreg/Tools/dreec_enrichment/enricherv2.4 /scratch/bis_klpoe/chsos/data/MOTIFS/feature_files/feature_comb_CNS_NH.tsv",
          paste0(in_dir,"/","PAN_",cats[[i]],".tsv"), '"-m" "2" "-f" "0.05" "-o"',
          paste0(res_dir,"/","PAN_",cats[[i]],"_CNS.tsv"))
  
}

################################################################################
#creating the bash file
final_list <- c("#!/bin/bash",
                "# -*- coding: utf-8 -*-",
                core_list1,
                core_list2,
                pan_list1,
                pan_list2)
###run for normal categories
#saving the bash file
write.table(unlist(final_list),
            paste0(xx_dir,"/","MOTIF_ENRICHMENT.sh"),sep = "\t",na = "",quote = F,col.names = F,row.names = F)


################################################################################
