
DEG_edgeR_func_comb <- function(combinefolders=FALSE,out_name=NULL,folder_name, pval,out_dir,met_dir,plot_MDS,numCores,group_vect){

  ##############################################################################
  #Loading libraries
  library(ggplot2);library(edgeR);
  library(tximport);library(parallel)
  ##############################################################################
  #checking that combinefolders is not null
  if(!is.null(combinefolders)){
    if(!combinefolders %in% c(TRUE,FALSE)){
      stop("combinefolders is invalid. It must be TRUE or FALSE")
    }
  } else{
    stop("combinefolders is invalid. It must be TRUE or FALSE")
  }
  
  
  ##############################################################################
  #checking if it combining several folders in one run
  if(combinefolders==T){
    if(length(folder_name)>1){
      message(paste("Processing multifolders:",folder_name,"\n"))
      message(paste("Using ",out_name," to name output folder"))
    } else {
      stop("combinefolders is activated but only one ")
    }
  } else {
      message(paste("Processing folder:",folder_name))
    }

  ##############################################################################
 
  # Detecting number of cores for parallelizing
  
  x_det <- NULL
  x_det <- parallel::detectCores()
  message(paste("Detecting cores, total cores are: ",x_det))
  ##############################################################################
  if (numCores > x_det) {
    stop("Number of cores exceed the maximum allowed by the machine,
         use a coherent number of cores such as four")
  }
  
  ##############################################################################
  #Creating output dir folder
  if(combinefolders==T){
    out_dir_1 <- paste0(out_dir,"/",out_name)
    if(!dir.exists(out_dir_1)){
      dir.create(paste0(out_dir,"/",out_name))
    }
  } else {
    out_dir_1 <- paste0(out_dir,"/",folder_name)
    if(!dir.exists(out_dir_1)){
      dir.create(paste0(out_dir,"/",folder_name))
    }
  }

  ##############################################################################
  # Defining plot dirs for contrasts
  plt_dir <- paste0(out_dir_1,"/graphics")
  if(!dir.exists(plt_dir)){
    dir.create(plt_dir)
  }
  
  # Defining csv dir for contrasts
  csv_dir <- paste0(out_dir_1,"/csv")
  if(!dir.exists(csv_dir)){
    dir.create(csv_dir)
  }
  
  ##############################################################################
  #Reading salmon files
  
  if(combinefolders==T){
    message(paste("Reading salmon files for",folder_name,"\n"))
  } else {
    message("Reading salmon files")
  }
  
  
  #reading and combining salmon_tx2gene.tsv files
  tx2gene <- list()
  for(i in 1:length(folder_name)){
    tx2gene[[i]] <- read.table(paste0(data_dir[[i]],"/","salmon_tx2gene.tsv"))
    };rm(i)
  tx2gene <- do.call(cbind,tx2gene)
  colnames(tx2gene) <- paste("V",1:ncol(tx2gene))

    #reading metadata
  #tx2gene
  meta <- list()
  for(i in 1:length(folder_name)){
    meta[[i]] <- read.csv(paste0(met_dir,"/","run",folder_name[[i]],".csv"), header = TRUE, sep=",")
    #meta[[i]]$folder_name <- folder_name[[i]]
  };rm(i)
  meta <- do.call(rbind,meta)
  meta$trt <- sapply(strsplit(meta$sample,"_"),"[[",1)
  
  
  #reading quant files
  list1 = list.files(path = data_dir,pattern ="quant.sf",recursive = T,full.names = T)
  #original approach
  #list1 = paste0(data_dir,"/", meta$sample,"/", "quant.sf")
  
  #read salmon tximport
  txi <- tximport::tximport(files=list1, type = "salmon", tx2gene = tx2gene)
  data <- as.data.frame(txi$counts)
  ab <-  as.data.frame(txi$abundance)
  df <- colSums(data)
  
  
  write.table(df, paste0(out_dir_1,"/","library_size.txt"), col.names = F, quote = F)
  
  ##############################################################################
  #Testing that groups have more than one sample. (this avoid NA dispersion values)
  trt_summary <- as.data.frame(tapply(meta$trt,meta$trt,length))
  trt_summary$group <- row.names(trt_summary)
  trt_summary <- trt_summary[,c(2,1)]
  colnames(trt_summary) <- c("group","count")
  
  
  #Using group file if it is available
  if(!is.null(group_vect)){
    if(length(meta$trt)!=length(group_vect)){
      stop("You are providing a group_vect object that not match with the number of samples.
           Please cheack and resubmit")
    } else {
      meta$trt <-   group_vect
    }
    
  } else {
    #If a group have only one sample the function stops!
    if(any(trt_summary$count==1)){
      stop("Each group must have at least two samples, please check and provide
           an object group_vect with the group and run again.")
    }
  }
  
  
  
  ##############################################################################
  message("Defining groups with the metadata provided and levels")
  
  # Defining treatments
  colnames(data) <- meta$sample
  colnames(ab) <- meta$sample
  trt <- factor(meta$trt)#,levels = levels)
  
  #Defining possible combinations
  cmb_un <- as.data.frame(t(combn(levels(trt),2)))
  colnames(cmb_un) <- c("SOURCE","TARGET")
  cmb_un$COMP <- paste0(cmb_un$TARGET, "-", cmb_un$SOURCE)
  
  message("Possible combinations available:",nrow(cmb_un))
  print(cmb_un)
  # Save normalized (but not filtered) CPM
  y <- edgeR::DGEList(counts = data,group = trt)
  y$samples$lib.size <- colSums(y$counts)
  y <- edgeR::calcNormFactors(y)
  
  
  message("Saving raw data obtained (tpm and cpm)")
  
  #saving cpm
  write.table(cpm(y), paste0(out_dir_1,"/","CPM_normalized.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  #save tpm
  write.table(ab, paste0(out_dir_1,"/","abundance_drought_tpm_genelevel.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  ##############################################################################
  # Filtering
  message("Filtering using cpm 2")
  keep <- rowSums(cpm(data)>1)>=2 
  data_kept <- data[keep,]
  ##############################################################################
  #normal factors
  message("Normalizing")
  # Normalization and MDS plot
  d <- edgeR::DGEList(data_kept,group =trt)
  d$samples$lib.size <- colSums(d$counts)
  d <- edgeR::calcNormFactors(d)
  
  ##############################################################################
  #Plotting MDS
  t <- limma::plotMDS(d,plot = F)
  MDS_df <- data.frame(t$x, t$y)
  
  #saving MDS with shapes (it plots numbers)
  
  MDS<- ggplot(data = MDS_df) +
    geom_point(aes(x = t.x, y =t.y,shape = meta$trt),size = 4)  +
    theme_bw() +
    #    aes(label = meta$sample)+ #allpoints) +
    theme(legend.box="horizontal") +
    xlab("Leading logFC dim1") +
    ylab("Leading logFC dim2") +
    ggtitle("") +
    scale_shape_manual(values = c(1:length(trt), letters, LETTERS, "0", "1"))+
    theme(panel.grid =element_blank())
  
  MDS <- MDS +
    guides(shape = guide_legend(title = "Groups"))
  
  ggsave(paste0(out_dir_1,"/","MDS.pdf"), MDS,height = 8, width = 8)
  
  if(plot_MDS==TRUE){
    plot(MDS)  
  }
  ##############################################################################
  #Refining design to create contrasts
  #design
  trt <- trt
  #using treatments and 0 to do easily
  design <- model.matrix( ~0 + trt )
  colnames( design ) <- levels( trt )
  
  ##############################################################################
  message("Estimate dispersion")
  #estimateDisp
  d <- edgeR::estimateDisp(y = d, design = design,robust = T)
  pdf(paste0(out_dir_1,"/","plotBCV.pdf"),height = 8, width = 8)
  plotBCV(d)
  dev.off() 
  
  ##############################################################################
  #contrasts
  message("Creating contrasts")
  
  message("Using automatical combinations, analyze carefully!")
  if(nrow(cmb_un)>1){
    contrasts <- list()
    for(i in 1:nrow(cmb_un)){
      contrasts[[i]] <- makeContrasts(contrasts =  cmb_un$COMP[[i]],
                                      levels=levels(trt))
    };rm(i)
  } else {
    contrasts <- list(makeContrasts( contrasts = cmb_un$COMP[[1]] ,
                                     levels=levels(trt)))
  }
  
  message(paste("Contrasts provided: ",length(contrasts)))
  ##############################################################################
  #obtaining DEG per contrast 
  message(paste0("Calculating DE genes in parallel, using ",numCores," cores"))
  
  if(length(contrasts)==1){
    warning("Only one contrast will be run")
  }
  
  if(length(contrasts)>0){
    
    #saving contrasts design
    x_cont <- do.call(cbind,contrasts)
    write.csv(x_cont, file= paste0(out_dir_1,"/","contrasts.csv"),na = "",row.names = T)
    
    
    #Using parallel approach
    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("contrasts","d","design",
                                          "pval","cmb_un","plt_dir","out_dir_1",
                                          "csv_dir"),envir=environment())
    
    cont_results <- parallel::parLapplyLB(cl,
                                          X = seq_len(length(contrasts)),
                                          fun = function (i){
                                            #glmQLFit
                                            fit <- edgeR::glmQLFit(y = d, design = design,contrast=contrasts[[i]]) 
                                            pdf(paste0(plt_dir,"/",cmb_un$COMP[[i]],"_BCV.pdf"),height = 8, width = 8)
                                            edgeR::plotQLDisp(fit)
                                            dev.off()
                                            
                                            #testing glmQFTest
                                            qlf.2vs1 <- edgeR::glmQLFTest(glmfit = fit,contrast = contrasts[[i]])
                                            qlf.2vs1$table$FDR <- p.adjust(qlf.2vs1$table$PValue,method = "BH")
                                            qlf.2vs1$table$status <- as.factor(sign(qlf.2vs1$table$logFC))
                                            #defining names for values 1 and -1 
                                            qlf.2vs1$table$status_name <- NA
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- cmb_un$TARGET[[i]]
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- cmb_un$SOURCE[[i]]
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- "UP"
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- "DOWN"
                                            #subsetting and saving
                                            qlf.2vs1_final <- qlf.2vs1$table
                                            qlf.2vs1_final2 <- qlf.2vs1_final[which(qlf.2vs1_final$FDR < pval),]
                                            
                                            message(paste0("Genes kept after filtering: ",nrow(qlf.2vs1_final2)), "printing groups...")
                                            print(tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            
                                            #saving in csv tables (full and filtered)
                                            write.csv(qlf.2vs1_final, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                   cmb_un$COMP[[i]],"_","pval_",
                                                                                   as.character(pval),"_","fulltable.csv"),
                                                      row.names = T,quote = F)
                                            write.csv(qlf.2vs1_final2, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                    cmb_un$COMP[[i]],"_","pval_",
                                                                                    as.character(pval),"_","filtered.csv"),
                                                      row.names = T,quote = F) 
                                            
                                            
                                            #plot significant genes
                                            
                                            pdf(paste0(plt_dir,"/","plotMD_glmQLFTest_",cmb_un$COMP[[i]],".pdf"),height = 8, width = 8)
                                            limma::plotMD(qlf.2vs1)
                                            #abline(h=c(-1, 1), col="blue")
                                            dev.off()
                                            
                                            
                                            q1 <- data.frame(count = tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            q1$groups <- row.names(q1)
                                            q1 <- q1[,c(2,1)]
                                            print(q1)
                                            #saving in csv tables (full and filtered)
                                            write.csv(q1, file= paste0(out_dir_1,"/","glmQLFTest_",
                                                                       cmb_un$COMP[[i]],"_","pval_",
                                                                       as.character(pval),"_","summary.csv"),
                                                      row.names = F)
                                            
                                            #Do not move it is need to return qlf.2vs1_final
                                            qlf.2vs1_final2
                                          })
    
    parallel::stopCluster(cl)
  } else {
    stop("no contrasts to run.") #NEW
    cont_results <- NULL # NEW
  }
  
  final_list <- list(contrasts = contrasts,
                     contrasts_results=cont_results)
  message("DONE!")
  return(final_list)
  
}

##############################################################################
#examples
#group_vect <- c("DT6H","CK0H","ABA6H","CK0H","DT6H","ABA6H","CK0H","DT6H","ABA6H")
#group_vect <- NULL
##############################################################################
#parameters
#folder_name <- c(1,2,3) #folder name
pval = 0.05 #p value for filtering
numCores <- 4 #number of cores to use in parallel
plot_MDS <- TRUE #if plot should be appears in R session
#groups if it not easy to get from the headers
##############################################################################
#SRP098756 #done
# folder_name <- c(1,2,3) #folder name
# group_vect <- c("CROWN_C","CROWN_C","CROWN_D","CROWN_D", #1 folder
#                 "CROWN_C","LEAVES_C","LEAVES_C","LEAVES_C", #2 folder
#                 "CROWN_D","LEAVES_D","LEAVES_D","LEAVES_D", #2 folder
#                 "ROOTS_C","ROOTS_C","ROOTS_C",#3 folder
#                 "ROOTS_D","ROOTS_D","ROOTS_D") #3 folder
# out_name <- "SRP098756"
# combinefolders <- T
##############################################################################
#SRP237462 # done
# folder_name <- 4 #folder name
# group_vect <- c("D0H","D6H","D0H","D6H","D0H","D6H","UNKNOWN","UNKNOWN") #4 folder
# combinefolders <- F
# out_name <- "SRP237462"
##############################################################################
#SRP356530 # done
# folder_name <- 5 #folder name
# group_vect <- c("LEAVES_D","LEAVES_D","LEAVES_D",
#                 "LEAVES_C","LEAVES_C","LEAVES_C") #5 folder
# combinefolders <- F
# out_name <- "SRP356530"
##############################################################################
#SRP257474 #
# folder_name <- 6 #folder name
# combinefolders <- F
# group_vect <- c("DT6H","CK0H","ABA6H",
#                 "CK0H","DT6H","ABA6H",
#                 "CK0H","DT6H","ABA6H") #6 hours
# out_name <- "SRP257474"
##############################################################################
#SRP072216 #done
folder_name <- c(7,8,9,10) #folder name
group_vect <- c("LEAVES_TABA","LEAVES_TABA","LEAVES_TABA","LEAVES_D","LEAVES_D", #7 folder
                "LEAVES_D","LEAVES_C","LEAVES_C","LEAVES_C", #8 folder
                "LEAVES_TABA","LEAVES_TABA","LEAVES_TABA","LEAVES_D","LEAVES_D", #9 folder
                "LEAVES_D","LEAVES_C","LEAVES_C","LEAVES_C") #10 folder

# group_vect <- c("WTLEAVES_TABA","WTLEAVES_TABA","WTLEAVES_TABA","WTLEAVES_D","WTLEAVES_D", #7 folder
#                 "WTLEAVES_D","WTLEAVES_C","WTLEAVES_C","WTLEAVES_C", #8 folder
#                 "TaPYLEAVES_TABA","TaPYLEAVES_TABA","TaPYLEAVES_TABA","TaPYLEAVES_D","TaPYLEAVES_D", #9 folder
#                 "TaPYLEAVES_D","TaPYLEAVES_C","TaPYLEAVES_C","TaPYLEAVES_C") #10 folder
out_name <- "SRP072216"
combinefolders <- T
##############################################################################
##############################################################################
##############################################################################
##############################################################################
#Directories
#where are the salmon files
data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
#defining output folder
out_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
#path where the nf core metadata is available
met_dir <- "/scratch/bis_klpoe/chsos/data/sample_files/DONE"
##############################################################################
x <- DEG_edgeR_func_comb(combinefolders,out_name,folder_name, pval,out_dir,met_dir,plot_MDS,numCores,group_vect)
