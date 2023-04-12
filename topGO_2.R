topGO_gen <- function(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background){
  
  ##############################################################################
  #loading reducereduntGO_ALL.R function to reduce redundancy in enrichments using DAG
  source("/scratch/bis_klpoe/chsos/analysis/script/reducereduntGO_ALL.R")  
  ##############################################################################
  # #calling libraries  
  require(data.table);require(topGO);require(dplyr);require(gProfileR);require(GO.db)
  #require(parallel);require(biomaRt);
  
  ##############################################################################
  ##reading goa file downloaded from plaza
  message("Reading GOA file")
  ta_go_file <- fread(ta_go_file)
  ######################
  if(custom_format==F){
    #changing "#gene_id" colname to "gene_id"
    #  colnames(ta_go_file)[1] <- unlist(strsplit(colnames(ta_go_file)[1],"#"))[[2]]
    #subsetting to reduce size
    ta_go_file <- ta_go_file[,c(2,1)]#.(gene_id,go,description)]
    colnames(ta_go_file) <- c("gene_id","go")
  } else {
    #  colnames(ta_go_file)<- c("go","gene_id","na")
    #subsetting to reduce size
    ta_go_file <- ta_go_file[,c(2,1)]#.(gene_id,go)]  
    colnames(ta_go_file) <- c("gene_id","go")
  }
  
  #indexing background genes to subset and do not use the complete genome
  index_b <- ta_go_file$gene_id %in% unique(ta_go_file$`gene_id`)
  x_back_vect <- ta_go_file[index_b,]
  
  
  if(is.null(background)){
    message("NO BACKGROUND ADDED. USING FILE ANNOTATION UNIQUE GENES PROVIDED IN THE ta_go_file")
  } else {
    x_back_vect <- x_back_vect[x_back_vect$gene_id %in% background,] 
  }
  ##############################################################################
  #if results are written, omitting run again!
  if(!file.exists(paste0(out_dir,"/",output_name,".tsv"))){
    
    #########
    #GO:MF  
    #########
    
    message("Running GO:MF enrichment using gprofiler2")
    x_s2 <-  gprofiler2::gost(query = genes_id,
                              organism = "taestivum", ordered_query = FALSE,
                              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                              measure_underrepresentation = FALSE, evcodes = FALSE,
                              user_threshold = pval, correction_method = "g_SCS",
                              domain_scope = "annotated", custom_bg = unique(x_back_vect$gene_id),
                              numeric_ns = "", sources = "GO:MF", as_short_link = FALSE)
    if(!is.null(x_s2)){
      all_res1_MF <- reducereduntGOALL(df=x_s2$result,
                                       GO.ID="term_id",FDR_col="p_value",mode="MF")
      if(nrow(all_res1_MF)>0){
        all_res1_MF$TYPE <- "GO:MF:GPROFILER"
      } else {
        all_res1_MF <- NULL
      }
      
    } else {
      message("No significant enriched results!")
      all_res1_MF <- NULL
    }
    
    ##############################################################################
    message("Running GO:BP using plaza files and topGO")
    
    
    #creating a list of gene ids and their GO terms     
    gene_2_GO <- unstack(x_back_vect[,c(2,1)])
    
    message("Formatting to topGO format")
    #remover genes sin anotacion
    keep <- genes_id %in% unique(x_back_vect$`gene_id`)
    keep <- which(keep==TRUE)
    candidate_list <- genes_id[keep]
    #converting gene list to factors 
    geneList <- factor(as.integer(unique(x_back_vect$gene_id) %in% candidate_list),levels = c(0,1))
    names(geneList) <- unique(x_back_vect$`gene_id`)
    
    message("Running topGO")
    #Create a topGO object
    GOdata <- new('topGOdata',
                  ontology='BP', 
                  allGenes = geneList, 
                  #annot = annFUN.GO2genes,#
                  annot = topGO::annFUN.gene2GO,
                  gene2GO = gene_2_GO)
    
    ##new (get DAG level for each GO term)
    # x_level <- buildLevels(GOdata@graph)
    # x_level_l <- data.frame(GO=names(as.list(x_level$nodes2level)),
    #                         LEVEL=as.numeric(as.list(x_level$nodes2level)))
    
    #Run statistical tests
    resultFisher <- topGO::runTest(GOdata,algorithm = "classic",
                                   statistic = "fisher")
    
    #returning all GO terms used
    allGO1 <- topGO::usedGO(GOdata)
    allGO1a <- allGO1
    #allGO1a <-allGO1[allGO1 %in% ta_go_file[,2]]
    
    
    message("Summarizing  topGO results")
    
    #summarizing results in one table
    all_res1 <- topGO::GenTable(GOdata, classicFisher=resultFisher,
                                orderBy='classicFisher', topNodes=length(allGO1))
    
    #getting scores to avoid mistakes in GenTable
    pValue.classic <- topGO::score(resultFisher)
    #creating a data.frame to join with all_res1
    pValue.classic <- data.frame(GO.ID=names(pValue.classic),p_val=pValue.classic)
    #joining accurate p values
    all_res1 <- dplyr::left_join(all_res1,pValue.classic,c("GO.ID"))
    #fitting to goa file
    ###############################
    message("Pre-filtering step: keeping only GO Terms used in the GOA input file")
    all_res1 <-all_res1[all_res1$GO.ID %in% allGO1a,]
    ###############################
    #adjust p value using fisher test
    all_res1$FDR <- stats::p.adjust(all_res1$p_val,method = "BH")
    
    #getting -log FDR for a Manhattan plot
    #all_res1$logpv <- -log10(all_res1$FDR)
    #creating categories significant and no significant
    all_res1$status_name <- NA
    all_res1$status_name[which(all_res1$FDR<pval)] <- "significant"
    all_res1$status_name[which(all_res1$FDR>=pval)] <- "no significant"
    #removing classicFisher column (this is replaced by p_val )
    all_res1$classicFisher <- NULL
    
    message("Subsetting only significant results")
    #subsetting only significant
    all_res1_s <- all_res1[which(all_res1$status_name=="significant"),]
    
    
    #Only apply filtering steps if there are enriched results!
    
    if(nrow(all_res1_s)>0){
      all_res1_BP <- reducereduntGOALL(df=all_res1_s,
                                       GO.ID="GO.ID",FDR_col="FDR",mode="BP")
      
      all_res1_BP <- all_res1_BP[order(all_res1_BP$FDR,decreasing = F),]
      all_res1_BP$TYPE <- "GO:BP:TopGO"
      
    } else {
      message("No significant enriched results for BP!")
      all_res1_BP <- NULL
    }
    
    ############################################################################
    #joining GO:MF and GO:BP in one file
    
    
    if(!is.null(all_res1_MF) & !is.null(all_res1_BP)){
      message("joining GO:MF and GO:BP in one file and saving results")
      
      all_res1_MF<- all_res1_MF[,c("term_id","term_name","p_value","TYPE")]
      all_res1_BP<- all_res1_BP[,c("GO.ID","Term","FDR","TYPE")]
      colnames(all_res1_BP) <- colnames(all_res1_MF)
      all_res1_s <- rbind(all_res1_MF,all_res1_BP)
      write.table(all_res1_s,paste0(out_dir,"/",output_name,".tsv"),row.names = F,na = "",sep = "\t")
      
    } else if(!is.null(all_res1_MF) & is.null(all_res1_BP)){
      message("Saving results only for GO:MF")
      
      all_res1_s<- all_res1_MF[,c("term_id","term_name","p_value","TYPE")]
      write.table(all_res1_s,paste0(out_dir,"/",output_name,".tsv"),row.names = F,na = "",sep = "\t")
      
    } else if(is.null(all_res1_MF) & !is.null(all_res1_BP)){
      message("Saving results only for GO:BP")
      
      all_res1_s<- all_res1_BP[,c("GO.ID","Term","FDR","TYPE")]
      colnames(all_res1_s) <- c("term_id","term_name","p_value","TYPE")
      write.table(all_res1_s,paste0(out_dir,"/",output_name,".tsv"),row.names = F,na = "",sep = "\t")
      
    } else if(is.null(all_res1_MF) & is.null(all_res1_BP)){
      
      warning("NO RESULTS TO BE REPORTED!")
      
    }
  } else {
    message("Already process, loading resulting if exists")
    all_res1_s <- fread(paste0(out_dir,"/",output_name,".tsv"),header=T)
  }
  message("DONE!")
  return(all_res1_s)
  
}

############################################################################
#getting DEG file
DEG_file <- data.table::fread("/scratch/bis_klpoe/chsos/analysis/DEG/LFC_PAN.tsv")

CATS <- unique(DEG_file$RES_HC)

top_results <- list()
for(i in 1:length(CATS)){
  #outdir
  out_dir <- "/scratch/bis_klpoe/chsos/analysis/GO_DEG/SUBSETS"
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  #file for GO:BP
  ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
  #format
  custom_format =F #format
  #p val 
  pval <- 0.05
  #getting gene ids
  genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_pan.tsv"
  genes_id_file <- data.table::fread(genes_id_file)
  background <- genes_id_file$V1

  genes_id <- DEG_file$`#`[which(DEG_file$RES_HC==CATS[[i]])]
  output_name <- paste0("GO_PAN_",CATS[[i]]) #outcome filename
  top_results[[i]] <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background)
};rm(i)

############################################################################
#getting DEG file
DEG_file <- data.table::fread("/scratch/bis_klpoe/chsos/analysis/DEG/LFC_CORE.tsv")

CATS <- unique(DEG_file$RES_HC)

top_results2 <- list()
for(i in 1:length(CATS)){
  #outdir
  out_dir <- "/scratch/bis_klpoe/chsos/analysis/GO_DEG/SUBSETS"
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  #file for GO:BP
  ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
  #format
  custom_format =F #format
  #p val 
  pval <- 0.05
  #getting gene ids
  genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_pan.tsv"
  genes_id_file <- data.table::fread(genes_id_file)
  background <- genes_id_file$V1
  
  genes_id <- DEG_file$`#`[which(DEG_file$RES_HC==CATS[[i]])]
  output_name <- paste0("GO_CORE_",CATS[[i]]) #outcome filename
  top_results[[i]] <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background)
};rm(i)

############################################################################
ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
#format
custom_format =F #format
#outdir
out_dir <- "/scratch/bis_klpoe/chsos/analysis/GO_DEG/SUBSETS"
#p val 
pval <- 0.05
#getting gene ids
genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_pan.tsv"
genes_id_file <- data.table::fread(genes_id_file)
genes_id <- genes_id_file$V1
output_name <- "PAN"

x <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background=NULL)


############################################################################
##CORE 

DEG_file <- data.table::fread("/scratch/bis_klpoe/chsos/analysis/DEG/LFC_CORE.tsv")

ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
#format
custom_format =F #format
#outdir
out_dir <- "/scratch/bis_klpoe/chsos/analysis/GO_DEG/SUBSETS"
#p val 
pval <- 0.05
#getting gene ids
genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_pan.tsv"
genes_id_file <- data.table::fread(genes_id_file)
background <- genes_id_file$V1
#getting gene ids
genes_id <- DEG_file$`#`

########################################
#NULL BACKGROUND
output_name <- "CORE_BNULL"
x <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background=NULL)
########################################
#USING PAN AS BACKGROUND
ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
output_name <- "CORE"
x <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name,background=background)

# go_gene_file <- data.table::fread(paste0(out_dir,"/","go_gene.tsv"),header = F)
# go_gene_file_sub <- go_gene_file[go_gene_file$V2 %in% genes_id]
# #write.table(go_gene_file_sub,paste0(out_dir,"/","go_gene_52.tsv"),row.names = F,na = "",sep = "\t")
# #write.table(go_gene_file_sub,paste0(out_dir,"/","go_gene_419.tsv"),row.names = F,na = "",sep = "\t",col.names = F)
# write.table(go_gene_file_sub,paste0(out_dir,"/","go_gene_pan.tsv"),row.names = F,na = "",sep = "\t",col.names = F)
# 
# ##
# #subset GO2GENE
# 
