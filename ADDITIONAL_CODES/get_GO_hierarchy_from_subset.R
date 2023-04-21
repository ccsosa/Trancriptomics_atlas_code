
# 
do_custom_files <- function(ta_go_file,
                            custom_format,
                            genes_id_file_pan,
                            genes_id_core,
                            out_dir,
                            output_name,
                            go_gene_file,
                            custom_go_pan,
                            custom_go_core){
  
  #calling libraries
  require(data.table);require(topGO);require(dplyr);
  require(GO.db);require(parallel);require(igraph);library(visNetwork)
  
  message("Reading GOA file")
  #reading goa file downloaded from plaza
  ta_go_file <- fread(ta_go_file)
  
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

  x <- ta_go_file[,c("go","gene_id")]
  
  x1 <- data.frame(go=custom_go_pan,gene_id=genes_id_pan)
  x2 <- data.frame(go=custom_go_core,gene_id=genes_id_core)
   
  x <- rbind(x,x1)
  x <- rbind(x,x2)
  
  write.table(x,paste0(out_dir,"/","go_gene_mod.tsv"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  
}


    
    #     
#     #indexing background genes to subset and do not use the complete genome
#     index_b <- ta_go_file$gene_id %in% genes_id
#     x_back_vect <- ta_go_file[index_b,]
#     
#     #creating a list of gene ids and their GO terms     
#     gene_2_GO <- unstack(x_back_vect[,c(2,1)])
#     #G2g <- inverseList(gene_2_GO) #convert from gene2GO to GO2gene
#     go2genes <- annFUN.gene2GO(whichOnto = "BP", gene2GO = gene_2_GO)
#     go2genes_un_list <- unique(names(go2genes))
#     message("Formatting to topGO format")
#     #remover genes sin anotacion
#     keep <- genes_id %in% unique(ta_go_file$`gene_id,go`)
#     keep <- which(keep==TRUE)
#     candidate_list <- genes_id[keep]
#     #converting gene list to factors 
#     geneList <- factor(as.integer(unique(ta_go_file$`gene_id,go`) %in% candidate_list),levels = c(0,1))
#     names(geneList) <- unique(ta_go_file$`gene_id,go`)
#     
#     message("Running topGO")
#     #Create a topGO object
#     GOdata <- new('topGOdata',
#                   ontology='BP', 
#                   allGenes = geneList, 
#                   #annot = annFUN.GO2genes,#
#                   annot = topGO::annFUN.gene2GO,
#                   gene2GO = gene_2_GO)
#     #new (get DAG level for each GO term)
#     x_level <- buildLevels(GOdata@graph)
#     x_level_l <- data.frame(GO=names(as.list(x_level$nodes2level)),
#                             LEVEL=as.numeric(as.list(x_level$nodes2level)))
#     
#     
#     x_level_l$original <- NA
#     x_level_l$keep <- NA
#     for(i in 1:length(go2genes_un_list)){
#       
#       s <- any(x_level_l$GO==go2genes_un_list[[i]])
#       if(isTRUE(s)){
#         x_level_l$original[which(x_level_l$GO==go2genes_un_list[[i]])] <- TRUE 
#       } else {
#         x_level_l$original[which(x_level_l$GO==go2genes_un_list[[i]])] <- FALSE 
#       }
#     };rm(i)
#       
#     x_level_l$keep[ x_level_l$LEVEL %in% c(max(x_level_l$LEVEL)-2, max(x_level_l$LEVEL))]   <- TRUE
#     x_level_l$keep[ !x_level_l$LEVEL %in% c(max(x_level_l$LEVEL)-2, max(x_level_l$LEVEL))]   <- FALSE
#     
#       go2genes_un_list
#     
#     x_level_l$original[ go2genes_un_list %in% x_level_l$GO ] <- TRUE
#     
#     x_level_l[which()]
#     go2genes_un_list
#     
#     g2 <- graph_from_graphnel(GOdata@graph)
#     igv <- visIgraph(igraph = g2, layout = "layout_with_sugiyama")
#     visNodes(igv, shape = "circle", value = 30)
#     
#     
#     return(all_res1_s)
#   

# as.list(GOBPCHILDREN[go2genes_un_list[[1]]])
# as.list(GOBPPARENTS[go2genes_un_list[[1]]])
#
ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
custom_format =F #format


custom_go_pan <- "GO:2001318"
custom_go_core <- "GO:2001319"

#getting gene ids
genes_id_file_pan <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_pan.tsv"
genes_id_file_pan <- data.table::fread(genes_id_file_pan)
genes_id_pan <- genes_id_file_pan$V1

genes_id_file_core <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_419.tsv"
genes_id_file_core <- data.table::fread(genes_id_file_core)
genes_id_core <- genes_id_file_core$V1

out_dir <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input" #outpur dir

x <-        do_custom_files(ta_go_file,
                            custom_format,
                            genes_id_file_pan,
                            genes_id_core,
                            out_dir,
                            output_name,
                            go_gene_file,
                            custom_go_pan,
                            custom_go_core)