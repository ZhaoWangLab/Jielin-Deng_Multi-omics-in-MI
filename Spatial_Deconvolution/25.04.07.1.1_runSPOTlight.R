library(SingleCellExperiment) #Necessary for metadata(dec).



#This file: follow https://www.bioconductor.org/packages/release/bioc/vignettes/SPOTlight/inst/doc/SPOTlight_kidney.html.



runSPOTlight<-function(spatialSample){
  # spatialSample<-"p16"
  
  cat("\nRunning SPOTlight on ",spatialSample," heart...\n",sep="")
  timeStart<-Sys.time()
  
  if(spatialSample=="p16"){
    dataSpatial<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_p16.rds") #63.8 MB.
  }else if(spatialSample=="spider"){
    dataSpatial<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_spider.rds") #78.3 MB.
  }
  dataSnRNASeq<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSnRNASeq.rds") #2.6 GB.
  
  
  
  spe<-Seurat::as.SingleCellExperiment(dataSpatial) #Warning messages: 1: Layer ‘data’ is empty. 2: Layer ‘scale.data’ is empty.
  dim(spe) #19465  1808.
  
  if(spatialSample=="p16"){
    dataSnRNASeqSub<-subset(dataSnRNASeq,sample=="MI_4")
  }else if(spatialSample=="spider"){
    dataSnRNASeqSub<-subset(dataSnRNASeq,sample=="MI_1")
  }
  sce<-Seurat::as.SingleCellExperiment(dataSnRNASeqSub)
  dim(sce) #32285 8518.
  #sce@assays@data@listData[["counts"]]: 1, 1, 1, 1, 3, 2...
  #sce@assays@data@listData[["logcounts"]]: 1.38, 1.38, 1.38, 1.38, 2.29, 1.94...
  
  
  
  sce<-scater::logNormCounts(sce)
  #sce@assays@data@listData[["counts"]]: 1, 1, 1, 1, 3, 2...
  #sce@assays@data@listData[["logcounts"]]: 1.01, 1.01, 1.01, 1.01, 2.02, 1.60...
  
  #Variance modeling.
  genes<-!grepl(pattern="^Rp[l|s]|Mt",x=rownames(sce)) #Get vector indicating which genes are neither ribosomal or mitochondrial.
  # which(genes==FALSE)
  
  dec<-scran::modelGeneVar(sce,subset.row=genes)
  plot(dec$mean,dec$total,xlab="Mean log-expression",ylab="Variance")
  curve(metadata(dec)$trend(x),col="blue",add=TRUE)
  
  hvg<-scran::getTopHVGs(dec,n=3000) #Get the top 3000 genes.
  
  #Obtain the marker genes for each cell identity.
  colLabels(sce)<-colData(sce)$celltype
  mgs<-scran::scoreMarkers(sce,subset.row=genes) #Compute marker genes. This takes a minute or so.
  
  #Keep only those genes that are relevant for each cell identity.
  mgs_fil<-lapply(names(mgs),function(i){
    x<-mgs[[i]]
    x<-x[x$mean.AUC>0.8,] #Filter and keep relevant marker genes, those with AUC>0.8.
    x<-x[order(x$mean.AUC,decreasing=TRUE),] #Sort the genes from highest to lowest weight.
    x$gene< rownames(x) #Add gene and cluster id to the dataframe.
    x$cluster<-i
    data.frame(x)
  })
  mgs_df<-do.call(rbind,mgs_fil)
  
  
  
  #Cell Downsampling.
  idx<-split(seq(ncol(sce)),sce$celltype) #Split cell indices by identity.
  n_cells<-100
  cs_keep<-lapply(idx,function(i){
    n<-length(i)
    if(n<n_cells)
      n_cells<-n
    sample(i,n_cells)
  })
  sce<-sce[,unlist(cs_keep)]
  dim(sce) #32285  866.
  
  
  
  #Deconvolution.
  mgs_df$gene<-rownames(mgs_df) #This is necessary. Otherwise: ids %in% names(mgs) are not all TRUE.
  res<-SPOTlight::SPOTlight(
    x=sce,
    y=spe,
    groups=as.character(sce$celltype),
    mgs=mgs_df,
    hvg=hvg,
    weight_id="mean.AUC",
    group_id="cluster",
    gene_id="gene")
  
  
  
  timeEnd<-Sys.time()
  print(timeEnd-timeStart) #Time difference of ~4 mins.
  
  
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.1_resultSPOTlight/",spatialSample,".rds")
  saveRDS(res,filename)
}



results<-parallel::mclapply(c("p16","spider"),FUN=function(spatialSample){
  # spatialSample<-"p16"
  
  runSPOTlight(spatialSample=spatialSample)
  
  return(0)
},mc.cores=2)
cat("\n") #Include a new line here to emphasize the print message.
print(unlist(results)) #0 means success.



#Rscript 25.04.07.1.1_runSPOTlight.R > 1.1_resultSPOTlight/_log.txt 2> 1.1_resultSPOTlight/_log2.txt



result<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.1_resultSPOTlight/p16.rds")
temp<-result$mat
range(rowSums(temp)) #Row sums are 1.








