library(Seurat) #Necessary for this script to run in terminal.



#This file: follow https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html.



runRCTD<-function(spatialSample,numOfCores){
  # spatialSample<-"p16"
  
  cat("\nRunning RCTD on ",spatialSample," heart...\n",sep="")
  timeStart<-Sys.time()
  
  if(spatialSample=="p16"){
    dataSpatial<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_p16.rds") #63.8 MB.
  }else if(spatialSample=="spider"){
    dataSpatial<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_spider.rds") #78.3 MB.
  }
  dataSnRNASeq<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSnRNASeq.rds") #2.6 GB.
  
  
  
  #Prepare snRNASeq data.
  counts<-dataSnRNASeq@assays[["RNA"]]@layers[["counts"]]
  rownames(counts)<-rownames(dataSnRNASeq)
  colnames(counts)<-colnames(dataSnRNASeq)
  
  # temp1<-as.vector(dataSnRNASeq$celltype)
  # temp2<-as.character(dataSnRNASeq@active.ident)
  # identical(temp1,temp2) #TRUE is good.
  cell_types<-dataSnRNASeq@active.ident
  reference<-spacexr::Reference(counts,cell_types) #Warning: number of cells per cell type is 18689, larger than maximum allowable of 10000. Downsampling number of cells to: 10000.
  
  
  
  #Prepare spatial data.
  coords<-data.frame(cells=dataSpatial@images$slice1$centroids@cells,
                     spatial_coord_x=dataSpatial@images$slice1$centroids@coords[,1],
                     spatial_coord_y=dataSpatial@images$slice1$centroids@coords[,2])
  rownames(coords)<-coords$cells
  coords$cells<-NULL
  
  counts<-dataSpatial@assays[["Spatial"]]@layers[["counts"]]
  rownames(counts)<-rownames(dataSpatial)
  colnames(counts)<-colnames(dataSpatial)
  puck<-spacexr::SpatialRNA(coords,counts) #Create SpatialRNA object.
  
  
  
  #Run RCTD.
  myRCTD<-spacexr::create.RCTD(puck,reference,max_cores=numOfCores)
  myRCTD<-spacexr::run.RCTD(myRCTD,doublet_mode="full") #Full mode: no restrictions on number of cell types per spot.
  
  
  
  timeEnd<-Sys.time()
  print(timeEnd-timeStart) #Time difference of ~20 mins.
  
  
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.2_resultRCTD/",spatialSample,".rds")
  saveRDS(myRCTD,filename)
}



runRCTD(spatialSample="p16",numOfCores=8) #6 cores seems to be about as efficient.
runRCTD(spatialSample="spider",numOfCores=8) #6 cores seems to be about as efficient.



#Rscript 25.04.07.1.2_runRCTD.R > 1.2_resultRCTD/_log.txt 2> 1.2_resultRCTD/_log2.txt



result<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.2_resultRCTD/p16.rds")

weights<-as.matrix(result@results$weights)
range(rowSums(weights)) #Row sums are not all 1.

weights<-as.data.frame(spacexr::normalize_weights(result@results$weights)) #Make row sums 1.
range(rowSums(weights)) #Row sums are 1.









