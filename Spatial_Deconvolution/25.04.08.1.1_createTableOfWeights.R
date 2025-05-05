library(dplyr)
library(tidyr)



cellTypesOrdered<-c("Fibroblasts","Macrophages","Epicardial cells","Cardiomyocytes","Endothelial cells",
                    "Pericytes","Smooth muscle cells","Immune-T cells","Immune-B cells","Dentritic cells")
#Spelling of "dendritic cells" is corrected in this file.



createTableOfWeights<-function(spatialSample){
  # spatialSample<-"p16"
  
  #Process SPOTlight result.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.1_resultSPOTlight/",spatialSample,".rds")
    result<-readRDS(filename)
    resultSPOTlight<-as.data.frame(result$mat) #1808*10.
    range(rowSums(resultSPOTlight)) #Row sums are 1.
    
    resultSPOTlight$cellName<-rownames(resultSPOTlight) #1808*11.
    rownames(resultSPOTlight)<-NULL
    resultSPOTlight$method<-"SPOTlight" #1808*12.
    resultSPOTlight<-resultSPOTlight%>%select(cellName,method,all_of(cellTypesOrdered)) #1808*12.
  }
  
  #Process RCTD result.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.2_resultRCTD/",spatialSample,".rds")
    result<-readRDS(filename)
    resultRCTD<-as.data.frame(spacexr::normalize_weights(result@results$weights)) #1808*10.
    range(rowSums(resultRCTD)) #Row sums are 1.
    
    resultRCTD$cellName<-rownames(resultRCTD) #1808*11.
    rownames(resultRCTD)<-NULL
    resultRCTD$method<-"RCTD" #1808*12.
    resultRCTD<-resultRCTD%>%select(cellName,method,all_of(cellTypesOrdered)) #1808*12.
  }
  
  #Process CellDART result.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.3_resultCellDART/",spatialSample,"/cellfraction.csv")
    resultCellDART<-readr::read_csv(filename)
    colnames(resultCellDART)[1]<-"cellName"
    range(rowSums(resultCellDART[,-1])) #Row sums are approximately 1.
    
    #Edit column names.
    temp<-colnames(resultCellDART)[-1]
    temp2<-substr(temp,1,nchar(temp)-6) #Remove suffixes.
    colnames(resultCellDART)[-1]<-temp2
    
    resultCellDART$method<-"CellDART" #1808*12.
    resultCellDART<-resultCellDART%>%select(cellName,method,all_of(cellTypesOrdered)) #1808*12.
  }

  #Process Tangram result.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.4_resultTangram/",spatialSample,".csv")
    resultTangram<-readr::read_csv(filename)
    colnames(resultTangram)[1]<-"cellName"
    range(rowSums(resultTangram[,-1])) #Row sums are between 0 and 140.
    
    #Convert to proportions within spot (see Chunk 9 plot in https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb).
    temp<-resultTangram[,-1]
    temp2<-temp/rowSums(temp)
    resultTangram[,-1]<-temp2
    
    resultTangram$method<-"Tangram" #1808*12.
    resultTangram<-resultTangram%>%select(cellName,method,all_of(cellTypesOrdered)) #1808*12.
  }
  
  #Process cell2location result.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/1.5_resultCell2location/",spatialSample,".csv")
    resultCell2location<-readr::read_csv(filename)
    colnames(resultCell2location)[1]<-"cellName"
    range(rowSums(resultCell2location[,-1])) #Row sums are between 0 and 40.
    
    #Edit column names.
    temp<-colnames(resultCell2location)[-1]
    temp2<-substr(temp,24,nchar(temp)) #Remove prefixes.
    colnames(resultCell2location)[-1]<-temp2
    
    resultCell2location$method<-"cell2location" #1808*12.
    resultCell2location<-resultCell2location%>%select(cellName,method,all_of(cellTypesOrdered)) #1808*12.
  }
  
  tableOfWeights<-rbind(resultSPOTlight,resultRCTD,resultCellDART,resultTangram,resultCell2location)
  tableOfWeights<-tableOfWeights%>%rename("Dendritic cells"="Dentritic cells") #Correct spelling.
  
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.1_tableOfWeights_",spatialSample,".rds")
  saveRDS(tableOfWeights,filename)
}



createTableOfWeights(spatialSample="p16")
createTableOfWeights(spatialSample="spider")



#For future use:
spatialSample<-"p16"
filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.1_tableOfWeights_",spatialSample,".rds")
tableOfWeights<-readRDS(filename)








