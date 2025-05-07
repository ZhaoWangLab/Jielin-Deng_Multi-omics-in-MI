library(dplyr)
library(tidyr)
library(ggplot2)



plotHeatmap<-function(spatialSample){
  # spatialSample<-"p16"
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.1_tableOfWeights_",spatialSample,".rds")
  tableOfWeights<-readRDS(filename)
  
  #Convert Cell2location numbers to proportions within spots.
  if(TRUE){
    rowIndices<-which(tableOfWeights$method=="Cell2location")
    temp<-tableOfWeights[rowIndices,-c(1:2)]
    temp2<-temp/rowSums(temp)
    tableOfWeights[rowIndices,-c(1:2)]<-temp2
    # range(rowSums(tableOfWeights[,-c(1:2)])) #Row sums are all approximately 1.
  }
  
  #Filter for positive spots only.
  if(TRUE){
    filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/0.1_positiveSpots_",spatialSample,".csv")
    positiveSpots<-read.csv(filename)
    # sum(positiveSpots$Barcode%in%tableOfWeights$cellName)
    tableOfWeights<-tableOfWeights%>%filter(cellName%in%(positiveSpots$Barcode)) #1000*12.
  }
  
  #Prepare dataPlot.
  if(TRUE){
    tableOfWeights_melted<-tableOfWeights%>%reshape2::melt(id.vars=c("cellName","method"),
                                                           variable.name="cellType",
                                                           value.name="proportion")
    tableOfWeights_melted$method<-factor(tableOfWeights_melted$method,
                                         levels=unique(tableOfWeights_melted$method))
    
    dataPlot<-tableOfWeights_melted%>%group_by(method,cellType)%>%
      summarise(n=n(),mean=mean(proportion),sd=sd(proportion),.groups="drop")%>% #Use .groups="drop" to avoid warning message.
      mutate(se=sd/sqrt(n))
    
    dataPlot2<-dataPlot%>%select(method,cellType,mean)
    dataPlot3<-dataPlot2%>%spread(key=method,value=mean)
  }
  
  # cor(dataPlot3$RCTD,dataPlot3$Cell2location)
  # dataPlot3%>%ggplot(aes(x=RCTD,y=Cell2location))+
  #   geom_point()
  
  corMatrix<-cor(dataPlot3[,-1],dataPlot3[,-1])
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.3.2_corMatrix_",spatialSample,".rds")
  saveRDS(corMatrix,filename)
  
  corMatrixMelted<-reshape2::melt(corMatrix)
  colnames(corMatrixMelted)<-c("method1","method2","PCC")
  
  #Calculate average correlation with other methods.
  corMatrix2<-corMatrix
  diag(corMatrix2)<-0
  cat("Average correlation with the other methods is:\n")
  print(rowSums(corMatrix2)/(nrow(corMatrix2)-1))
  
  
  
  corMatrixMelted%>%ggplot(aes(x=method1,
                               y=method2,
                               fill=PCC,
                               color=PCC))+
    geom_tile()+
    coord_fixed()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          text=element_text(size=20),
          legend.text=element_text(size=12),
          legend.title=element_text(size=15),
          panel.background=element_blank(), #Panel background.
          panel.grid.major.y=element_blank(), #Horizontal grid lines.
          panel.grid.major.x=element_blank(), #Vertical grid lines.
          axis.text.x=element_text(angle=45,hjust=1,vjust=1), #See https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot.
    )
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.3_plots/heatmap_",
                   spatialSample,".pdf")
  ggsave(filename,width=15,height=13,unit="cm")
}



plotHeatmap(spatialSample="p16")
plotHeatmap(spatialSample="spider")









