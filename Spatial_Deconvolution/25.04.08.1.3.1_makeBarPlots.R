library(dplyr)
library(tidyr)
library(ggplot2)



makeBarPlot<-function(spatialSample){
  # spatialSample<-"p16"
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.1_tableOfWeights_",spatialSample,".rds")
  tableOfWeights<-readRDS(filename)
  
  #Convert cell2location numbers to proportions within spots.
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
    
    dataPlot$cellTypeIndex<-match(dataPlot$cellType,unique(dataPlot$cellType))
    dataPlot$cellTypeIndex_method<-paste0(dataPlot$cellTypeIndex,"_",dataPlot$method)
    dataPlot<-dataPlot%>%arrange(cellTypeIndex,method)
    dataPlot$cellTypeIndex_method<-factor(dataPlot$cellTypeIndex_method,
                                         levels=unique(dataPlot$cellTypeIndex_method))
  }
  
  
  
  dataPlot%>%ggplot(aes(x=cellTypeIndex_method,y=mean,fill=cellType))+
    geom_col(width=0.7,position=position_dodge(width=0.8))+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),
                  width=0.7,
                  position=position_dodge(0.8),
                  linewidth=0.5)+
    labs(y="Proportion within spot")+
    scale_x_discrete(labels=dataPlot$method)+
    scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))+
    theme(axis.title.x=element_blank(),
          legend.title=element_blank(),
          text=element_text(size=36),
          panel.background=element_blank(), #Panel background.
          panel.grid.major.y=element_blank(), #Horizontal grid lines.
          panel.grid.major.x=element_blank(), #Vertical grid lines.
          axis.line=element_line(linewidth=0.75), #Axis lines.
          plot.margin=margin(t=2,r=1.5,b=1.5,l=2,unit="lines"), #Plot margin (default is 5.5 points for t, r, b, and l).
          axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16), #See https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot.
          legend.key.height=unit(2.2,"lines"), #Legend key height.
          legend.key.width=unit(2,"lines"), #Legend key width.
    )
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.3_plots/barplot_",
                   spatialSample,".pdf")
  ggsave(filename,width=60,height=25,unit="cm")
}



makeBarPlot(spatialSample="p16")
makeBarPlot(spatialSample="spider")









