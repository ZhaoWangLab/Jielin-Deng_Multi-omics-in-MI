library(dplyr)
library(tidyr)
library(ggplot2)


cellTypesOrdered<-c("Fibroblasts","Macrophages","Epicardial cells","Cardiomyocytes","Endothelial cells",
                    "Pericytes","Smooth muscle cells","Immune-T cells","Immune-B cells","Dendritic cells")


plotMap<-function(spatialSample,methodCurr,indexOfMethod,scale){
  # spatialSample<-"p16"
  # methodCurr<-"cell2location"
  # indexOfMethod<-5
  # scale<-"commonScale"
  
  #Prepare tableOfWeightsSub.
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.1_tableOfWeights_",spatialSample,".rds")
  tableOfWeights<-readRDS(filename)
  tableOfWeightsSub<-tableOfWeights%>%filter(method==methodCurr)
  upperLimit<-max(tableOfWeightsSub[,-c(1:2)])
  upperLimit<-ceiling(upperLimit/0.1)*0.1 #Round up to the nearest 0.1. Changing to 0.2 doesn't make a difference in the plots.
  
  #Prepare dataSpatial.
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_",spatialSample,".rds")
  dataSpatial<-readRDS(filename)
  metaData<-dataSpatial@meta.data
  metaData$cellName<-rownames(metaData)
  
  metaData2<-metaData%>%left_join(tableOfWeightsSub,by="cellName")
  rownames(metaData2)<-metaData2$cellName #Necessary for plotting to work properly.
  dataSpatial@meta.data<-metaData2
  
  # Seurat::SpatialFeaturePlot(dataSpatial,features="Fibroblasts")+theme(legend.position="right") #Check plotting.
  
  plotList<-vector(mode="list",length=length(cellTypesOrdered))
  
  for(i in 1:length(plotList)){
    # i<-1
    
    cellType<-cellTypesOrdered[i]
    p<-Seurat::SpatialFeaturePlot(dataSpatial,features=cellType,
                                  image.scale="hires",image.alpha=0
    )+
      labs(title=cellType)+
      theme(legend.position="right",
            plot.title=element_text(hjust=0.5),
            legend.title=element_blank(),
            text=element_text(size=24)
      )
    
    if(spatialSample=="p16"){
      p<-p+theme(aspect.ratio=1.2)
    }else if(spatialSample=="spider"){
      p<-p+theme(aspect.ratio=1.7)
    }
    
    if(scale=="individualScales"){
      limits<-NULL
    }else if(scale=="commonScale"){
      limits<-c(0,upperLimit)
    }
    
    if(methodCurr=="SPOTlight"){
      p<-p+scico::scale_fill_scico(limits=limits,palette="oslo") #Blue.
    }else if(methodCurr=="RCTD"){
      p<-p+scico::scale_fill_scico(limits=limits,palette="lajolla") #Yellow.
    }else if(methodCurr=="CellDART"){
      p<-p+scale_fill_viridis_c(limits=limits,option="viridis") #Green.
    }else if(methodCurr=="Tangram"){
      p<-p+scale_fill_viridis_c(limits=limits,option="cividis") #Light yellow.
    }else if(methodCurr=="cell2location"){
      p<-p+scale_fill_viridis_c(limits=limits,option="magma") #Red.
    }
    
    plotList[[i]]<-p
  } #End of for(i in 1:length(plotList))
  
  
  
  pCombined<-ggpubr::ggarrange(plotlist=plotList,nrow=2,ncol=5)
  
  filename<-paste0("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_spatialCellTypes/9_downstreamAnalysis/1.2_plots/",
                   scale,"_",spatialSample,"_",indexOfMethod,"_",methodCurr,".pdf")
  if(spatialSample=="p16"){
    ggsave(filename,width=60,height=25,unit="cm")
  }else if(spatialSample=="spider"){
    ggsave(filename,width=60,height=30,unit="cm")
  }
}


spatialSamples<-c("p16","spider")
methods<-c("SPOTlight","RCTD","CellDART","Tangram","cell2location")
scales<-c("commonScale","individualScales")



for(spatialSample in spatialSamples){
  for(indexOfMethod in 1:length(methods)){
    # indexOfMethod<-1 #For naming the pdf files.
    
    for(scale in scales){
      plotMap(spatialSample=spatialSample,methodCurr=methods[indexOfMethod],indexOfMethod,scale=scale)
    }
  }
}










