library(dplyr)
library(ggplot2)


data<-Seurat::Load10X_Spatial(data.dir="/labs/zhaowang/Seq/231006_IGC-JD-21455/spidersample/outs/",
                              filename="filtered_feature_bc_matrix.h5") #167.6 MB. Takes about 10 seconds.
dim(data) #19465  3177. 19465 genes, 3177 "cells".

Seurat::SpatialDimPlot(data,pt.size.factor=0.5)


#Filter cells. From Didi.
if(TRUE){
  data[["percent.mt"]]<-Seurat::PercentageFeatureSet(data,
                                                     pattern="^mt-") #Calculate the proportion of transcripts mapping to mitochondrial genes per cell.
  
  metaData<-data@meta.data #Each row is a cell.
  unique(metaData$percent.mt) #0.
  hist(metaData$nCount_Spatial) #Up to 120,000.
  hist(metaData$nFeature_Spatial) #Up to 10,000.
  
  # Seurat::SpatialFeaturePlot(data,features="nCount_Spatial")+theme(legend.position="right")
  # Seurat::SpatialFeaturePlot(data,features="nFeature_Spatial")+theme(legend.position="right")
  
  data<-subset(data,
               subset=nCount_Spatial<50000 & nFeature_Spatial>1000 & nFeature_Spatial<8000 & percent.mt<50) #https://www.biostars.org/p/407036/. Warning messages: Not validating Centroids objects, etc.
  dim(data) #19465  3046. 19465 genes, 3046 "cells".
}


#Take the left half.
if(TRUE){
  #https://github.com/satijalab/seurat/issues/3704.
  data_coordinates<-data.frame(cells=data@images$slice1$centroids@cells,
                               spatial_coord_x=data@images$slice1$centroids@coords[,1],
                               spatial_coord_y=data@images$slice1$centroids@coords[,2])
  # Seurat::SpatialDimPlot(data,interactive=TRUE) #x-coodinate plays the role of y-coordinate, and vice versa.
  #Top: x-coordinate: 302. Bottom: x-coordinate: 8934.
  #Left: y-coordinate: 300. Right: y-coordinate: 8686.
  
  #Take the left half.
  range(data_coordinates$spatial_coord_y) #300 8686.
  cutoff<-0.5*(min(data_coordinates$spatial_coord_y)+max(data_coordinates$spatial_coord_y)) #4493.
  data_coordinatesSub<-data_coordinates%>%filter(spatial_coord_y<=cutoff)
  data_leftHalf<-subset(data,cells=data_coordinatesSub$cells)
  
  dim(data_leftHalf) #19465  1420. 19465 genes, 1420 "cells".
  Seurat::SpatialDimPlot(data_leftHalf,pt.size.factor=0.5)
}


saveRDS(data_leftHalf,"~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_spider.rds")


#Required by scCustomize::as.anndata().
reticulate::use_python("~/_Programs/anaconda3/envs/bioinfo/bin")
reticulate::use_condaenv("bioinfo")
reticulate::py_module_available(module="anndata") #TRUE is good.

scCustomize::as.anndata(x=data_leftHalf,
                        file_path="~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/",
                        file_name="dataSpatial_spider.h5ad",
                        main_layer="counts",
                        other_layers=NULL
)
# • Checking Seurat object validity
# • Extracting Data from Spatial assay to transfer to anndata.
# The following columns were removed as they contain identical values for all rows:
#   ℹ orig.ident and percent.mt
# • Creating anndata object.
# • Writing anndata file:
#   "/home/hezhou/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_spider.h5ad"
# AnnData object with n_obs × n_vars = 1420 × 19465
# obs: 'nCount_Spatial', 'nFeature_Spatial'
# var: 'names'


#For future use:
dataSpatial<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSpatial_spider.rds") #78.3 MB.
dim(dataSpatial) #19465  1420. 19465 genes, 1420 "cells".





