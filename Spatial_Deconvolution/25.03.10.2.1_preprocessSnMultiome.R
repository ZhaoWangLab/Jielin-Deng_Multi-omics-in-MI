library(dplyr)


load("/labs/zhaowang/data/Jielin_MI_12052024/R4.3.2_seurat5/SN3_MI14_subset_renameIdent.RData") #This takes a minute or two. 17 GB.
Seurat::DefaultAssay(SN3_MI14_subset_renameIdent) #RNA.
dim(SN3_MI14_subset_renameIdent) #32285 44573. 32285 genes, 44573 cells.


SN3_MI14_subset_renameIdent@assays[["ATAC"]]<-NULL
SN3_MI14_subset_renameIdent@assays[["peaks"]]<-NULL
SN3_MI14_subset_renameIdent@assays[["SCT"]]<-NULL #3 GB.


SN3_MI14_subset_renameIdent@graphs<-vector(mode="list",length=0)
SN3_MI14_subset_renameIdent@neighbors<-vector(mode="list",length=0)
SN3_MI14_subset_renameIdent@reductions<-vector(mode="list",length=0)
SN3_MI14_subset_renameIdent@commands<-vector(mode="list",length=0) #2.6 GB.


metaData<-SN3_MI14_subset_renameIdent@meta.data #44,573*32.
metaData2<-metaData%>%select(nCount_RNA,nFeature_RNA,percent.mt,sample,celltype) #44,573*5.
metaData2$sample<-as.character(metaData2$sample)
metaData2$celltype<-as.character(metaData2$celltype)
SN3_MI14_subset_renameIdent@meta.data<-metaData2


saveRDS(SN3_MI14_subset_renameIdent,"~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSnRNASeq.rds") #This takes a few minutes.



#Required by scCustomize::as.anndata().
reticulate::use_python("~/_Programs/anaconda3/envs/bioinfo/bin")
reticulate::use_condaenv("bioinfo")
reticulate::py_module_available(module="anndata") #TRUE is good.

scCustomize::as.anndata(x=SN3_MI14_subset_renameIdent,
                        file_path="~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/",
                        file_name="dataSnRNASeq.h5ad",
                        main_layer="counts",
                        other_layers=NULL
)
# • Checking Seurat object validity
# • Extracting Data from RNA assay to transfer to anndata.
# • Creating anndata object.
# • Writing anndata file:
#   "/home/hezhou/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSnRNASeq.h5ad"
# AnnData object with n_obs × n_vars = 44573 × 32285
# obs: 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'sample', 'celltype'
# var: 'names'


#For future use:
dataSnRNASeq<-readRDS("~/_Projects/25.02.27_Jielin_Senescence/25.04.07_prepareData/_Data/dataSnRNASeq.rds") #2.6 GB.
dim(dataSnRNASeq) #32285 44573. 32285 genes, 44573 cells.
table(dataSnRNASeq@meta.data[["sample"]]) #MI_1 is spider heart, MI_4 is p16 heart.





