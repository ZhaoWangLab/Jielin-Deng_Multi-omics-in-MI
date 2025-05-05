install.packages('Seurat')
library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)
#install.packages("hdf5r")
library(hdf5r)
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

spider_wholeheart  <- data
spider_wholeheart <- SCTransform(spider_wholeheart, assay = "Spatial", verbose = FALSE)
spider_wholeheart <- RunPCA(spider_wholeheart, assay = "SCT", verbose = FALSE)
spider_wholeheart <- FindNeighbors(spider_wholeheart, reduction = "pca", dims = 1:30)
spider_wholeheart <- FindClusters(spider_wholeheart, verbose = FALSE)
spider_wholeheart <- RunUMAP(spider_wholeheart, reduction = "pca", dims = 1:30, umap.method = "uwot")
p1 <- SpatialDimPlot(spider_wholeheart,  pt.size.factor=1)+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )+
  theme_void() +  # Removes grid and axes
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # Set the central plot area to white
    plot.background = element_rect(fill = "white", color = NA),   # Set outer plot background to white
    panel.grid.major = element_blank(),   # Remove major gridlines
    panel.grid.minor = element_blank())    # Remove minor gridlines
ggsave("spider_wholeheart_spatialplot_bycluster.tiff", plot = p1,limitsize = FALSE, width = 30, height = 28, units = "cm", scale =1)
#Use the  specific cell type biomarker in the snMultiome to perform annotation
dmarkers <- c("Actn2","Tnnt2", "Ryr2","Pecam1", "Tie1","Vegfc","Col1a1","Dcn", "Postn", "Cd68" , "C1qa", "C1qb",  "Msln", "Wt1", "Bnc1", "Myh11", "Tagln", "Acta2","Abcc9", "Rgs5", "Kcnj8",
              "Cd3e", "Cd8a","Cd3g",  "Igkc", "Cr2",  "Cd79a", "Xcr1", "Ccr7", "Flt3")
genelist <- c("Ryr2", "Postn", "C1qa", "Wt1", "Vegfc", "Kcnj8", "Myh11", "Cd3e", "Igkc", "Ccr7")
combined_genelist <- unique(dmarkers, genelist)
# Create output directory if it doesn't exist
output_dir <- "SpatialFeaturePlots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set the default assay
DefaultAssay(spider_wholeheart) <- "SCT"

# Get available features (genes) in the dataset
available_genes <- rownames(spider_wholeheart)

# Loop through each gene in dmarkers
for (gene in dmarkers) {

  # Check if the gene exists in the dataset
  if (!(gene %in% available_genes)) {
    warning(paste("Skipping:", gene, "- Not found in dataset"))
    next  # Skip this iteration
  }

  # Generate SpatialFeaturePlot
  p <- SpatialFeaturePlot(spider_wholeheart, features = gene, image.scale = "hires", image.alpha = 0, pt.size.factor = 0.5) +
    theme(
      text = element_text(family = "Arial"),  # Set all text to Arial
      plot.title = element_text(size = 14, face = "bold"),  # Title size and style
      axis.title = element_text(size = 12),  # Axis title size
      axis.text = element_text(size = 10),  # Axis tick labels size
      legend.text = element_text(size = 10),  # Legend text size
      legend.title = element_text(size = 12)  # Legend title size
    )

  # Define the file name
  filename <- paste0(output_dir, "/", gene, "_SpatialFeaturePlot.tiff")

  # Save the plot
  ggsave(filename, plot = p, width = 10, height = 8, units = "cm", dpi = 300)

  # Print status message
  print(paste("Saved plot for:", gene))
}

print("All SpatialFeaturePlots have been saved!")




#Take the left half.
if(TRUE){
  #https://github.com/satijalab/seurat/issues/3704.
  data_coordinates<-data.frame(cells=spider_wholeheart@images$slice1$centroids@cells,
                               spatial_coord_x=spider_wholeheart@images$slice1$centroids@coords[,1],
                               spatial_coord_y=spider_wholeheart@images$slice1$centroids@coords[,2])
  # Seurat::SpatialDimPlot(spider_wholeheart,interactive=TRUE)
  #Top: x-coordinate: 302. Bottom: x-coordinate: 8934.
  #Left: y-coordinate: 300. Right: y-coordinate: 8686.
  #Take the left half.
  range(data_coordinates$spatial_coord_y) #300 8686.
  cutoff<-0.5*(min(data_coordinates$spatial_coord_y)+max(data_coordinates$spatial_coord_y)) #4493.
  data_coordinatesSub<-data_coordinates%>%filter(spatial_coord_y<=cutoff)
  data_leftHalf<-subset(spider_wholeheart,cells=data_coordinatesSub$cells)
  dim(data_leftHalf) #19465  1420. 19465 genes, 1420 "cells".
  Seurat::SpatialDimPlot(data_leftHalf,pt.size.factor=0.5)
}

spider_halfheart <- data_leftHalf
saveRDS(spider_halfheart,"data_leftHalf_spider.rds")
#SpatialDimPlot(spider_halfheart,pt.size.factor=0.5)


data<-Seurat::Load10X_Spatial(data.dir="/labs/zhaowang/Seq/231006_IGC-JD-21455/gfpsamplenew/outs/",
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

p16_wholeheart  <- data
p16_wholeheart <- SCTransform(p16_wholeheart, assay = "Spatial", verbose = FALSE)
p16_wholeheart <- RunPCA(p16_wholeheart, assay = "SCT", verbose = FALSE)
p16_wholeheart <- FindNeighbors(p16_wholeheart, reduction = "pca", dims = 1:30)
p16_wholeheart <- FindClusters(p16_wholeheart, verbose = FALSE)
p16_wholeheart <- RunUMAP(p16_wholeheart, reduction = "pca", dims = 1:30, umap.method = "uwot")
p1 <- SpatialDimPlot(p16_wholeheart,  pt.size.factor=1)+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("p16_wholeheart_spatialplot_bycluster.tiff", plot = p1,limitsize = FALSE, width = 30, height = 28, units = "cm", scale =1)
#Take the left half.
if(TRUE){
  #https://github.com/satijalab/seurat/issues/3704.
  data_coordinates<-data.frame(cells=p16_wholeheart@images$slice1$centroids@cells,
                               spatial_coord_x=p16_wholeheart@images$slice1$centroids@coords[,1],
                               spatial_coord_y=p16_wholeheart@images$slice1$centroids@coords[,2])
  # Seurat::SpatialDimPlot(p16_wholeheart,interactive=TRUE)
  #Top: x-coordinate: 302. Bottom: x-coordinate: 8934.
  #Left: y-coordinate: 300. Right: y-coordinate: 8686.
  #Take the left half.
  range(data_coordinates$spatial_coord_y) #300 8686.
  cutoff<-0.5*(min(data_coordinates$spatial_coord_y)+max(data_coordinates$spatial_coord_y)) #4493.
  data_coordinatesSub<-data_coordinates%>%filter(spatial_coord_y<=cutoff)
  data_leftHalf<-subset(p16_wholeheart,cells=data_coordinatesSub$cells)
  dim(data_leftHalf) #19465  1420. 19465 genes, 1420 "cells".
  Seurat::SpatialDimPlot(data_leftHalf,pt.size.factor=0.5)
}

p16_halfheart <- data_leftHalf
saveRDS(p16_halfheart,"data_leftHalf_p16.rds")



#process
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
library(Signac)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)
BiocManager::install("biovizBase", force = TRUE)
library(biovizBase)

install.packages('tidyverse')
library(tidyverse)
#library(gridExtra)
BiocManager::install("multtest")
library(multtest)
install.packages('metap')
library(metap)
library(dplyr)

#input RNA and ATAC data from cell ranger aggr
counts <- Read10X_h5("SN3_MI14/outs/filtered_feature_bc_matrix.h5")
fragpath <- "SN3_MI14/outs/atac_fragments.tsv.gz"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA data
SN3_MI14 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  names.field = 2, names.delim = "-"
)

# create ATAC assay and add it to the object
SN3_MI14[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

#calculating nucleosomesignal and TSSEnrichment score
DefaultAssay(SN3_MI14) <- "ATAC"
SN3_MI14 <- NucleosomeSignal(SN3_MI14)
SN3_MI14 <- TSSEnrichment(SN3_MI14)

#create a sample column
SN3_MI14$sample <- rownames(SN3_MI14@meta.data)
view(SN3_MI14@meta.data)

#split sample column
SN3_MI14@meta.data <- separate(SN3_MI14@meta.data, col = "sample", into = c("barcode","samplecode"), sep = "-")
view(SN3_MI14@meta.data)
unique(SN3_MI14@meta.data$samplecode)
#rename group names
Idents(SN3_MI14) <- "samplecode"
table(SN3_MI14$orig.ident)

SN3_MI14 <- RenameIdents(SN3_MI14,"1"="SN3_1", "2"="SN3_2", "3"="MI_1","4"="MI_4")
SN3_MI14$sample <- Idents(SN3_MI14)
#SN3_MI14@meta.data <- separate(SN3_MI14@meta.data, col = "sample", into = c("group","sample_number"), sep = "_")
view(SN3_MI14@meta.data)

#obtain mt data
Idents(SN3_MI14) <- "sample"
DefaultAssay(SN3_MI14) <- "RNA"
SN3_MI14[["percent.mt"]] <- PercentageFeatureSet(SN3_MI14, pattern = "^mt-")

DefaultAssay(SN3_MI14) <- "ATAC"
#obtain peaks and SCT assays in all data SN3_MI14_subset
peaks <- CallPeaks(SN3_MI14, macs2.path="/home/diren/anaconda3/envs/macs3_env/bin/macs3")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
macs2_counts <- FeatureMatrix(
  fragments = Fragments(SN3_MI14),
  features = peaks,
  cells = colnames(SN3_MI14)
)
SN3_MI14[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
#quality of library
DefaultAssay(SN3_MI14) <- "ATAC"

p1 <- VlnPlot(
  object = SN3_MI14,
  features = c("nCount_RNA", "nCount_ATAC", "percent.mt", "TSS.enrichment", "nucleosome_signal"),
  ncol = 5,
  pt.size = 0,
)

ggsave("SN3_MI14_QC_beforefilter.tiff", plot = p1,limitsize = FALSE, width = 100, height = 28, units = "cm", scale = 0.5)
save(SN3_MI14, file = "SN3_MI14.RData")

#subset to remove low quality cells
SN3_MI14_subset <- subset(
  x = SN3_MI14,
  subset = nCount_ATAC < 150000 &
    nCount_RNA < 100000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 5
)
p1 <- VlnPlot(
  object = SN3_MI14_subset,
  features = c("nCount_RNA", "nCount_ATAC", "percent.mt", "TSS.enrichment", "nucleosome_signal"),
  ncol = 5,
  pt.size = 0,
)
ggsave("SN3_MI14_QC_afterfilter.tiff", plot = p1,limitsize = FALSE, width = 100, height = 28, units = "cm", scale = 0.5)

BiocManager::install('glmGamPoi')
library(glmGamPoi)

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(patchwork)
DefaultAssay(SN3_MI14_subset) <- "RNA"
options(future.globals.maxSize = 8000 * 1024^3)
future::plan("sequential")
SN3_MI14_subset <- NormalizeData(SN3_MI14_subset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(dims = 1:20) %>% FindNeighbors() %>% FindClusters() %>% RunUMAP(dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#pK Identification (no ground-truth)
sweep.res.list<- paramSweep(SN3_MI14_subset, PCs= 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate

annotations <- SN3_MI14_subset@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*nrow(SN3_MI14_subset@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#run doubletFinder
SN3_MI14_subset <- doubletFinder(SN3_MI14_subset, PCs = 1:20, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = TRUE)
view(SN3_MI14_subset@meta.data)
#visualize doublets
p1 <- DimPlot(SN3_MI14_subset,  group.by = "DF.classifications_0.25_0.3_3666", repel = TRUE, raster = FALSE) + ggtitle("Doublet proportion")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_doublet.tiff", plot = p1,limitsize = FALSE, width = 30, height = 28, units = "cm", scale = 0.5)

#number of singlets and doublets
table(SN3_MI14_subset@meta.data$DF.classifications_0.25_0.3_3666)
SN3_MI14_subset <- subset(x = SN3_MI14_subset, subset= DF.classifications_0.25_0.3_3666 == "Singlet" )
save(SN3_MI14_subset, file="SN3_MI14_subset.RData")
#Run normalization with RNA data
DefaultAssay(SN3_MI14_subset) <- "RNA"
options(future.globals.maxSize = 8000 * 1024^3)
future::plan("sequential")
SN3_MI14_subset <- SCTransform(SN3_MI14_subset, verbose = TRUE)
SN3_MI14_subset <- RunPCA(SN3_MI14_subset)
DefaultAssay(SN3_MI14_subset) <- "peaks"
SN3_MI14_subset <- FindTopFeatures(SN3_MI14_subset, min.cutoff = 30)
SN3_MI14_subset <- RunTFIDF(SN3_MI14_subset)
SN3_MI14_subset <- RunSVD(SN3_MI14_subset)

#strategy one to integrate RNA and ATAC assays
SN3_MI14_subset <- FindMultiModalNeighbors(
  object = SN3_MI14_subset,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
SN3_MI14_subset <- RunUMAP(
  object = SN3_MI14_subset,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)
#assign the chromatin assay with a new name so the second stratey of integration will overwrite the peaks assay
SN3_MI14_subset[["peaks_strategy1"]] <- SN3_MI14_subset[["peaks"]]
#plot the integrated umap from strategy one
p1 <- DimPlot(SN3_MI14_subset, reduction = "umap",  label = TRUE,label.size = 4 ,repel = TRUE, raster = FALSE) + ggtitle("UMAP")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_umap_label.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)



#Strategy two: WNN analysis to integrate RNA and ATAC data together, the benefit is we still can umap ATAC and RNA separately to study them if we want
DefaultAssay(SN3_MI14_subset) <- "RNA"
options(future.globals.maxSize = 8000 * 1024^3)
future::plan("sequential")
SN3_MI14_subset <- SCTransform(SN3_MI14_subset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(SN3_MI14_subset) <- "ATAC"
SN3_MI14_subset <- RunTFIDF(SN3_MI14_subset)
SN3_MI14_subset <- FindTopFeatures(SN3_MI14_subset, min.cutoff = 'q0')
SN3_MI14_subset <- RunSVD(SN3_MI14_subset)
SN3_MI14_subset <- RunUMAP(SN3_MI14_subset, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
SN3_MI14_subset <- FindMultiModalNeighbors(SN3_MI14_subset, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
SN3_MI14_subset <- RunUMAP(SN3_MI14_subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
SN3_MI14_subset <- FindClusters(SN3_MI14_subset, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

save(SN3_MI14_subset, file = "SN3_MI14_subset.RData")

p1 <- DimPlot(SN3_MI14_subset, reduction = "umap.rna", label = FALSE, label.size = 2.5, repel = TRUE, raster = FALSE ) + ggtitle("RNA")+
theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("SN3_MI14_subset_wnn_rna.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

p1 <- DimPlot(SN3_MI14_subset, reduction = "umap.atac", label = FALSE, label.size = 2.5, repel = TRUE, raster = FALSE) + ggtitle("ATAC")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_wnn_atac.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

p1 <- DimPlot(SN3_MI14_subset, reduction = "wnn.umap",  label = TRUE, label.size = 4, repel = TRUE, raster = FALSE) + ggtitle("WNN")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_wnn_wnn.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)
p1 <- DimPlot(SN3_MI14_subset, reduction = "wnn.umap",  label = FALSE, repel = TRUE, raster = FALSE) + ggtitle("WNN")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_wnn_wnn_nolabel.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

#celltype annotation with specific markers
Idents(SN3_MI14_subset) <- SN3_MI14_subset$seurat_clusters
dmarkers <- c("Actn2","Tnnt2", "Ryr2","Pecam1", "Tie1","Vegfc","Col1a1","Dcn", "Postn", "Cd68" , "C1qa", "C1qb",  "Msln", "Wt1", "Bnc1", "Myh11", "Tagln", "Acta2","Abcc9", "Rgs5", "Kcnj8",
              "Cd3e", "Cd8a","Cd3g",  "Igkc", "Cr2",  "Cd79a", "Xcr1", "Ccr7", "Flt3")
DefaultAssay(SN3_MI14_subset) <- "SCT"

p0 <- DotPlot(SN3_MI14_subset, features = dmarkers, group.by = "seurat_clusters", cols = c("blue", "red"))+RotatedAxis()+coord_flip()+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_celltype_marker_dotplot_bycluster1.tiff", plot = p0, width = 52, height = 20, units = "cm", scale = 1)


SN3_MI14_subset_renameIdent <- SN3_MI14_subset
SN3_MI14_subset_renameIdent <- RenameIdents(SN3_MI14_subset_renameIdent, "6"= "Cardiomyocytes","14"= "Cardiomyocytes",
                                            "3"= "Endothelial cells" ,"8"= "Endothelial cells", "10"= "Endothelial cells", "19"= "Endothelial cells","23"= "Endothelial cells","29"= "Endothelial cells",
                                            "0"= "Fibroblasts", "5"= "Fibroblasts" ,"7"= "Fibroblasts",  "9"= "Fibroblasts","11"= "Fibroblasts","12"= "Fibroblasts" ,"15"= "Fibroblasts",  "18"= "Fibroblasts","26"= "Fibroblasts",
                                            "1"= "Macrophages","2"= "Macrophages", "13"= "Macrophages" , "17"= "Macrophages" , "20"= "Macrophages" , "22"= "Macrophages" , "24"= "Macrophages" , "27"= "Macrophages" ,
                                            "16"= "Epicardial cells","28"= "Smooth muscle cells" ,
                                            "4"= "Pericytes","21"= "Pericytes",
                                            "25"= "Immune-T cells","30"= "Immune-B cells",
                                            "31"= "Dentritic cells", "32"= "Dentritic cells")

SN3_MI14_subset_renameIdent$celltype <- Idents(SN3_MI14_subset_renameIdent)
p1 <- DimPlot(SN3_MI14_subset_renameIdent,  reduction = "umap",group.by="celltype",label = FALSE, repel = TRUE, raster = FALSE )+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_annotation_celltype.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)
Idents(SN3_MI14_subset_renameIdent) <- "celltype"
DefaultAssay(SN3_MI14_subset_renameIdent) <- "RNA"
p1 <- DotPlot(SN3_MI14_subset_renameIdent, features = dmarkers, group.by = "celltype",cols = c("blue", "red"))+RotatedAxis()+coord_flip()+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_renameIdent_celltype_marker_dotplot_bycelltype.tiff", plot = p1, width = 30, height = 20, units = "cm", scale = 1)
save(SN3_MI14_subset_renameIdent, file = "SN3_MI14_subset_renameIdent.RData")
DimPlot(SN3_MI14_subset_renameIdent,  reduction = "wnn.umap",group.by="celltype",label = FALSE, repel = TRUE, raster = FALSE )
#Get p16 expression score to identify senescence cells
install.packages("wesanderson")
library(wesanderson)
pal <- wes_palette("Zissou1",100, type = "continuous")

DefaultAssay(SN3_MI14_subset_renameIdent) <- "SCT"
p1 <- FeaturePlot(SN3_MI14_subset_renameIdent, feature = "Cdkn2a", reduction = "umap", pt.size = 0.5,  max.cutoff = "q80")+ scale_color_gradientn(colors = pal)+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_renameIdent_Featurep16.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)
#vlnplot of p16 expression in celltypes
p1 <- VlnPlot(SN3_MI14_subset_renameIdent, features = "Cdkn2a", group.by="celltype",pt.size = 0, log = TRUE)+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_renameIdent_vlnplotp16.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

FeaturePlot(SN3_MI14_subset_renameIdent, feature = "Cdkn2a", reduction = "umap", pt.size = 0.5,  max.cutoff = "q80")+ scale_color_gradientn(colors = pal)

#get senesence enrichment score by p16 expression
library(scales)
senescence_marker= "Cdkn2a"
SN3_MI14_subset_renameIdent=AddModuleScore(SN3_MI14_subset_renameIdent, features = list(senescence_marker), assay = "RNA")
SN3_MI14_subset_renameIdent$Cluster1-> SN3_MI14_subset_renameIdent$Senescence_signature
colors <- rainbow(16, start=0.1, end=0.9)
senenorm<-scales::rescale(x=SN3_MI14_subset_renameIdent$Senescence_signature, c(0,1))
SN3_MI14_subset_renameIdent$Senescence_score<-senenorm
##########  Identification of senescent cells
# add a column to separate a scale in metadata
poscell <- subset(SN3_MI14_subset_renameIdent, Senescence_score>=0.5)
cellnames <- rownames(poscell@meta.data)
SN3_MI14_subset_renameIdent$barcode <-rownames(SN3_MI14_subset_renameIdent@meta.data)
SN3_MI14_subset_renameIdent@meta.data <- SN3_MI14_subset_renameIdent@meta.data %>% mutate(senescence_status = ifelse((SN3_MI14_subset_renameIdent$barcode %in% cellnames), "p16-Positive",  "p16-Negative"))
p1 <-DimPlot(SN3_MI14_subset_renameIdent,  reduction = "umap", pt.size = 0.5, label = FALSE, repel = TRUE, raster = FALSE, group.by = "senescence_status",cols= c( "p16-Positive"="red","p16-Negative"="gray" ))+ ggtitle("p16+ cells")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_wnn_senescence_p16+cells_identification.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

#senescence cell proportion

Idents(SN3_MI14_subset_renameIdent) <- "celltype"
pt <- table(Idents(SN3_MI14_subset_renameIdent), SN3_MI14_subset_renameIdent$senescence_status)
write.csv(pt, file = "SN3_MI14_subset_renameIdent_p16_senescence_proportion.csv")
#pt <- read.csv(file = "SN3_MI14_subset_renameIdent_senescence_proportion.csv", header=TRUE)
pt <- as.data.frame(pt)
row_sub = apply(pt, 1, function(row) all(row !=0 ))
pt[row_sub,]
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1, levels=unique(pt$Var1))

level_order <- c("Epicardial cells",
                 "Fibroblasts",
                 "Pericytes",
                 "Dentritic cells",
                 "Endothelial cells",
                 "Macrophages",
                 "Immune-B cells",
                 "Smooth muscle cells",
                 "Immune-T cells",
                 "Cardiomyocytes"
)
p1 <- ggplot(pt, aes(x=factor(Var1, levels = level_order), y=Freq, fill=Var2)) +
  theme_bw(base_size = 15) + theme(aspect.ratio = 2/1) +
  geom_col(position = "fill", width = 0.5) +
  xlab("cell type") +
  ylab("p16+ cell proportion")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )

ggsave("SN3_MI14_subset_renameIdent_p16_senescence_proportion.tiff", plot = p1, width = 18, height = 8, units = "cm", scale = 3)


#visulize the distribution of p16+ cells
#library(ggplot2)
# Extract metadata including cell type and score
df <- data.frame(
  CellType = SN3_MI14_subset_renameIdent$celltype,  # Replace with your metadata column name
  Score = SN3_MI14_subset_renameIdent$Senescence_score         # Replace with the actual score column
)

# Create scatter plot
p1 <- ggplot(df, aes(x = CellType, y = Score, color = Score >= 0.5)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1) +  # Add jitter to spread points
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "green")) +  # Define colors
  labs(x = "Cell type", y = "p16 Score", title = "P16 expression score in Fibroblasts") +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +  # Add horizontal separator
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12,  color = "black"),  # Make axis labels bold and bigger
    axis.title = element_text(size = 14, color = "black"),  # Make axis titles bold and bigger
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title and increase size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    text = element_text(family = "Arial")  # Set all text to Arial
  )
ggsave(
  filename = "SN3_MI14_subset_renameIdent_p16+score distribution.tiff",
  plot = p1,
  width = 12,
  height = 6,
  dpi = 900
)






#We identified that Fibroblasts to have most p16+ cells using snMultiome dataset only, now try to further analyze Fb
view(SN3_MI14_subset_renameIdent@meta.data)
Idents(SN3_MI14_subset_renameIdent) <- "celltype"
SN3_MI14_Fb_subset <- subset(SN3_MI14_subset_renameIdent, idents = "Fibroblasts")
DefaultAssay(SN3_MI14_Fb_subset) <- "SCT"
# Run PCA on the subsetted data
SN3_MI14_Fb_subset <- RunPCA(SN3_MI14_Fb_subset)
# Run UMAP on the subsetted data
SN3_MI14_Fb_subset <- RunUMAP(SN3_MI14_Fb_subset, dims = 1:20)
# Re-cluster the subset
SN3_MI14_Fb_subset <- FindNeighbors(SN3_MI14_Fb_subset, dims = 1:20)
SN3_MI14_Fb_subset <- FindClusters(SN3_MI14_Fb_subset, resolution = 0.4)
view(SN3_MI14_Fb_subset@meta.data)

Idents(SN3_MI14_subset_renameIdent) <- "SCT_snn_res.0.4"

SN3_MI14_Fb_subset$subcluster <- SN3_MI14_Fb_subset$SCT_snn_res.0.4
p0 <- DimPlot(SN3_MI14_Fb_subset, reduction = "umap", pt.size = 0.5, group.by="subcluster", label = TRUE, repel = TRUE, raster = FALSE)+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("SN3_MI14_Fb_subset_Dimplot.tiff", plot = p0, width = 24, height = 20, units = "cm", scale = 1)

dmarkers <- c("Hmcn2","Col6a6","Mt2", "Ror1","Lars2", "Camk1d", "Top2a", "Neil3",   "Cacna1c", "Postn","Itgbl1", "Nox4","Col8a1",  "Prkca")

p0 <- DotPlot(SN3_MI14_Fb_subset, features = dmarkers, group.by = "subcluster", cols = c("blue", "red"))+RotatedAxis()+coord_flip()+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
p0
ggsave("SN3_MI14_Fb_subset_dotplot_bysubcluster.tiff", plot = p0, width = 20, height = 15, units = "cm", scale = 1)


SN3_MI14_Fb_subset_renameIdent <- SN3_MI14_Fb_subset
SN3_MI14_Fb_subset_renameIdent <- RenameIdents(SN3_MI14_Fb_subset_renameIdent,  "5"="Uninjured QuiCF","7"="Acute QuiCF",  "8"="Early ActCF","4"="Early MF","0"="Late MF","1"="Late MF","3"="Late MF", "6"="MaF","9"="MaF", "2"="Chronic QuiCF")
SN3_MI14_Fb_subset_renameIdent$subcelltype <- Idents(SN3_MI14_Fb_subset_renameIdent)
p0 <- DimPlot(SN3_MI14_Fb_subset_renameIdent, reduction = "umap",group.by = "subcelltype" ,pt.size = 0.5, label = FALSE, repel = TRUE, raster = FALSE) +theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("SN3_MI14_Fb_subset_renameIdent_Dimplot_bysubcelltype.tiff", plot = p0, width = 24, height = 20, units = "cm", scale = 1)

#annotation marker
p0 <- DotPlot(SN3_MI14_Fb_subset_renameIdent, features = dmarkers, group.by = "subcelltype", cols = c("blue", "red"))+RotatedAxis()+coord_flip()+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
p0
ggsave("SN3_MI14_Fb_subset_renameIdent_dotplot_bysubcelltype.tiff", plot = p0, width = 20, height = 15, units = "cm", scale = 1)


#visulize the p16+cells distribution in Fb
p0 <- DimPlot(SN3_MI14_Fb_subset_renameIdent, reduction = "umap", pt.size = 0.5, group.by = "senescence_status",cols= c( "p16-Positive"="red","p16-Negative"="gray" ))+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("SN3_MI14_Fb_subset_renameIdent_senescence__Dimplot.tiff", plot = p0, width = 24, height = 20, units = "cm", scale = 1)
save(SN3_MI14_Fb_subset_renameIdent, file = "SN3_MI14_Fb_subset_renameIdent.RData")
#visulize the distribution of p16+ cells in Fb subcelltypes
library(ggplot2)

# Extract metadata including cell type and score
df <- data.frame(
  CellType = SN3_MI14_Fb_subset_renameIdent$subcelltype,  # Replace with your metadata column name
  Score = SN3_MI14_Fb_subset_renameIdent$Senescence_score         # Replace with the actual score column
)

# Create scatter plot
p1 <- ggplot(df, aes(x = CellType, y = Score, color = Score >= 0.5)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1) +  # Add jitter to spread points
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "green")) +  # Define colors
  labs(x = "Fibroblasts subtype", y = "p16 Score", title = "P16 expression score in Fibroblasts") +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +  # Add horizontal separator
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12,  color = "black"),  # Make axis labels bold and bigger
    axis.title = element_text(size = 14, color = "black"),  # Make axis titles bold and bigger
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title and increase size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    text = element_text(family = "Arial")  # Set all text to Arial
    )
ggsave(
  filename = "p16+score distribution.tiff",
  plot = p1,
  width = 12,
  height = 6,
  dpi = 900
)

#senescence cell proportion
pt <- table(Idents(SN3_MI14_Fb_subset_renameIdent), SN3_MI14_Fb_subset_renameIdent$senescence_status)
write.csv(pt, file = "SN3_MI14_Fb_subset_renameIdent_senescence_proportion.csv")



#Now we identified that late MF, early MF are the Fb subcelltype that have most p16+ marked senesence cells, we next investigate the mechanism
#Get DEG and pathway enrichment for p16+ vs p16- (lateMF, early MF, total MF)
Idents(SN3_MI14_Fb_subset_renameIdent) <- "senescence_status"
current_levels <- levels(Idents(SN3_MI14_Fb_subset_renameIdent))
# Reverse the order of the levels
reversed_levels <- rev(current_levels)
# Set the new levels in the Seurat object
Idents(SN3_MI14_Fb_subset_renameIdent) <- factor(Idents(SN3_MI14_Fb_subset_renameIdent), levels = reversed_levels)
# Check the new order of levels
levels(Idents(SN3_MI14_Fb_subset_renameIdent))

# Define the subcell types
subcell_types <- c("Late MF", "Early MF", "Total MF")

# Create output directory if it doesn't exist
output_dir <- "DEG_Enrichment_MF"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each subcell type and run FindMarkers
for (subcell in subcell_types) {

  # Create a unique variable name for the subset
  subset_name <- paste0(gsub(" ", "", subcell), "_Fb_subset")
  Idents(SN3_MI14_Fb_subset_renameIdent) <- "subcelltype"

  # Subset the data based on subcell type
  if (subcell == "Total MF") {
    subcell_subset <- subset(SN3_MI14_Fb_subset_renameIdent, idents = c("Late MF", "Early MF"))
  } else {
    subcell_subset <- subset(SN3_MI14_Fb_subset_renameIdent, idents = subcell)
  }

  # Assign the subset to a unique variable name in the environment
  assign(subset_name, subcell_subset)

  # Set RNA as the default assay
  Idents(subcell_subset) <- "senescence_status"
  DefaultAssay(subcell_subset) <- "RNA"

  # Define output filename
  output_filename <- paste0(output_dir, "/", gsub(" ", "", subcell), "_p16pvsp16n_markers.csv")

  # Run FindMarkers
  markers <- FindMarkers(
    subcell_subset,
    ident.1 = "p16-Positive",
    ident.2 = "p16-Negative"
  )

  # Convert row names (gene names) into a proper column
  markers <- cbind(gene = rownames(markers), markers)

  # Save results to CSV with gene names as the first column
  write.csv(markers, file = output_filename, row.names = FALSE)

  # Print status message
  print(paste("Finished processing:", subcell, "Results saved to:", output_filename))

  # Print message indicating subset object was created
  print(paste("Subset saved as:", subset_name))
}



#Now performing GO enrichment of the DEG from above FindMarker results
library(clusterProfiler)
library(dplyr)
library(ggplot2)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# Create output directory if it doesn't exist
output_dir <- "DEG_Enrichment_MF"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to perform GO enrichment analysis
getGOEnrichmentResult <- function(filename, FCFilter) {

  cat("Processing file:", filename, ", FCFilter=", FCFilter, "...\n", sep="")

  # Read the DEG file
  data_FDRFiltered <- read.csv(filename)

  # Apply log2FC filtering
  if (FCFilter == "FCAbove1.5") {
    data_FCFiltered <- data_FDRFiltered %>% filter(p_val < 0.05 & avg_log2FC >= 0.5)
  } else if (FCFilter == "FCUnderNeg1.5") {
    data_FCFiltered <- data_FDRFiltered %>% filter(p_val < 0.05 & avg_log2FC <= -0.5)
  } else if (FCFilter == "FCBothDirections") {
    data_FCFiltered <- data_FDRFiltered %>% filter(p_val < 0.05 & abs(avg_log2FC) >= 0.5)
  }

  # Perform GO enrichment
  resultEnrichGO <- clusterProfiler::enrichGO(
    gene = data_FCFiltered$gene,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  result <- resultEnrichGO@result
  result2 <- result %>% filter(p.adjust <= 0.05)

  # Generate output filenames based on input file
  file_prefix <- gsub("_p16pvsp16n_markers.csv", "", basename(filename))

  output_csv <- paste0(output_dir, "/", file_prefix, "_EnrichedGOList_", FCFilter, ".csv")
  write.csv(result2, file = output_csv, row.names = FALSE)

  # Save dotplot of the top 20 enriched GO terms
  output_pdf <- paste0(output_dir, "/", file_prefix, "_Top20EnrichedGO_", FCFilter, ".pdf")

  clusterProfiler::dotplot(resultEnrichGO, orderBy = "p.adjust", x = "p.adjust", showCategory = 20, font.size = 9)
  ggsave(output_pdf, width = 25, height = 30, unit = "cm")

  cat("Finished processing:", filename, "\nResults saved to:", output_csv, "and", output_pdf, "\n\n")
}

# Get all marker files from the output directory
marker_files <- list.files(output_dir, pattern = "_p16pvsp16n_markers.csv", full.names = TRUE)

# Loop through each marker file and apply all FCFilter options
for (file in marker_files) {
  getGOEnrichmentResult(file, FCFilter = "FCAbove1.5")
  getGOEnrichmentResult(file, FCFilter = "FCUnderNeg1.5")
  getGOEnrichmentResult(file, FCFilter = "FCBothDirections")
}


#plot customized dotplot from clusterProfiler


#Strategy 1: projecting spatial spider+ and spider- from half heart to snMultiome dataset
#spider.halfheart <- readRDS("data_leftHalf_spider.rds")
#p16.halfheart <- readRDS("data_leftHalf_p16.rds")

#processing spatial spider
spiderposition <- read.csv("spiderp_spidern.csv", header = TRUE)

spider_halfheart@meta.data$spiderstatus <- ifelse(rownames(spider_halfheart@meta.data) %in% spiderposition$Barcode, "SPiDER+", "SPiDER-")
view(spider_halfheart@meta.data)
# Set the DefaultAssay to this new assay
DefaultAssay(SN3_MI14_subset_renameIdent) <- "SCT"
# Identify variable features for both reference and query
spider_halfheart <- FindVariableFeatures(spider_halfheart, selection.method = "vst", nfeatures = 2000)
SN3_MI14_subset_renameIdent <- FindVariableFeatures(SN3_MI14_subset_renameIdent, selection.method = "vst", nfeatures = 2000)
# Get common features
shared_features <- intersect(VariableFeatures(spider_halfheart), VariableFeatures(SN3_MI14_subset_renameIdent))

# Check if there are enough shared features
length(shared_features)  # Should return a number > 0

anchors <- FindTransferAnchors(reference = spider_halfheart, query = SN3_MI14_subset_renameIdent, dims = 1:30,
                               reference.reduction = "pca", query.assay = "SCT", k.filter = NA,  features = shared_features )
predictions_spider <- TransferData(anchorset = anchors, refdata = spider_halfheart$spiderstatus, dims = 1:30)
SN3_MI14_subset_renameIdent <- AddMetaData(SN3_MI14_subset_renameIdent, metadata = predictions_spider)
colnames(SN3_MI14_subset_renameIdent@meta.data)[colnames(SN3_MI14_subset_renameIdent@meta.data) == "predicted.id"] <- "predicted.spider.status"
SN3_MI14_subset_renameIdent@meta.data$predicted.spider.status <- factor(
  SN3_MI14_subset_renameIdent@meta.data$predicted.spider.status,
  levels = c("SPiDER+", "SPiDER-") # Define order
)
#for the predicted spider.status, it almost show a lot of SPiDER+ cells now try to use spider score >=0.5 as indicator to identify SPiDER+ cells
# add a column to separate a scale in metadata，

poscell <- subset(SN3_MI14_subset_renameIdent, prediction.score.SPiDER..1>=0.5)
cellnames <- rownames(poscell@meta.data)
SN3_MI14_subset_renameIdent$barcode <-rownames(SN3_MI14_subset_renameIdent@meta.data)
SN3_MI14_subset_renameIdent@meta.data <- SN3_MI14_subset_renameIdent@meta.data %>% mutate(spider_status = ifelse((SN3_MI14_subset_renameIdent$barcode %in% cellnames), "SPiDER+",  "SPiDER-"))
p1 <-DimPlot(SN3_MI14_subset_renameIdent,  reduction = "umap", pt.size = 0.5, label = FALSE, repel = TRUE, raster = FALSE, group.by = "spider_status",cols= c( "SPiDER+"="red","SPiDER-"="gray" ))+ ggtitle("SPiDER+ cells")+
  theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )
ggsave("SN3_MI14_subset_SPiDER+cells_identification.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)

p1 <- DimPlot(SN3_MI14_subset_renameIdent,  reduction = "umap",group.by="predicted.spider.status"  , label = FALSE, repel = TRUE, raster = FALSE, pt.size = 0.4)
ggsave("SN3_MI14_subset_wnn_predicted_spider_usehalfheart.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)
p1 <- DimPlot(SN3_MI14_subset_renameIdent, reduction = "umap", pt.size = 0.5, group.by = "predicted.spider.status", label = FALSE, repel = TRUE, raster = FALSE)+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("SN3_MI14_subset_renameIdent_SPIDER_projection.tiff", plot = p1, width = 24, height = 20, units = "cm", scale = 1)
view(SN3_MI14_subset_renameIdent@meta.data)
# Extract metadata including cell type and score
df <- data.frame(
  CellType = SN3_MI14_subset_renameIdent$celltype,  # Replace with your metadata column name
  Score = SN3_MI14_subset_renameIdent$prediction.score.SPiDER..1         # Replace with the actual score column
)
write.csv(df, file = "spidersocre.csv")
# Create scatter plot
p1 <- ggplot(df, aes(x = CellType, y = Score, color = Score >= 0.5)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1) +  # Add jitter to spread points
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "green")) +  # Define colors
  labs(x = "Celltype", y = "SPiDER+ Score", title = "SPiDER+ score projection in snMultiome") +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +  # Add horizontal separator
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12,  color = "black"),  # Make axis labels bold and bigger
    axis.title = element_text(size = 14, color = "black"),  # Make axis titles bold and bigger
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title and increase size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    text = element_text(family = "Arial")  # Set all text to Arial
  )
ggsave(
  filename = "SPiDER+score distribution.tiff",
  plot = p1,
  width = 12,
  height = 6,
  dpi = 900
)
SN3_MI14_subset_renameIdent@meta.data$spider_status <- factor(
  SN3_MI14_subset_renameIdent@meta.data$spider_status,
  levels = c("SPiDER+", "SPiDER-") # Define order
)
Idents(SN3_MI14_subset_renameIdent) <- "celltype"
pt <- table(Idents(SN3_MI14_subset_renameIdent), SN3_MI14_subset_renameIdent$spider_status)
write.csv(pt, file = "SN3_MI14_subset_renameIdent_spider_status_proportion.csv")
#pt <- read.csv(file = "SN3_MI14_subset_renameIdent_senescence_proportion.csv", header=TRUE)
pt <- as.data.frame(pt)
row_sub = apply(pt, 1, function(row) all(row !=0 ))
pt[row_sub,]
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1, levels=unique(pt$Var1))

level_order <- c(
                 "Fibroblasts",
                 "Epicardial cells","Macrophages","Smooth muscle cells","Dentritic cells","Immune-T cells","Immune-B cells",
                 "Pericytes",
                 "Endothelial cells",
                 "Cardiomyocytes"
)
p1 <- ggplot(pt, aes(x=factor(Var1, levels = level_order), y=Freq, fill=Var2)) +
  theme_bw(base_size = 15) + theme(aspect.ratio = 2/1) +
  geom_col(position = "fill", width = 0.5) +
  xlab("cell type") +
  ylab("SPiDER+ cell proportion")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(
    text = element_text(family = "Arial"),  # Set all text to Arial
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),  # Axis tick labels size
    legend.text = element_text(size = 12),  # Legend text size
    legend.title = element_text(size = 14)  # Legend title size
  )

ggsave("SN3_MI14_subset_renameIdent_spider_status_proportion.tiff", plot = p1, width = 18, height = 8, units = "cm", scale = 3)
save(SN3_MI14_subset_renameIdent, file = "SN3_MI14_subset_renameIdent.RData")
save(SN3_MI14_Fb_subset_renameIdent, file = "SN3_MI14_Fb_subset_renameIdent.RData")

#Stategy2: project snMultiome cell type annotation to spatial transcriptome to identify the spider+ cell cell type composition
#Software1: RCTD
library(spacexr)
#processing the snMultiome reference
counts <- SN3_MI14_subset_renameIdent@assays$SCT$counts
meta_data <- SN3_MI14_subset_renameIdent@meta.data
cell_types <- meta_data$celltype; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nCount_SCT; names(nUMI) <- meta_data$barcode # create nUMI named list

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
#> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
#> colSums of counts. If this is unintended, please correct this discrepancy. If
#> this is intended, there is no problem.

## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type
save(reference, file="reference.RData")
load("reference.RData")
#processing the spatial query, now spider
counts_spatial <- spider_halfheart@assays$Spatial$counts
spider_cells <- colnames(spider_halfheart@assays$Spatial)
coords <- GetTissueCoordinates(spider_halfheart)[spider_cells, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_spatial, colSums(counts_spatial))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
spider_halfheart <- AddMetaData(spider_halfheart, metadata = RCTD@results$results_df)
save(spider_halfheart, file = "spider_halfheart.RData")
#visulize RCTD results
Idents(spider_halfheart) <- "first_type"
# now we can spatially map the location of any scRNA-seq cell type
cells <- CellsByIdentities(spider_halfheart)
p <- SpatialDimPlot(spider_halfheart, cells.highlight = cells, cols.highlight = c("#FFFF00", "grey50"),pt.size.factor = 1, facet.highlight = T, combine = T, ncol = 4)+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("spider_halfheart_RCTD.tiff", plot = p, width = 30, height = 30, units = "cm", scale = 1)


#processing the spatial query, now p16
counts_spatial <- p16_halfheart@assays$Spatial$counts
p16_cells <- colnames(p16_halfheart@assays$Spatial)
coords <- GetTissueCoordinates(p16_halfheart)[p16_cells, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_spatial, colSums(counts_spatial))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
p16_halfheart <- AddMetaData(p16_halfheart, metadata = RCTD@results$results_df)
save(p16_halfheart, file = "p16_halfheart.RData")

#visulize RCTD results
Idents(p16_halfheart) <- "first_type"
# now we can spatially map the location of any scRNA-seq cell type
cells <- CellsByIdentities(p16_halfheart)
p <- SpatialDimPlot(p16_halfheart, cells.highlight = cells, cols.highlight = c("#FFFF00", "grey50"),pt.size.factor = 1, facet.highlight = T, combine = T, ncol = 4)+theme(
  text = element_text(family = "Arial"),  # Set all text to Arial
  plot.title = element_text(size = 16, face = "bold"),  # Title size and style
  axis.title = element_text(size = 14),  # Axis title size
  axis.text = element_text(size = 12),  # Axis tick labels size
  legend.text = element_text(size = 12),  # Legend text size
  legend.title = element_text(size = 14)  # Legend title size
)
ggsave("p16_halfheart_RCTD.tiff", plot = p, width = 30, height = 30, units = "cm", scale = 1)



#chromatin assay analysis to invetgiate TF enrichment for the senescence LateMF population


