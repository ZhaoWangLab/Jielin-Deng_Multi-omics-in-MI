library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)
library(hdf5r)

#input 10xvisium spatial dataset
data<-Seurat::Load10X_Spatial(data.dir="path/outs/",
                              filename="filtered_feature_bc_matrix.h5") 


#Filter cells. 
if(TRUE){
  data[["percent.mt"]]<-Seurat::PercentageFeatureSet(data,
                                                     pattern="^mt-") 
  metaData<-data@meta.data 
  data<-subset(data,
               subset=nCount_Spatial<50000 & nFeature_Spatial>1000 & nFeature_Spatial<8000 & percent.mt<50) 
  dim(data) 
}

wholeheart  <- data
wholeheart <- SCTransform(wholeheart, assay = "Spatial", verbose = FALSE)
wholeheart <- RunPCA(wholeheart, assay = "SCT", verbose = FALSE)
wholeheart <- FindNeighbors(wholeheart, reduction = "pca", dims = 1:30)
wholeheart <- FindClusters(wholeheart, verbose = FALSE)
wholeheart <- RunUMAP(wholeheart, reduction = "pca", dims = 1:30, umap.method = "uwot")


#Take the infarct region of heart
DefaultAssay(spidersample_subset) <-"Spatial"
view(wholeheart@images$slice1)
SpatialDimPlot(wholeheart, cells.highlight = WhichCells(wholeheart, expression = slice1_imagerow > 600 | slice1_imagecol < 300))
spidersample_infarction <- subset(wholeheart, slice1_imagerow > 600 | slice1_imagecol < 300, invert = FALSE)


#process
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(biovizBase)
library(tidyverse)
library(gridExtra)
library(multtest)
library(metap)
library(dplyr)

#input RNA and ATAC data from cell ranger aggr
counts <- Read10X_h5("path/outs/filtered_feature_bc_matrix.h5")
fragpath <- "path/outs/atac_fragments.tsv.gz"

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

#split sample column
SN3_MI14@meta.data <- separate(SN3_MI14@meta.data, col = "sample", into = c("barcode","samplecode"), sep = "-")

#rename group names
Idents(SN3_MI14) <- "samplecode"
table(SN3_MI14$orig.ident)
SN3_MI14 <- RenameIdents(SN3_MI14,"1"="SN3_1", "2"="SN3_2", "3"="MI_1","4"="MI_4")
SN3_MI14$sample <- Idents(SN3_MI14)
SN3_MI14@meta.data <- separate(SN3_MI14@meta.data, col = "sample", into = c("group","sample_number"), sep = "_")

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
library(glmGamPoi)
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
#number of singlets and doublets
table(SN3_MI14_subset@meta.data$DF.classifications_0.25_0.3_3666)
SN3_MI14_subset <- subset(x = SN3_MI14_subset, subset= DF.classifications_0.25_0.3_3666 == "Singlet" )
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

#cell type annotation
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



#We identified that Fibroblasts to have most p16+ cells using snMultiome dataset only, now try to further analyze Fibroblasts
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
#annotation of Fibroblasts subtypes
SN3_MI14_Fb_subset_renameIdent <- SN3_MI14_Fb_subset
SN3_MI14_Fb_subset_renameIdent <- RenameIdents(SN3_MI14_Fb_subset_renameIdent,  "5"="Uninjured QuiCF","7"="Acute QuiCF",  "8"="Early ActCF","4"="Early MF","0"="Late MF","1"="Late MF","3"="Late MF", "6"="MaF","9"="MaF", "2"="Chronic QuiCF")
SN3_MI14_Fb_subset_renameIdent$subcelltype <- Idents(SN3_MI14_Fb_subset_renameIdent)

#visulize the p16+cells distribution in Fb

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



#Strategy 1: projecting spatial spider+ and spider- from half heart to snMultiome dataset

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

#ATAC analysis
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
DefaultAssay(SN3_MI14_subset_renameIdent)<-"peaks"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

SN3_MI14_subset_renameIdent<- AddMotifs(
  object = SN3_MI14_subset_renameIdent,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

library(chromVAR)
SN3_MI14_subset_renameIdent <- RunChromVAR(
  object = SN3_MI14_subset_renameIdent,
  genome = BSgenome.Mmusculus.UCSC.mm10
)



DefaultAssay(SN3_MI14_subset_renameIdent)<-"peaks"
Idents(SN3_MI14_subset_renameIdent) <- "senescence_status"

da_peaks <- FindMarkers(
  object = SN3_MI14_subset_renameIdent,
  ident.1 = "p16-Positive",
  ident.2 = "p16-Negative",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
write.csv(da_peaks, file = "Fb_p16pvsp16n_peaks_markers.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

SN3_MI14_subset_renameIdent.enriched.motifs <- FindMotifs(
  object = SN3_MI14_subset_renameIdent,
  features = top.da.peak
)
SN3_MI14_subset_renameIdent.enriched.motifs <-SN3_MI14_subset_renameIdent.enriched.motifs[rev(order(SN3_MI14_subset_renameIdent.enriched.motifs$fold.enrichment)),]

MotifPlot(
  object = SN3_MI14_subset_renameIdent,
  motifs = head(rownames(SN3_MI14_subset_renameIdent.enriched.motifs), n = 30),
)
