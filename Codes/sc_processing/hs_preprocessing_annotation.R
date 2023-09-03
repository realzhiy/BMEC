setwd("~/Desktop/Applied_Genomics/BMEC/control/Control_Samples_PD")
setwd("/media/sdc1-2TB/zw4419/Control_Samples_PD")
setwd("/rds/general/user/zw4419/home/control/Control_Samples_PD")

library(Seurat)
library(dplyr)
library(patchwork)
library(multtest)
library(metap)

C36.OE <- readRDS("C36.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
C48.OE <- readRDS("C48.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
PDC05.OE <- readRDS("PDC05.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
PDC22.OE <- readRDS("PDC22.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
PDC34.OE <- readRDS("PDC34.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
PDC87.OE <- readRDS("PDC87.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
PDC91.OE <- readRDS("PDC91.OE.NonNormalised.filteredcounts.dgCMatrix.rds")

# Initialize the Seurat object with the raw (non-normalized data)
# Filter mt

data1 <- CreateSeuratObject(counts = C36.OE, project = "C36",
                            min.cells = 3, min.features = 200)

data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^mt-")
data1 <- subset(data1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data2 <- CreateSeuratObject(counts = C48.OE, project = "C48",
                            min.cells = 3, min.features = 200)

data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^mt-")
data2 <- subset(data2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data3 <- CreateSeuratObject(counts = PDC05.OE, project = "PDC05",
                            min.cells = 3, min.features = 200)

data3[["percent.mt"]] <- PercentageFeatureSet(data3, pattern = "^mt-")
data3 <- subset(data3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data4 <- CreateSeuratObject(counts = PDC22.OE, project = "PDC22",
                            min.cells = 3, min.features = 200)

data4[["percent.mt"]] <- PercentageFeatureSet(data4, pattern = "^mt-")
data4 <- subset(data4, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data5 <- CreateSeuratObject(counts = PDC34.OE, project = "PDC34",
                            min.cells = 3, min.features = 200)

data5[["percent.mt"]] <- PercentageFeatureSet(data5, pattern = "^mt-")
data5 <- subset(data5, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data6 <- CreateSeuratObject(counts =PDC87.OE, project = "PDC87",
                            min.cells = 3, min.features = 200)

data6[["percent.mt"]] <- PercentageFeatureSet(data6, pattern = "^mt-")
data6 <- subset(data6, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

data7 <- CreateSeuratObject(counts = PDC91.OE, project = "PDC91",
                            min.cells = 3, min.features = 200)

data7[["percent.mt"]] <- PercentageFeatureSet(data7, pattern = "^mt-")
data7 <- subset(data7, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)



control.list <- list(data1, data2, data3, data4, data5, data6, data7)


# normalize and identify variable features for each dataset independently
control.list <- lapply(X = control.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = control.list)

control.anchors <- FindIntegrationAnchors(object.list = control.list, anchor.features = features)

saveRDS(control.anchors, file = "/rds/general/user/zw4419/home/control/Control_Samples_PD/control_anchors_seurat_object.rds")
saveRDS(control.anchors, file = "~/Desktop/Applied_Genomics/BMEC/control/Control_Samples_PD/control_anchors_seurat_object.rds")

# this command creates an 'integrated' data assay
control.combined <- IntegrateData(anchorset = control.anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(control.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
control.combined <- ScaleData(control.combined, verbose = FALSE)
control.combined <- RunPCA(control.combined, npcs = 30, verbose = FALSE)
control.combined <- RunUMAP(control.combined, reduction = "pca", dims = 1:30)
control.combined <- FindNeighbors(control.combined, reduction = "pca", dims = 1:30)
control.combined <- FindClusters(control.combined, resolution = 2)

DimPlot(control.combined, reduction = "umap", split.by = "orig.ident")

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(control.combined) <- "RNA"

saveRDS(control.combined, file = "/rds/general/user/zw4419/home/control/Control_Samples_PD/control_combined_seurat_object.rds")

#Finding differential expressed features (biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
control.combined.markers <- FindAllMarkers(control.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
control.combined.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
control.combined.markers.LIST=unstack(control.combined.markers, control.combined.markers$gene ~ control.combined.markers$cluster)

saveRDS(control.combined.markers, file = "/rds/general/user/zw4419/home/control/Control_Samples_PD/control_combined_markers.rds")

DER.UMI <- read.csv("DER-21_Single_cell_markergenes_UMI.csv", sep = ",", check.names = FALSE)
DER.UMI.LIST=unstack(DER.UMI, (DER.UMI$Gene) ~ DER.UMI$Cluster)


# Create an object for wang annotation
control.combined.wang=control.combined
Filtered.Genes=unique(unlist(control.combined.markers.LIST))

source("UseFulRfunction.R")

control.combined.FisherTest.res.wang=lapply(control.combined.markers.LIST,FisherTest.wang)

TMP.wang = matrix(unlist(control.combined.FisherTest.res.wang), ncol = 25, byrow = TRUE)
colnames(TMP.wang)=names(DER.UMI.LIST)

TMP.LABELS.wang=CellTypeTest.wang(TMP.wang)
names(TMP.LABELS.wang)=names(control.combined.FisherTest.res.wang)
wang.cluster.combined.ids = TMP.LABELS.wang
names(wang.cluster.combined.ids) = levels(control.combined.wang)
control.combined.wang = RenameIdents(control.combined.wang,wang.cluster.combined.ids)
levels(control.combined.wang)

pdf("annotated_umap.pdf")
DimPlot(control.combined.wang, reduction = "umap", label = TRUE, pt.size = 1,label.size=5)
dev.off()

table(Idents(control.combined.wang))


t2dm <- subset(control.combined.wang, subset = orig.ident == "PDC05")

# Subset Seurat object for "cardiac" containing PDC34, PDC87, and PDC22
cardiac <- subset(control.combined.wang, subset = orig.ident %in% c("PDC34", "PDC87", "PDC22"))

# Subset Seurat object for "noncardiac" containing PDC91, C36, and C48
noncardiac <- subset(control.combined.wang, subset = orig.ident %in% c("PDC91", "C36", "C48"))

## Use FindMarkers to identify DEGs between the two conditions
## Here, we are comparing "dbm" to "dbdb"

# in the original data dbdb is labelled as -1

