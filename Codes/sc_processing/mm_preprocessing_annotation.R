setwd("~/Desktop/Applied_Genomics/BMEC/brain_ldsc")

library(Seurat)
library(dplyr)
library(patchwork)
library(multtest)
library(metap)

data3.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data)
# Filter mt

data3 <- CreateSeuratObject(counts = data3.data, project = "db",
                            min.cells = 3, min.features = 200)

data3[["percent.mt"]] <- PercentageFeatureSet(data3, pattern = "^mt-")
data3 <- subset(data3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
                & percent.mt < 20 & nCount_RNA > 400)

dbdb_cells=colnames(data3)[(grep ("-1",colnames(data3)))]
dbm_cells=colnames(data3)[(grep ("-2",colnames(data3)))]

data3_dbdb <- subset(x = data3, cells = dbdb_cells)
data3_dbm <- subset(x = data3, cells = dbm_cells)

dbdbbrain.list <- list(data3_dbdb, data3_dbm)


# normalize and identify variable features for each dataset independently
dbdbbrain.list <- lapply(X = dbdbbrain.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dbdbbrain.list)

brain.anchors <- FindIntegrationAnchors(object.list = dbdbbrain.list, anchor.features = features)

# this command creates an 'integrated' data assay
brain.combined <- IntegrateData(anchorset = brain.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(brain.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
brain.combined <- ScaleData(brain.combined, verbose = FALSE)
brain.combined <- RunPCA(brain.combined, npcs = 30, verbose = FALSE)
brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:30)
brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:30)
brain.combined <- FindClusters(brain.combined, resolution = 0.75)

# Visualization
p1 <- DimPlot(brain.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(brain.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(brain.combined, reduction = "umap", split.by = "orig.ident")

VlnPlot(brain.combined, features = c("Gpx1", "Cat"))

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(brain.combined) <- "RNA"


#Finding differential expressed features (biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
brain.combined.markers <- FindAllMarkers(brain.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
brain.combined.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
brain.combined.markers.LIST=unstack(brain.combined.markers, brain.combined.markers$gene ~ brain.combined.markers$cluster)


DER.UMI <- read.csv("DER-21_Single_cell_markergenes_UMI.csv", sep = ",", check.names = FALSE)
DER.UMI.LIST=unstack(DER.UMI, (DER.UMI$Gene) ~ DER.UMI$Cluster)
DER.UMI.LIST <- DER.UMI.LIST[-1]

# Create an object for wang annotation
brain.combined.wang=brain.combined
Filtered.Genes=unique(unlist(brain.combined.markers.LIST))

source("UseFulRfunction.R")

brain.combined.FisherTest.res.wang=lapply(brain.combined.markers.LIST,FisherTest.wang)

TMP.wang = matrix(unlist(brain.combined.FisherTest.res.wang), ncol = 25, byrow = TRUE)
colnames(TMP.wang)=names(DER.UMI.LIST)

TMP.LABELS.wang=CellTypeTest.wang(TMP.wang)
names(TMP.LABELS.wang)=names(brain.combined.FisherTest.res.wang)
wang.cluster.combined.ids = TMP.LABELS.wang
names(wang.cluster.combined.ids) = levels(brain.combined.wang)
brain.combined.wang = RenameIdents(brain.combined.wang,wang.cluster.combined.ids)
levels(brain.combined.wang)

pdf("annotated_umap")
DimPlot(brain.combined.wang, reduction = "umap", label = TRUE, pt.size = 1,label.size=5)
dev.off()

table(Idents(brain.combined.wang))

## Use FindMarkers to identify DEGs between the two conditions
## Here, we are comparing "dbm" to "dbdb"

# in the original data dbdb is labelled as -1
Idents(brain.combined.wang) <- grepl("-1",colnames(brain.combined.wang))    

i_dbdb_cells=colnames(brain.combined.wang)[(grep ("-1",colnames(brain.combined.wang)))]
i_dbm_cells=colnames(brain.combined.wang)[(grep ("-2",colnames(brain.combined.wang)))]                     

i_cellchat_dbdb <- subset(x = brain.combined.wang, cells = i_dbdb_cells)
i_cellchat_dbm <- subset(x = brain.combined.wang, cells = i_dbm_cells)

brain_db_DEG<- FindMarkers(object = brain.combined.wang, ident.1 = dbdb_cells, ident.2 = dbm_cells, logfc.threshold = 0.25, test.use = "bimod")
write.table(brain_db_DEG, file = "brain_DEG.txt", sep = "\t")


# Get all expressed gene
expressed_genes <- rownames(GetAssa\?yData(brain.combined.wang, slot = "data"))[colSums(GetAssayDat /b...xva(brain.combined.wang, slot = "data")) > 0]

write.table(expressed_genes, file = "brain_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)



saveRDS(brain.combined, file = "~/Desktop/Applied_Genomics/BMEC/brain_ldsc/mouse_combined_unannotated.rds")
saveRDS(brain.combined.wang, file = "~/Desktop/Applied_Genomics/BMEC/brain_ldsc/mouse_combined_annotated.rds")




