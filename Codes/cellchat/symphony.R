setwd("/rds/general/user/zw4419/home/symphony")


#load required packages
library(Seurat)
library(harmony)
library(symphony)
library(ggplot2)
library(dplyr)
library(biomaRt) #change gene name
# loading required functions
source("/rds/general/user/zw4419/home/symphony/utils_seurat.R")


# Load the human and mouse Seurat object from the RDS file
mm_obj <- readRDS("mouse_combined_annotated.rds")
hs_obj <- readRDS("human_combined_annotated1.rds")

seurat_obj <- mm_obj


seurat_obj@meta.data$Names <- seurat_obj@active.ident
 
#deal with meta information of seurat object
mm_ref_seurat_obj=seurat_obj[ ,colnames(mm_obj)]
mm_ref_seurat_obj@meta.data$sub_celltype=mm_obj@active.ident
mm_ref_seurat_obj@meta.data$sample_id=as.character(mm_ref_seurat_obj@meta.data$orig.ident)
 


.verbose=FALSE
mm_ref_harmony <- mm_ref_seurat_obj %>%
                NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
                FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
                ScaleData(verbose = .verbose) %>% 
                RunPCA(verbose = .verbose) %>% 
                RunHarmony.Seurat("sample_id",verbose = .verbose) %>% 
                FindNeighbors(dims = 1:20, reduction = "harmony", verbose = .verbose) %>%  
                FindClusters(resolution = 0.5, verbose = .verbose)

 
mm_ref_harmony[["umap"]]=RunUMAP2(Embeddings(mm_ref_harmony, "harmony")[, 1:20], 
                          assay="RNA", verbose=FALSE, umap.method="uwot", return.model=TRUE)


 
 
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(mm_ref_harmony, reduction='umap', group.by='sub_celltype', shuffle = TRUE)
 
 
ref=buildReferenceFromSeurat(mm_ref_harmony, assay="RNA", verbose = TRUE, save_umap = TRUE, save_uwot_path = 'cache_symphony.uwot')
 
saveRDS(ref, '/rds/general/user/zw4419/home/symphony/testing_reference1.rds')
 
 
 
Map query (human)
 
 
#change human genes homologous to mouse
convertHumanGeneList=function(x){
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
 
 
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
mousex = genesV2
 
# Print the first 6 genes found to the screen
print(head(mousex))
return(mousex)
}
 
 
hsg_to_mmg=convertHumanGeneList(rownames(hs_obj@assays$RNA@counts))
 
map_count=hs_obj@assays$RNA@counts[hsg_to_mmg$HGNC.symbol, ]
rownames(map_count)=hsg_to_mmg$MGI.symbol
 
 
query=mapQuery(
    map_count, 
    hs_obj@meta.data, 
    ref,
    vars="orig.ident", 
    return_type="Seurat"
)
 
pdf("1.pdf")
DimPlot(mm_ref_harmony, reduction = 'umap', group.by='sub_celltype', label = TRUE, shuffle = TRUE) + labs(title = 'Original Reference (Mouse - Clusters)') 
DimPlot(query, reduction = 'umap', group.by='celltype', label = TRUE, shuffle = TRUE, raster=FALSE) + labs(title = 'Human - Mapped Query (Day)')
dev.off()
 
 
 
saveRDS(query, file = "/rds/general/user/zw4419/home/symphony/query.data.rds")
 
query <- readRDS("/rds/general/user/zw4419/home/symphony/query.data.rds") 
 
 
ref@meta.data$Names <- ref@meta.data$celltype
query <- knnPredict.Seurat(query, ref, 'Names', confidence = TRUE)
 
saveRDS(query, file = "/rds/general/user/zw4419/home/symphony/query_transfer_label_1.data.rds")
 
query_transfer_label_1 <- readRDS("/rds/general/user/zw4419/home/symphony/query_transfer_label_1.data.rds") 

