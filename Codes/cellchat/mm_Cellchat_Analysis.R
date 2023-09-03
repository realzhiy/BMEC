library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(multtest)
library(metap)

setwd("~/Desktop/Applied_Genomics/BMEC/symphony")

# load mouse seurat object
mouse <- readRDS("mouse_combined_annotated.rds")

# in the original data dbdb is labelled as -1
## Idents(mouse) <- grepl("-1",colnames(mouse))    

dbdb_cells=colnames(mouse)[(grep ("-1",colnames(mouse)))]
dbm_cells=colnames(mouse)[(grep ("-2",colnames(mouse)))]                     

cellchat_dbdb <- subset(x = mouse, cells = dbdb_cells)
cellchat_dbm <- subset(x = mouse, cells = dbm_cells)

# Cell chat
## Create cellchat object
DefaultAssay(cellchat_dbm) <- "RNA"
cellchat.dbm <- createCellChat(cellchat_dbm, group.by = "ident")

# Setting ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB# simply use the default CellChatDB
# set the used database in the object
cellchat.dbm@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat.dbm <- subsetData(cellchat.dbm) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat.dbm <- identifyOverExpressedGenes(cellchat.dbm)
cellchat.dbm <- identifyOverExpressedInteractions(cellchat.dbm)

## Compute the communication probability and infer cellular communication network
cellchat.dbm <- computeCommunProb(cellchat.dbm)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.dbm <- filterCommunication(cellchat.dbm, min.cells = 10)

# a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat.dbm)

# NB: The inferred intercellular communication network of each ligand-receptor 
# pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cellchat.dbm <- computeCommunProbPathway(cellchat.dbm)

## Calculate the aggregated cell-cell communication network
cellchat.dbm <- aggregateNet(cellchat.dbm)

setwd("~/Desktop/Applied_Genomics/BMEC/symphony/cc_output")
pdf("dbm_interaction_countandweight.pdf")
groupSize <- as.numeric(table(cellchat.dbm@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.dbm@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.dbm@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("dbm_edgeweights.pdf")
mat <- cellchat.dbm@net$weight
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2), xpd = TRUE)  # Adjusted plot margins
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

saveRDS(cellchat.dbm, file = "cellchat_dbm.rds")


## Create cellchat object
DefaultAssay(cellchat_dbdb) <- "RNA"
cellchat.dbdb <- createCellChat(cellchat_dbdb, group.by = "ident")

# Setting ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB# simply use the default CellChatDB
# set the used database in the object
cellchat.dbdb@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat.dbdb <- subsetData(cellchat.dbdb) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat.dbdb <- identifyOverExpressedGenes(cellchat.dbdb)
cellchat.dbdb <- identifyOverExpressedInteractions(cellchat.dbdb)

## Compute the communication probability and infer cellular communication network
cellchat.dbdb <- computeCommunProb(cellchat.dbdb)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.dbdb <- filterCommunication(cellchat.dbdb, min.cells = 10)

# a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat.dbdb)

# NB: The inferred intercellular communication network of each ligand-receptor 
# pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cellchat.dbdb <- computeCommunProbPathway(cellchat.dbdb)

## Calculate the aggregated cell-cell communication network
cellchat.dbdb <- aggregateNet(cellchat.dbdb)

pdf("dbdb_circle_interaction_countandweight.pdf")
groupSize <- as.numeric(table(cellchat.dbdb@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.dbdb@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.dbdb@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("dbdb_edgeweights.pdf")
mat <- cellchat.dbdb@net$weight
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2), xpd = TRUE)  # Adjusted plot margins
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat.dbdb <- netAnalysis_computeCentrality(cellchat.dbdb)
cellchat.dbm <- netAnalysis_computeCentrality(cellchat.dbm)

saveRDS(cellchat.dbdb, file = "cellchat_dbdb.rds")



dbdb <- netAnalysis_computeCentrality(dbdb)
dbm <- netAnalysis_computeCentrality(dbm)



