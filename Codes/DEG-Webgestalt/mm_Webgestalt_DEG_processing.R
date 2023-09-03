## Subset Seurat object to include only cells of the cell type of interest
microglia_subset <- subset(x = brain.combined.wang, idents = "Microglia")

Idents(microglia_subset) <- grepl("-1",colnames(microglia_subset))

microglia_db_markers <- FindMarkers(object = microglia_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(microglia_db_markers, file = "Microglia_db_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(microglia_subset, slot = "data"))[colSums(GetAssayData(microglia_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Microglia_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)




## Subset Seurat object to include only cells of the cell type of interest
Endo_subset <- subset(x = brain.combined.wang, idents = "Endo")

Idents(Endo_subset) <- grepl("-1",colnames(Endo_subset))

Endo_db_markers <- FindMarkers(object = Endo_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Endo_db_markers, file = "Endo_db_markers.txt", sep = "\t")

# Get all expressed gene
endo_expressed_genes <- rownames(GetAssayData(Endo_subset, slot = "data"))[colSums(GetAssayData(Endo_subset, slot = "data")) > 0]

write.table(endo_expressed_genes, file = "Endo_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)




## Subset Seurat object to include only cells of the cell type of interest
Oligo_subset <- subset(x = brain.combined.wang, idents = "Oligo")

Idents(Oligo_subset) <- grepl("-1",colnames(Oligo_subset))

Oligo_db_markers <- FindMarkers(object = Oligo_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Oligo_db_markers, file = "Oligo_db_markers.txt", sep = "\t")

# Get all expressed gene
oli_expressed_genes <- rownames(GetAssayData(Oligo_subset, slot = "data"))[colSums(GetAssayData(Oligo_subset, slot = "data")) > 0]

write.table(oli_expressed_genes, file = "Oligo_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)




## Subset Seurat object to include only cells of the cell type of interest
Excitatory_subset <- subset(x = brain.combined.wang, idents = "Excitatory")

Idents(Excitatory_subset) <- grepl("-1",colnames(Excitatory_subset))

Excitatory_db_markers <- FindMarkers(object = Excitatory_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Excitatory_db_markers, file = "Excitatory_db_markers.txt", sep = "\t")

# Get all expressed gene
ex_expressed_genes <- rownames(GetAssayData(Excitatory_subset, slot = "data"))[colSums(GetAssayData(Excitatory_subset, slot = "data")) > 0]

write.table(ex_expressed_genes, file = "Excitatory_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
Inhibitory_subset <- subset(x = brain.combined.wang, idents = "Inhibitory")

Idents(Inhibitory_subset) <- grepl("-1",colnames(Inhibitory_subset))

Inhibitory_db_markers <- FindMarkers(object = Inhibitory_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Inhibitory_db_markers, file = "Inhibitory_db_markers.txt", sep = "\t")

# Get all expressed gene
in_expressed_genes <- rownames(GetAssayData(Inhibitory_subset, slot = "data"))[colSums(GetAssayData(Inhibitory_subset, slot = "data")) > 0]

write.table(in_expressed_genes, file = "Inhibitory_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)





## Subset Seurat object to include only cells of the cell type of interest
Astro_subset <- subset(x = brain.combined.wang, idents = "Astro")

Idents(Astro_subset) <- grepl("-1",colnames(Astro_subset))

Astro_db_markers <- FindMarkers(object = Astro_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Astro_db_markers, file = "Astro_db_markers.txt", sep = "\t")

# Get all expressed gene
as_expressed_genes <- rownames(GetAssayData(Astro_subset, slot = "data"))[colSums(GetAssayData(Astro_subset, slot = "data")) > 0]

write.table(as_expressed_genes, file = "Astro_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)





## Subset Seurat object to include only cells of the cell type of interest
Per_subset <- subset(x = brain.combined.wang, idents = "Per")

Idents(Per_subset) <- grepl("-1",colnames(Per_subset))

Per_db_markers <- FindMarkers(object = Per_subset, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.25, test.use = "bimod")


# Print the list of markers
write.table(Per_db_markers, file = "Per_db_markers.txt", sep = "\t")

# Get all expressed gene
per_expressed_genes <- rownames(GetAssayData(Per_subset, slot = "data"))[colSums(GetAssayData(Per_subset, slot = "data")) > 0]

write.table(per_expressed_genes, file = "Per_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)













