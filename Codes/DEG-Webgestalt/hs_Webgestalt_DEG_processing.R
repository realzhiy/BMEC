## Subset Seurat object to include only cells of the cell type of interest
microglia_subset <- subset(x = control.combined.wang, idents = "Microglia")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(microglia_subset) <- microglia_subset$orig.ident


microglia_cardiac_markers <- FindMarkers(microglia_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

microglia_noncardiac_markers <- FindMarkers(microglia_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(microglia_cardiac_markers, 
	file = "Microglia_cardiac_markers.txt", sep = "\t")
write.table(microglia_noncardiac_markers, 
	file = "Microglia_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(microglia_subset, slot = "data"))[colSums(GetAssayData(microglia_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Microglia_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
Oligo_subset <- subset(x = control.combined.wang, idents = "Oligo")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Oligo_subset) <- Oligo_subset$orig.ident


Oligo_cardiac_markers <- FindMarkers(Oligo_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Oligo_noncardiac_markers <- FindMarkers(Oligo_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Oligo_cardiac_markers, 
	file = "Oligo_cardiac_markers.txt", sep = "\t")
write.table(Oligo_noncardiac_markers, 
	file = "Oligo_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Oligo_subset, slot = "data"))[colSums(GetAssayData(Oligo_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Oligo_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
Excitatory_subset <- subset(x = control.combined.wang, idents = "Excitatory")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Excitatory_subset) <- Excitatory_subset$orig.ident


Excitatory_cardiac_markers <- FindMarkers(Excitatory_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Excitatory_noncardiac_markers <- FindMarkers(Excitatory_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Excitatory_cardiac_markers, 
	file = "Excitatory_cardiac_markers.txt", sep = "\t")
write.table(Excitatory_noncardiac_markers, 
	file = "Excitatory_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Excitatory_subset, slot = "data"))[colSums(GetAssayData(Excitatory_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Excitatory_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
Astro_subset <- subset(x = control.combined.wang, idents = "Astro")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Astro_subset) <- Astro_subset$orig.ident


Astro_cardiac_markers <- FindMarkers(Astro_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Astro_noncardiac_markers <- FindMarkers(Astro_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Astro_cardiac_markers, 
	file = "Astro_cardiac_markers.txt", sep = "\t")
write.table(Astro_noncardiac_markers, 
	file = "Astro_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Astro_subset, slot = "data"))[colSums(GetAssayData(Astro_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Astro_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)



## Subset Seurat object to include only cells of the cell type of interest
Inhibitory_subset <- subset(x = control.combined.wang, idents = "Inhibitory")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Inhibitory_subset) <- Inhibitory_subset$orig.ident


Inhibitory_cardiac_markers <- FindMarkers(Inhibitory_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Inhibitory_noncardiac_markers <- FindMarkers(Inhibitory_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Inhibitory_cardiac_markers, 
	file = "Inhibitory_cardiac_markers.txt", sep = "\t")
write.table(Inhibitory_noncardiac_markers, 
	file = "Inhibitory_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Inhibitory_subset, slot = "data"))[colSums(GetAssayData(Inhibitory_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Inhibitory_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
OPC_subset <- subset(x = control.combined.wang, idents = "OPC")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(OPC_subset) <- OPC_subset$orig.ident


OPC_cardiac_markers <- FindMarkers(OPC_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

OPC_noncardiac_markers <- FindMarkers(OPC_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(OPC_cardiac_markers, 
	file = "OPC_cardiac_markers.txt", sep = "\t")
write.table(OPC_noncardiac_markers, 
	file = "OPC_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(OPC_subset, slot = "data"))[colSums(GetAssayData(OPC_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "OPC_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)



## Subset Seurat object to include only cells of the cell type of interest
Per_subset <- subset(x = control.combined.wang, idents = "Per")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Per_subset) <- Per_subset$orig.ident


Per_cardiac_markers <- FindMarkers(Per_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Per_noncardiac_markers <- FindMarkers(Per_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Per_cardiac_markers, 
	file = "Per_cardiac_markers.txt", sep = "\t")
write.table(Per_noncardiac_markers, 
	file = "Per_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Per_subset, slot = "data"))[colSums(GetAssayData(Per_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Per_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Subset Seurat object to include only cells of the cell type of interest
Endo_subset <- subset(x = control.combined.wang, idents = "Endo")

t2dm <- c("PDC05")
cardiac <- c("PDC34", "PDC87", "PDC91")
noncardiac <- c("PDC22", "C36", "C48")

Idents(Endo_subset) <- Endo_subset$orig.ident


Endo_cardiac_markers <- FindMarkers(Endo_subset, 
	ident.1 = cardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)

Endo_noncardiac_markers <- FindMarkers(Endo_subset, 
	ident.1 = noncardiac,
	ident.2 = t2dm,
	test.use = "bimod",
	subset.ident = "orig.ident",
	logfc.threshold = 0.25,
	verbose = FALSE)


# Print the list of markers
write.table(Endo_cardiac_markers, 
	file = "Endo_cardiac_markers.txt", sep = "\t")
write.table(Endo_noncardiac_markers, 
	file = "Endo_noncardiac_markers.txt", sep = "\t")

# Get all expressed gene
expressed_genes <- rownames(GetAssayData(Endo_subset, slot = "data"))[colSums(GetAssayData(Endo_subset, slot = "data")) > 0]

write.table(expressed_genes, file = "Endo_allgene.txt", sep = "\t", quote = FALSE, row.names = FALSE)




# Read the original txt file
Endo_cardiac_markers <- read.table("Endo_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Endo_up_filtered_genes <- Endo_cardiac_markers[Endo_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Endo_up_filtered_genes, file = "Endo_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Endo_down_filtered_genes <- Endo_cardiac_markers[Endo_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Endo_down_filtered_genes, file = "Endo_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Endo_noncardiac_markers <- read.table("Endo_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Endo_up_filtered_genes <- Endo_noncardiac_markers[Endo_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Endo_up_filtered_genes, file = "Endo_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Endo_down_filtered_genes <- Endo_noncardiac_markers[Endo_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Endo_down_filtered_genes, file = "Endo_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)








# Read the original txt file
Astro_cardiac_markers <- read.table("Astro_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Astro_up_filtered_genes <- Astro_cardiac_markers[Astro_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Astro_up_filtered_genes, file = "Astro_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Astro_down_filtered_genes <- Astro_cardiac_markers[Astro_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Astro_down_filtered_genes, file = "Astro_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Astro_noncardiac_markers <- read.table("Astro_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Astro_up_filtered_genes <- Astro_noncardiac_markers[Astro_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Astro_up_filtered_genes, file = "Astro_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Astro_down_filtered_genes <- Astro_noncardiac_markers[Astro_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Astro_down_filtered_genes, file = "Astro_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)










# Read the original txt file
Excitatory_cardiac_markers <- read.table("Excitatory_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Excitatory_up_filtered_genes <- Excitatory_cardiac_markers[Excitatory_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Excitatory_up_filtered_genes, file = "Excitatory_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Excitatory_down_filtered_genes <- Excitatory_cardiac_markers[Excitatory_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Excitatory_down_filtered_genes, file = "Excitatory_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Excitatory_noncardiac_markers <- read.table("Excitatory_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Excitatory_up_filtered_genes <- Excitatory_noncardiac_markers[Excitatory_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Excitatory_up_filtered_genes, file = "Excitatory_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Excitatory_down_filtered_genes <- Excitatory_noncardiac_markers[Excitatory_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Excitatory_down_filtered_genes, file = "Excitatory_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)








# Read the original txt file
Inhibitory_cardiac_markers <- read.table("Inhibitory_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Inhibitory_up_filtered_genes <- Inhibitory_cardiac_markers[Inhibitory_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Inhibitory_up_filtered_genes, file = "Inhibitory_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Inhibitory_down_filtered_genes <- Inhibitory_cardiac_markers[Inhibitory_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Inhibitory_down_filtered_genes, file = "Inhibitory_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Inhibitory_noncardiac_markers <- read.table("Inhibitory_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Inhibitory_up_filtered_genes <- Inhibitory_noncardiac_markers[Inhibitory_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Inhibitory_up_filtered_genes, file = "Inhibitory_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Inhibitory_down_filtered_genes <- Inhibitory_noncardiac_markers[Inhibitory_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Inhibitory_down_filtered_genes, file = "Inhibitory_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)








# Read the original txt file
Microglia_cardiac_markers <- read.table("Microglia_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Microglia_up_filtered_genes <- Microglia_cardiac_markers[Microglia_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Microglia_up_filtered_genes, file = "Microglia_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Microglia_down_filtered_genes <- Microglia_cardiac_markers[Microglia_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Microglia_down_filtered_genes, file = "Microglia_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Microglia_noncardiac_markers <- read.table("Microglia_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Microglia_up_filtered_genes <- Microglia_noncardiac_markers[Microglia_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Microglia_up_filtered_genes, file = "Microglia_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Microglia_down_filtered_genes <- Microglia_noncardiac_markers[Microglia_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Microglia_down_filtered_genes, file = "Microglia_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)









# Read the original txt file
Oligo_cardiac_markers <- read.table("Oligo_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Oligo_up_filtered_genes <- Oligo_cardiac_markers[Oligo_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Oligo_up_filtered_genes, file = "Oligo_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Oligo_down_filtered_genes <- Oligo_cardiac_markers[Oligo_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Oligo_down_filtered_genes, file = "Oligo_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Oligo_noncardiac_markers <- read.table("Oligo_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Oligo_up_filtered_genes <- Oligo_noncardiac_markers[Oligo_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Oligo_up_filtered_genes, file = "Oligo_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Oligo_down_filtered_genes <- Oligo_noncardiac_markers[Oligo_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Oligo_down_filtered_genes, file = "Oligo_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)






# Read the original txt file
OPC_cardiac_markers <- read.table("OPC_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
OPC_up_filtered_genes <- OPC_cardiac_markers[OPC_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(OPC_up_filtered_genes, file = "OPC_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
OPC_down_filtered_genes <- OPC_cardiac_markers[OPC_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(OPC_down_filtered_genes, file = "OPC_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
OPC_noncardiac_markers <- read.table("OPC_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
OPC_up_filtered_genes <- OPC_noncardiac_markers[OPC_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(OPC_up_filtered_genes, file = "OPC_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
OPC_down_filtered_genes <- OPC_noncardiac_markers[OPC_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(OPC_down_filtered_genes, file = "OPC_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)






# Read the original txt file
Per_cardiac_markers <- read.table("Per_cardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Per_up_filtered_genes <- Per_cardiac_markers[Per_cardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Per_up_filtered_genes, file = "Per_up_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Per_down_filtered_genes <- Per_cardiac_markers[Per_cardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Per_down_filtered_genes, file = "Per_down_cardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the original txt file
Per_noncardiac_markers <- read.table("Per_noncardiac_markers.txt", 
	header = TRUE)

# Filter genes with avg_log2FC > 0.58
Per_up_filtered_genes <- Per_noncardiac_markers[Per_noncardiac_markers$avg_log2FC > 0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Per_up_filtered_genes, file = "Per_up_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)


# Filter genes with avg_log2FC > 0.58
Per_down_filtered_genes <- Per_noncardiac_markers[Per_noncardiac_markers$avg_log2FC < -0.58, "x"]

# Write the filtered gene names to a new txt file
write.table(Per_down_filtered_genes, file = "Per_down_noncardiac_markers.txt", 
	row.names = FALSE, col.names = FALSE, quote = FALSE)







