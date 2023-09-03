neuro_down_bp <- c("neurogenesis", 
  "neuron differenciation", 
  "axon development", 
  "central nervous system myelination", 
  "gliogenesis", 
  "glial cell differenciation", 
  "glial cell development",
  "plasma membrane bounded cell projection organization",
  "inflammatory response",
  "response to toxic substance",
  "oxidative phosphorylation",
  "response to oxidative stress",
  "neuron development",
  "neuron projection morphogenesis",
  "immune effector process",
  "immune response"
  )



library(dplyr)

# Filter the original data table to include only selected pathways
neuro_down_bp_subset <- mouse_BP_down %>%
  filter(description %in% neuro_down_bp)

# Identify missing pathways
missing_pathways <- setdiff(neuro_down_bp, subset_table$description)

# Create a data frame with missing pathways and fill p-values with 1 for all cell types
missing_rows <- data.frame(description = missing_pathways)
for (cell_type in colnames(subset_table)[2:8]) {
  missing_rows[[cell_type]] <- 1
}

# Combine the subset_table with the missing_rows data frame
final_subset_table <- bind_rows(subset_table, missing_rows)




mouse_cell_types <- c("Astro", "Endo", "Excitatory", "Inhibitory", "Microglia", "Oligo", "Per")
human_cell_types <- c(mouse_cell_types, "OPC")
all_cell_types <- c(mouse_cell_types, human_cell_types, human_cell_types) # Once for the mouse and twice for the two human datasets.

mouse_BP_down[is.na(mouse_BP_down)] <- 1
mouse_BP_down[is.na(human_BP_down_cardiac)] <- 1
mouse_BP_down[is.na(human_BP_down_noncardiac)] <- 1

# List all the pathways you're interested in 
all_pathways <- unique(c(mouse_BP_down$description, 
                         human_BP_down_cardiac$description, 
                         human_BP_down_noncardiac$description))

# Create an empty matrix
heatmap_matrix <- matrix(1, nrow = length(all_pathways), ncol = length(all_cell_types),
                         dimnames = list(all_pathways, all_cell_types))

mouse_BP_down <- as.data.frame(mouse_BP_down)
human_BP_down_cardiac <- as.data.frame(human_BP_down_cardiac)
human_BP_down_noncardiac <- as.data.frame(human_BP_down_noncardiac)

datasets <- list(mouse_BP_down, human_BP_down_cardiac, human_BP_down_noncardiac)
dataset_idx <- 1
for (data in datasets) {
  for (i in 1:nrow(data)) {
    pathway <- data$description[i]
    for (cell_type in colnames(data)[-1]) {
      heatmap_matrix[pathway, all_cell_types[(dataset_idx - 1) * length(colnames(data)[-1]) + which(colnames(data)[-1] == cell_type)]] <- data[i, cell_type]
    }
  }
  dataset_idx <- dataset_idx + 1
}

















upregulated human cardiac, 
"neurogenesis", 
"axonogenesis",
"blood vessel morphogenesis",
"epithelium development",
"negativer regulation of amyloid precursor protein biosynthetic process",
"apoptotic process",
"immune system process",
"generation of neurons",
"vasculature development",
  "positive regulation of tau-protein kinase activity",
  "regulation of response to stress",
  "response to oxidative stress",
  "angiogenesis",
  "endothelial cell differentiation",
  "Ras protein signal transduction"


