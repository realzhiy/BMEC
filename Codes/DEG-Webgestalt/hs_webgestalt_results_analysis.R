library(data.table)

# load webgestalt results
Astro_up_noncardiac_Reactome <- fread("enrichment_results_Astro_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Endo_up_noncardiac_Reactome <- fread("enrichment_results_Endo_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Excitatory_up_noncardiac_Reactome <- fread("enrichment_results_Excitatory_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Inhibitory_up_noncardiac_Reactome <- fread("enrichment_results_Inhibitory_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Microglia_up_noncardiac_Reactome <- fread("enrichment_results_Microglia_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Oligo_up_noncardiac_Reactome <- fread("enrichment_results_Oligo_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
OPC_up_noncardiac_Reactome <- fread("enrichment_results_OPC_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")
Per_up_noncardiac_Reactome <- fread("enrichment_results_Per_up_noncardiac_Reactome.txt", header = TRUE, sep = "\t")

# organize the result table
Astro <- Astro_up_noncardiac_Reactome %>%
       group_by(description) %>%
       summarise(min_FDR = min(FDR))

Astro$cell_type <- "Astro"

Endo <- Endo_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Endo$cell_type <- "Endo"

Excitatory <- Excitatory_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Excitatory$cell_type <- "Excitatory"

Inhibitory <- Inhibitory_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Inhibitory$cell_type <- "Inhibitory"

Microglia <- Microglia_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Microglia$cell_type <- "Microglia"

Oligo <- Oligo_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Oligo$cell_type <- "Oligo"

OPC <- OPC_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

OPC$cell_type <- "OPC"

Per <- Per_up_noncardiac_Reactome %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Per$cell_type <- "Per"

# create a combined table with all cell type
combined_table <- bind_rows(Astro, 
                            Endo, 
                            Excitatory, 
                            Inhibitory, 
                            Microglia, 
                            Oligo, 
                            OPC, 
                            Per)

summary_table <- combined_table %>%
  group_by(description, cell_type) %>%
  summarise(min_FDR = min(min_FDR))

library(tidyr)
wide_summary_table <- summary_table %>%
       pivot_wider(names_from = cell_type, values_from = min_FDR, values_fill = NA)
write.csv(wide_summary_table, file = "up_noncardiac.csv", row.names = FALSE)


# Replace NA with 1
wide_summary_table[is.na(wide_summary_table)] <- 1

# create matrix
data <- as.matrix(wide_summary_table[,-1])

data <- -log10(data)

# Default Heatmap
pdf("up_noncardiac_heatmap.pdf")
heatmap(data, xecCol=0.9, Colv = NA, Rowv = NA)
dev.off()


