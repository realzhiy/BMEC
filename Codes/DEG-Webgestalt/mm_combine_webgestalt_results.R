library(data.table)
library(dplyr)
library(tidyr)

setwd("~/Desktop/Applied_Genomics/BMEC/WebGestalt/KEGG/up")

Astro_up_KEGG <- fread("Project_Astro_up_KEGG/enrichment_results_Astro_up_KEGG.txt", header = TRUE, sep = "\t")
Endo_up_KEGG <- fread("Project_Endo_up_KEGG/enrichment_results_Endo_up_KEGG.txt", header = TRUE, sep = "\t")
Excitatory_up_KEGG <- fread("Project_Excitatory_up_KEGG/enrichment_results_Excitatory_up_KEGG.txt", header = TRUE, sep = "\t")
Inhibitory_up_KEGG <- fread("Project_Inhibitory_up_KEGG/enrichment_results_Inhibitory_up_KEGG.txt", header = TRUE, sep = "\t")
Microglia_up_KEGG <- fread("Project_Microglia_up_KEGG/enrichment_results_Microglia_up_KEGG.txt", header = TRUE, sep = "\t")
Oligo_up_KEGG <- fread("Project_Oligo_up_KEGG/enrichment_results_Oligo_up_KEGG.txt", header = TRUE, sep = "\t")
Per_up_KEGG <- fread("Project_Per_up_KEGG/enrichment_results_Per_up_KEGG.txt", header = TRUE, sep = "\t")


Astro <- Astro_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Astro$cell_type <- "Astro"

Endo <- Endo_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Endo$cell_type <- "Endo"

Excitatory <- Excitatory_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Excitatory$cell_type <- "Excitatory"

Inhibitory <- Inhibitory_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Inhibitory$cell_type <- "Inhibitory"

Microglia <- Microglia_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Microglia$cell_type <- "Microglia"

Oligo <- Oligo_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Oligo$cell_type <- "Oligo"


Per <- Per_up_KEGG %>%
  group_by(description) %>%
  summarise(min_FDR = min(FDR))

Per$cell_type <- "Per"



combined_table <- bind_rows(Astro, 
                            Excitatory, 
                            Inhibitory, 
                            Microglia, 
                            Per)

summary_table <- combined_table %>%
  group_by(description, cell_type) %>%
  summarise(min_FDR = min(min_FDR))


wide_summary_table <- summary_table %>%
  pivot_wider(names_from = cell_type, values_from = min_FDR, values_fill = NA)
write.csv(wide_summary_table, file = "up.csv", row.names = FALSE)


# Replace NA with 1
wide_summary_table[is.na(wide_summary_table)] <- 1

# The mtcars dataset:
data <- as.matrix(wide_summary_table[,-1])

data <- -log10(data)

# Default Heatmap
pdf("up_heatmap.pdf")
heatmap(data, xecCol=0.9, Colv = NA, Rowv = NA)
dev.off()










