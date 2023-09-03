library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gplots)
setwd("~/Desktop/Applied_Genomics/BMEC/WebGestalt")

mouse_BP_down <- fread("BP_down.csv")
mouse_BP_up <- fread("BP_up.csv")
mouse_KEGG_down <- fread("KEGG_down.csv")
mouse_KEGG_up <- fread("KEGG_up.csv")
mouse_Reactome_down <- fread("Reactome_down.csv")
mouse_Reactome_up <- fread("Reactome_up.csv")

setwd("~/Desktop/Applied_Genomics/BMEC/control/DEG_celltypes/output")

human_BP_up_cardiac <- fread("bp/down_cardiac.csv")
human_BP_up_noncardiac <- fread("bp/down_noncardiac.csv")
human_BP_down_cardiac <- fread("bp/up_cardiac.csv")
human_BP_down_noncardiac <- fread("bp/up_noncardiac.csv")

human_KEGG_up_cardiac <- fread("KEGG/down_cardiac.csv")
human_KEGG_up_noncardiac <- fread("KEGG/down_noncardiac.csv")
human_KEGG_down_cardiac <- fread("KEGG/up_cardiac.csv")
human_KEGG_down_noncardiac <- fread("KEGG/up_noncardiac.csv")

human_Reactome_up_cardiac <- fread("Reactome/down_cardiac.csv")
human_Reactome_up_noncardiac <- fread("Reactome/down_noncardiac.csv")
human_Reactome_down_cardiac <- fread("Reactome/up_cardiac.csv")
human_Reactome_down_noncardiac <- fread("Reactome/up_noncardiac.csv")


data_frames <- c("human_BP_up_cardiac", 
                 "human_BP_up_noncardiac", 
                 "human_BP_down_cardiac", 
                 "human_BP_down_noncardiac",
                 "human_KEGG_up_cardiac",
                 "human_KEGG_up_noncardiac",
                 "human_KEGG_down_cardiac",
                 "human_KEGG_down_noncardiac",
                 "human_Reactome_up_cardiac",
                 "human_Reactome_up_noncardiac",
                 "human_Reactome_down_cardiac",
                 "human_Reactome_down_noncardiac",
                 "mouse_BP_down",
                 "mouse_BP_up",
                 "mouse_KEGG_down",
                 "mouse_KEGG_up",
                 "mouse_Reactome_down",
                 "mouse_Reactome_up")

# Loop through each data frame and perform the select operation
for (i in data_frames) {
  if ("Endo" %in% colnames(get(i))) {
    new_df_name <- paste("Endo_", i, sep = "")
    assign(new_df_name, get(i) %>% select(Endo, description))
  } else {
    cat("Skipped:", i, "because it doesn't contain 'Endo' column.\n")
  }
}

Endo_data_frames <- c("Endo_human_Reactome_up_cardiac", 
                 "Endo_human_Reactome_down_cardiac", 
                 "Endo_human_Reactome_down_noncardiac",
                 "Endo_human_KEGG_up_cardiac",
                 "Endo_human_KEGG_up_noncardiac",
                 "Endo_human_KEGG_down_cardiac",
                 "Endo_human_KEGG_down_noncardiac",
                 "Endo_human_BP_up_cardiac",
                 "Endo_human_BP_up_noncardiac",
                 "Endo_human_BP_down_cardiac",
                 "Endo_human_BP_down_noncardiac",
                 "Endo_mouse_BP_down",
                 "Endo_mouse_BP_up",
                 "Endo_mouse_KEGG_down",
                 "Endo_mouse_Reactome_down")

# Loop through each data frame name and rename "Endo" column
for (df_name in Endo_data_frames) {
  if (exists(df_name)) {  # Check if data frame exists
    column_to_rename <- "Endo"
    new_column_name <- "FDR"
    assign(df_name, get(df_name) %>% rename(!!new_column_name := !!column_to_rename))
  } else {
    cat("Data frame", df_name, "not found.\n")
  }
}

Endo_human_BP_up_cardiac$Sample <- "human_up_cardiac"
Endo_human_BP_up_noncardiac$Sample <- "human_BP_up_noncardiac"
Endo_mouse_BP_up$Sample <- "mouse_BP_up"

Endo_human_BP_down_cardiac$Sample <- "human_BP_down_cardiac"
Endo_human_BP_down_noncardiac$Sample <- "human_BP_down_noncardiac"
Endo_mouse_BP_down$Sample <- "mouse_BP_down"

Endo_human_Reactome_down_cardiac$Sample <- "human_Reactome_down_cardiac"
Endo_human_Reactome_down_noncardiac$Sample <- "human_Reactome_down_noncardiac"
Endo_mouse_Reactome_down$Sample <- "mouse_Reactome_down"

Endo_human_KEGG_up_cardiac$Sample <- "human_KEGG_up_cardiac"
Endo_human_KEGG_up_noncardiac$Sample <- "human_KEGG_up_noncardiac"

Endo_human_KEGG_down_cardiac$Sample <- "human_KEGG_down_cardiac"
Endo_human_KEGG_down_noncardiac$Sample <- "human_KEGG_down_noncardiac"
Endo_mouse_KEGG_down$Sample <- "mouse_KEGG_down"

Endo_BP_up_data_frames <- list(Endo_human_BP_up_cardiac,
                               Endo_human_BP_up_noncardiac,
                               Endo_mouse_BP_up)

Endo_BP_down_data_frames <- list(Endo_human_BP_down_cardiac,
                                 Endo_human_BP_down_noncardiac,
                                 Endo_mouse_BP_down)

Endo_Reactome_down_data_frames <- list(Endo_human_Reactome_down_cardiac,
                                       Endo_human_Reactome_down_noncardiac,
                                       Endo_mouse_Reactome_down)


Endo_KEGG_up_data_frames <- list(Endo_human_KEGG_up_cardiac,
                                 Endo_human_KEGG_up_noncardiac)

Endo_KEGG_down_data_frames <- list(Endo_human_KEGG_down_cardiac,
                                   Endo_human_KEGG_down_noncardiac,
                                   Endo_mouse_KEGG_down)




Endo_BP_up <- do.call(rbind, Endo_BP_up_data_frames)

Endo_BP_up[is.na(Endo_BP_up)] <- 1

Endo_BP_up$log_FDR <- -log10(Endo_BP_up$FDR)

Endo_BP_up$Significant <- Endo_BP_up$FDR < 0.05

ggplot(Endo_BP_up, aes(x = description, y = Sample, fill = log_FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "deescription",
       y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



Endo_BP_down <- do.call(rbind, Endo_BP_down_data_frames)

Endo_BP_down[is.na(Endo_BP_down)] <- 1

Endo_BP_down$log_FDR <- -log10(Endo_BP_down$FDR)

Endo_BP_down$Significant <- Endo_BP_down$FDR < 0.05

ggplot(Endo_BP_down, aes(x = description, y = Sample, fill = log_FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "deescription",
       y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



Endo_KEGG_down <- do.call(rbind, Endo_KEGG_down_data_frames)

Endo_KEGG_down[is.na(Endo_KEGG_down)] <- 1

Endo_KEGG_down$log_FDR <- -log10(Endo_KEGG_down$FDR)

Endo_KEGG_down$Significant <- Endo_KEGG_down$FDR < 0.05

ggplot(Endo_KEGG_down, aes(x = description, y = Sample, fill = log_FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "deescription",
       y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(Endo_KEGG_down, aes(x = description, y = Sample)) +
  geom_tile(aes(fill = ifelse(Significant, log_FDR, NA)), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "description",
       y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









combined_table <- Endo_human_BP_up_cardiac %>%
  full_join(Endo_human_BP_up_noncardiac, by = "description") %>%
  full_join(Endo_mouse_BP_up, by = "description")

filtered_table <- combined_table %>%
  filter(!is.na(Endo_human_BP_up_cardiac) | !is.na(Endo_human_BP_up_noncardiac) | !is.na(Endo_mouse_BP_up))

BP_up_combined <- filtered_table %>%
  select(description, everything())
write.csv(BP_up_combined, file = "Endo_BP_up_combined.csv", row.names = FALSE)

BP_up_combined_numeric <- BP_up_combined[, -1]  # Remove the description column

# Replace NA with a specific value for visualization purposes (e.g., 1 for grey)
BP_up_combined_numeric[is.na(BP_up_combined_numeric)] <- 1
BP_up_combined_numeric <- -log10(BP_up_combined_numeric)

# Replace infinite values with a large value for visualization purposes
BP_up_combined_numeric[is.infinite(BP_up_combined_numeric)] <- max(BP_up_combined_numeric[!is.infinite(BP_up_combined_numeric)])
# Replace Inf with 20 in the entire data frame
BP_up_combined_numeric[is.infinite(BP_up_combined_numeric)] <- 20


# Create the heatmap
heatmap.2(as.matrix(BP_up_combined_numeric),
          Colv = FALSE,     # Turn off column dendrogram
          Rowv = FALSE,     # Turn off row dendrogram
          col = colorRampPalette(c("grey", "darkblue"))(100),
          na.color = "grey",
          trace = "none",
          main = "Heatmap of log10(FDR) values")

# Annotate the heatmap with description values
heatmap_description <- final_table$description
text(x = 0, y = seq_along(heatmap_description), labels = heatmap_description, srt = 0, adj = c(1, 0.5), xpd = TRUE)










data_frames <- c("human_Reactome_up_cardiac", 
                 "human_Reactome_up_noncardiac", 
                 "human_Reactome_down_cardiac", 
                 "human_Reactome_down_noncardiac",
                 "human_KEGG_up_cardiac",
                 "human_KEGG_up_noncardiac",
                 "human_KEGG_down_cardiac",
                 "human_KEGG_down_noncardiac",
                 "human_Reactome_up_cardiac",
                 "human_Reactome_up_noncardiac",
                 "human_Reactome_down_cardiac",
                 "human_Reactome_down_noncardiac",
                 "mouse_BP_down",
                 "mouse_BP_up",
                 "mouse_KEGG_down",
                 "mouse_KEGG_up",
                 "mouse_Reactome_down",
                 "mouse_Reactome_up")

