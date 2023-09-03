library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("~/Desktop/Applied_Genomics/BMEC/brain_ldsc/LDSC_Results")

FOLDERS=list.files()
count=0
for(FOLD in FOLDERS){
  All.files=list.files(path=FOLD,pattern="")
  for(file in All.files){
    data=fread(sprintf("%s/%s",FOLD,file))
    tmp=data[1,]
    tmp$FOLD=FOLD
    tmp$file=file
    if(count==0){
      Results.Table=tmp
    }else{
      Results.Table=rbind(Results.Table,tmp)
    }
    count=count+1
  }
}
Results.Table$ct=gsub("_\\w[^_]*.results","",Results.Table$file,perl=TRUE)
Results.Table$ct=gsub(".results","",Results.Table$file,perl=TRUE)
Results.Table$ct=gsub(paste(FOLDERS,collapse="|"),"",Results.Table$ct)
Results.Table$ct=gsub("_$","",Results.Table$ct)
Results.Table$ct=gsub("\\.ENTID\\.|\\.sumstats\\.gz", "", Results.Table$ct)
Results.Table$ct <- gsub("^microglia_up", "", Results.Table$ct)
Results.Table$ct <- gsub("^microglia_down", "", Results.Table$ct)

summary_table <- Results.Table %>%
  select(FOLD, ct, Enrichment_p)

celltype = c("Astro_down",
             "Astro_up",
             "Endo_down",
             "Endo_up",
             "Excitatory_down",
             "Excitatory_up",
             "Inhibitory_down",
             "Inhibitory_up",
             "Microglia_down",
             "Microglia_up",
             "Oligo_down",
             "Oligo_up",
             "Per_down",
             "Per_up")

adjusted_data <- list()

for (i in celltype) {
  filtered_data <- summary_table %>%
    filter(grepl(paste0(i,"$"), FOLD))
  
  filtered_data$adjusted_p <- p.adjust(filtered_data$Enrichment_p,method="BH")
  
  assign(i, filtered_data)
  
}


data_frames <- list(Astro_up, 
                    Endo_up,  
                    Excitatory_up, 
                    Inhibitory_up, 
                    Microglia_up, 
                    Oligo_up, 
                    Per_up)
up <- do.call(rbind, data_frames)

up$log_Adjusted_p <- -log10(up$adjusted_p)

up$Significant <- up$adjusted_p < 0.05

up$FOLD <- gsub("_up$", "", up$FOLD)

pdf("up_heatmap.pdf", width = 13, height = 4)
ggplot(up, aes(x = ct, y = FOLD, fill = log_Adjusted_p)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of -log10(Adjusted p-values)",
       x = "ct",
       y = "FOLD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sig_up_heatmap.pdf", width = 13, height = 4)
ggplot(up, aes(x = ct, y = FOLD)) +
  geom_tile(aes(fill = ifelse(Significant, log_Adjusted_p, NA)), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(title = "Heatmap of -log10(Adjusted p-values)",
       x = "ct",
       y = "FOLD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###########

data_frames <- list(Astro_down, 
                    Endo_down,  
                    Excitatory_down, 
                    Inhibitory_down, 
                    Microglia_down, 
                    Oligo_down, 
                    Per_down)
down <- do.call(rbind, data_frames)

down$log_Adjusted_p <- -log10(down$adjusted_p)

down$Significant <- down$adjusted_p < 0.05

down$FOLD <- gsub("_down$", "", down$FOLD)

pdf("down_heatmap.pdf", width = 13, height = 4)
ggplot(down, aes(x = ct, y = FOLD, fill = log_Adjusted_p)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of -log10(Adjusted p-values)",
       x = "ct",
       y = "FOLD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sig_down_heatmap.pdf", width = 13, height = 4)
ggplot(down, aes(x = ct, y = FOLD)) +
  geom_tile(aes(fill = ifelse(Significant, log_Adjusted_p, NA)), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(title = "Heatmap of -log10(Adjusted p-values)",
       x = "ct",
       y = "FOLD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
