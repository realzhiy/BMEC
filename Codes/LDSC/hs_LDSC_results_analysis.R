library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)

setwd("~/Desktop/Applied_Genomics/BMEC/final_LDSC/LDSC_Results")

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

summary_table <- Results.Table %>%
  select(FOLD, ct, Enrichment_p)

celltype = c("Astro_cardiac",
             "Astro_noncardiac",
             "Endo_cardiac",
             "Endo_noncardiac",
             "Excitatory_cardiac",
             "Excitatory_noncardiac",
             "Inhibitory_cardiac",
             "Inhibitory_noncardiac",
             "Microglia_cardiac",
             "Microglia_noncardiac",
             "Oligo_cardiac",
             "Oligo_noncardiac",
             "OPC_cardiac",
             "OPC_noncardiac",
             "Per_cardiac",
             "Per_noncardiac",
             "Astro",
             "Endo",
             "Excitatory",
             "Inhibitory",
             "Microglia",
             "Oligo",
             "Per")

adjusted_data <- list()

for (i in celltype) {
  filtered_data <- summary_table %>%
    filter(grepl(paste0(i,"$"), FOLD))
  
  filtered_data$adjusted_p <- p.adjust(filtered_data$Enrichment_p,method="BH")
  
  assign(i, filtered_data)
  
}

desired_order <- c("AD_2019",
                   "PD_2019",
                   "ANX_CC_2016",
                   "ANX_FS_2016",
                   "AN_2017",
                   "INS_2019",
                   "NEUR_2018",
                   "ADHD_2019",
                   "SCZ_2021",
                   "EPI_ALL_2018",
                   "EPI_FD_2018",
                   "EPI_GEN_2018",
                   "ASD_2017",
                   "MDD_2018",
                   "FTD_2017",
                   "ALS_2018",
                   "HV_2017",
                   "SD_2019",
                   "BD_2019",
                   "DS_2016",
                   "SWB_2016",
                   "INT_2018",
                   "EA_2018",
                   "2hGlu_2021",
                   "Glu_UKB",
                   "HbA1c_2021",
                   "T2D",
                   "FI_2021",
                   "FG_2021",
                   "LDLC_UKB",
                   "HDLC_UKB",
                   "RTC_UKB",
                   "CAD",
                   "IHD",
                   "IHD_UKB",
                   "CHD_UKB",
                   "AF_UKB",
                   "systolic_bp",
                   "diastolic",
                   "UKBB_COVID19",
                   "BMI_UKB",
                   "UKB_460K_body_HEIGHTz",
                   "UKB_460K_cov_EDU_YEARS",
                   "UKB_460K_pigment_HAIR",
                   "UKB_460K_pigment_SKIN")

data_frames <- list(Astro_cardiac, 
                    Endo_cardiac,  
                    Excitatory_cardiac, 
                    Inhibitory_cardiac, 
                    Microglia_cardiac, 
                    Oligo_cardiac, 
                    OPC_cardiac, 
                    Per_cardiac)
cardiac <- do.call(rbind, data_frames)

cardiac$log_Adjusted_p <- -log10(cardiac$adjusted_p)

cardiac$Significant <- cardiac$adjusted_p < 0.05

cardiac$FOLD <- gsub("_cardiac$", "", cardiac$FOLD)

cardiac <- cardiac %>%
  filter(ct %in% desired_order)

cardiac <- cardiac %>%
  arrange(factor(ct, levels = desired_order))


pdf("cardiac_heatmap.pdf", width = 13, height = 4)
ggplot(cardiac, aes(x = factor(ct, levels = desired_order), y = FOLD, fill = log_Adjusted_p)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("ivory", "darkblue", "indianred3"),
                       values = c(0, 0.5, 1),  # Specify the position of the second color
                       breaks = c(0, -log10(0.05), -log10(0)),
                       labels = c("Non-significant", "Significant", "Significant"),
                       name = "-log10(FDR)") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()

sig_sumstats <- unique(significant_cardiac$ct)
cardiac_sig <- cardiac %>%
  filter(ct %in% sig_sumstats)

cardiac_sig <- cardiac_sig[cardiac_sig$FOLD %in% c("Per", "Oligo", "Endo"), ]


pdf("sig_cardiac_heatmap.pdf", width = 13, height = 4)
ggplot(cardiac_sig, aes(x = factor(ct, levels = desired_order), y = FOLD)) +
  geom_tile(aes(fill = ifelse(Significant, log_Adjusted_p, NA)), color = "black") +
  scale_fill_gradient(low = "darkblue", high = "indianred3", na.value = "ivory", name = "-log10(FDR)") +
  labs(x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()


###########


data_frames <- list(Astro_noncardiac, 
                    Endo_noncardiac,  
                    Excitatory_noncardiac, 
                    Inhibitory_noncardiac, 
                    Microglia_noncardiac, 
                    Oligo_noncardiac, 
                    OPC_noncardiac, 
                    Per_noncardiac)
noncardiac <- do.call(rbind, data_frames)

noncardiac$log_Adjusted_p <- -log10(noncardiac$adjusted_p)

noncardiac$Significant <- noncardiac$adjusted_p < 0.05

noncardiac$FOLD <- gsub("_noncardiac$", "", noncardiac$FOLD)

noncardiac <- noncardiac %>%
  filter(ct %in% desired_order)

noncardiac <- noncardiac %>%
  arrange(factor(ct, levels = desired_order))


pdf("noncardiac_heatmap.pdf", width = 13, height = 4)
ggplot(noncardiac, aes(x = factor(ct, levels = desired_order), y = FOLD, fill = log_Adjusted_p)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("ivory", "darkblue", "indianred3"),
                       values = c(0, 0.5, 1),  # Specify the position of the second color
                       breaks = c(0, -log10(0.05), -log10(0)),
                       labels = c("Non-significant", "Significant", "Significant"),
                       name = "-log10(FDR)") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()

significant_noncardiac <- noncardiac[noncardiac$Significant, ]
sig_sumstats <- unique(significant_noncardiac$ct)
noncardiac_sig <- noncardiac %>%
  filter(ct %in% sig_sumstats)

noncardiac_sig <- noncardiac_sig[noncardiac_sig$FOLD %in% c("Per", "Endo"), ]


pdf("sig_noncardiac_heatmap.pdf", width = 13, height = 4)
ggplot(noncardiac_sig, aes(x = factor(ct, levels = desired_order), y = FOLD)) +
  geom_tile(aes(fill = ifelse(Significant, log_Adjusted_p, NA)), color = "black") +
  scale_fill_gradient(low = "darkblue", high = "indianred3", na.value = "ivory", name = "-log10(FDR)") +
  labs(x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()

########

data_frames <- list(Astro, 
                    Endo,  
                    Excitatory, 
                    Inhibitory, 
                    Microglia, 
                    Oligo,
                    Per)
mouse <- do.call(rbind, data_frames)

mouse$log_Adjusted_p <- -log10(mouse$adjusted_p)

mouse$Significant <- mouse$adjusted_p < 0.05

mouse$FOLD <- gsub("$", "", mouse$FOLD)

mouse <- mouse %>%
  filter(ct %in% desired_order)

mouse <- mouse %>%
  arrange(factor(ct, levels = desired_order))


pdf("mouse_heatmap.pdf", width = 13, height = 4)
ggplot(mouse, aes(x = factor(ct, levels = desired_order), y = FOLD, fill = log_Adjusted_p)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("ivory", "darkblue", "indianred3"),
                       values = c(0, 0.5, 1),  # Specify the position of the second color
                       breaks = c(0, -log10(0.05), -log10(0)),
                       labels = c("Non-significant", "Significant", "Significant"),
                       name = "-log10(FDR)") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()

significant_mouse <- mouse[mouse$Significant, ]
sig_sumstats <- unique(significant_mouse$ct)
mouse_sig <- mouse %>%
  filter(ct %in% sig_sumstats)

mouse_sig <- mouse_sig[mouse_sig$FOLD %in% c("Per", "Oligo", 
                                           "Inhibitory", "Excitatory", "Endo"), ]

pdf("sig_mouse_heatmap.pdf", width = 13, height = 4)
ggplot(mouse_sig, aes(x = factor(ct, levels = desired_order), y = FOLD)) +
  geom_tile(aes(fill = ifelse(Significant, log_Adjusted_p, NA)), color = "black") +
  scale_fill_gradient(low = "darkblue", high = "indianred3", na.value = "ivory", name = "-log10(FDR)") +
  labs(title = "Heatmap of -log10(FDR)",
       x = "GWAS SumStats",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank())
dev.off()
