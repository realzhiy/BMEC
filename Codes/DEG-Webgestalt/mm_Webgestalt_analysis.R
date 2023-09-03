#load packages
library(WebGestaltR)
#setting 
setwd("/Users/zhiy/Desktop/Applied_Genomics/BMEC/brain_ldsc/DEG_celltypes/GO")



celltype=c(
"Astro_up",       
"Astro_down",
"Endo_up",    
"Endo_down",
"Excitatory_up",      
"Excitatory_down",        
"Inhibitory_up",
"Inhibitory_down",         
"Microglia_up",
"Microglia_down",           
"Oligo_up",
"Oligo_down",
"Per_up",
"Per_down"         
)


# Empty list to store the data frames
data_list_mouse <- list()

# Loop through each celltype to read in the data
for (name in celltype) {
  # Construct the filename
  file_name <- paste0(name, "_markers.txt")
  
  # Load the data
  data_list_mouse[[name]] <- read.table(file_name, header=TRUE, stringsAsFactors=FALSE)
}


for (i in 1:length(celltype))
{

#parameter
genedic="./genelist/"
genefile=paste0(genedic,celltype[i],"_markers.txt")

backdic="./background/"
backfile=paste0(backdic,celltype[i],"_allgene.txt")

cluster=celltype[i]

#run enrichemnt anlaysis
outputDirectory="/Users/zhiy/Desktop/Applied_Genomics/BMEC/brain_ldsc/DEG_celltypes/GO/output"

WebGestaltR(enrichMethod="ORA", 
    organism="mmusculus", 
    enrichDatabase="geneontology_Biological_Process", 
    interestGeneFile=genefile, 
    interestGeneType="genesymbol", 
    referenceGeneFile=backfile, 
    referenceGeneType="genesymbol", 
    minNum=5, 
    maxNum=2000, 
    sigMethod="fdr",
    isOutput=TRUE, 
    outputDirectory=outputDirectory, 
    projectName=paste0(cluster,'_bp'))

WebGestaltR(enrichMethod="ORA", 
    organism="mmusculus", 
    enrichDatabase="pathway_Reactome", 
    interestGeneFile=genefile, 
    interestGeneType="genesymbol", 
    referenceGeneFile=backfile, 
    referenceGeneType="genesymbol", 
    minNum=5, 
    maxNum=2000, 
    sigMethod="fdr",
    isOutput=TRUE, 
    outputDirectory=outputDirectory, 
    projectName=paste0(cluster,'_Reactome'))

WebGestaltR(enrichMethod="ORA", 
    organism="mmusculus", 
    enrichDatabase="pathway_KEGG", 
    interestGeneFile=genefile, 
    interestGeneType="genesymbol", 
    referenceGeneFile=backfile, 
    referenceGeneType="genesymbol", 
    minNum=5, 
    maxNum=2000, 
    sigMethod="fdr",
    isOutput=TRUE, 
    outputDirectory=outputDirectory, 
    projectName=paste0(cluster,'_KEGG'))

}


