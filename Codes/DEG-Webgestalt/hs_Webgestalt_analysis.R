#load packages
library(WebGestaltR)
#setting 
setwd("/rds/general/user/zw4419/home/control/Control_Samples_PD/DEG_celltypes")



celltype=c(
"Astro_up_cardiac",       
"Astro_down_cardiac",
"Endo_up_cardiac",    
"Endo_down_cardiac",
"Excitatory_up_cardiac",      
"Excitatory_down_cardiac",        
"Inhibitory_up_cardiac",
"Inhibitory_down_cardiac",         
"Microglia_up_cardiac",
"Microglia_down_cardiac",           
"Oligo_up_cardiac",
"Oligo_down_cardiac",
"Per_up_cardiac",
"Per_down_cardiac",
"OPC_up_cardiac",
"OPC_down_cardiac",
"Astro_up_noncardiac",       
"Astro_down_noncardiac",
"Endo_up_noncardiac",    
"Endo_down_noncardiac",
"Excitatory_up_noncardiac",      
"Excitatory_down_noncardiac",        
"Inhibitory_up_noncardiac",
"Inhibitory_down_noncardiac",         
"Microglia_up_noncardiac",
"Microglia_down_noncardiac",           
"Oligo_up_noncardiac",
"Oligo_down_noncardiac",
"Per_up_noncardiac",
"Per_down_noncardiac",
"OPC_up_noncardiac",
"OPC_down_noncardiac"
)

# Empty list to store the data frames
data_list <- list()

# Loop through each celltype to read in the data
for (name in celltype) {
  # Construct the filename
  file_name <- paste0(name, "_markers.txt")
  
  # Load the data
  data_list[[name]] <- read.table(file_name, header=TRUE, stringsAsFactors=FALSE)
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
outputDirectory="~/Desktop/Applied_Genomics/BMEC/control/DEG_celltypes/output"

WebGestaltR(enrichMethod="ORA", 
    organism="hsapiens", 
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
    organism="hsapiens", 
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
    organism="hsapiens", 
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


