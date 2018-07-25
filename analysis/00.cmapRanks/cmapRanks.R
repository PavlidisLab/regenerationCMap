library(rhdf5)
library(data.table)
library(cmapR)
library(magrittr)
library(dplyr)
cmapData = cmapR::parse.gctx('/space/scratch/nlim/LINCS/GSE92742-Level5/Level5_Data.gctx')
print('read DE data')
gc()
print('garbage collection complete')
inst = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Sig_Info.txt')
L1000geneAnnots = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Gene_Info.txt')
print('read metadata')

instanceOrder = colnames(cmapData@mat)
geneOrder = rownames(cmapData@mat)

L1000geneAnnots = L1000geneAnnots[match(geneOrder,L1000geneAnnots$pr_gene_id),]
inst = inst[match(instanceOrder, inst$sig_id),]
# measured or nicely imputed
print('fiddling')


goodGenes = L1000geneAnnots$pr_is_bing==1
L1000geneAnnots = L1000geneAnnots[goodGenes,]
deMatrix = cmapData@mat[goodGenes,]

print('good Genes filtered')

goodChem = inst$pert_iname %>% table %>% {.>1} %>% which %>% names
goodInstances = inst$pert_iname %in% goodChem

inst = inst[goodInstances,]
deMatrix = deMatrix[,goodInstances]

print('Removed unreliable experiments')


rankMatrix = deMatrix %>% apply(2,frankv,order = -1)

print('ranks calcualted')

saveRDS(rankMatrix,file = 'analysis/00.cmapRanks/rankMatrix.rds',compress= 'bzip2')
saveRDS(L1000geneAnnots,file = 'analysis/00.cmapRanks/L1000geneAnnots.rds')
saveRDS(inst,file = 'analysis/00.cmapRanks/instances.rds')

print('data saved')

# data.table::frankv




# nathaniel's version
"/home/omancarci/allDirs/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/fc.matrix.RDS.XZ"
hede = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/final.assign.mapping.RDS.XZ")
hodo = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/sample.mapping.RDS.XZ")