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

rownames(rankMatrix) = rownames(deMatrix)

print('ranks calcualted')

saveRDS(rankMatrix,file = 'analysis/00.cmapRanks/rankMatrix.rds',compress= 'bzip2')
saveRDS(L1000geneAnnots,file = 'analysis/00.cmapRanks/L1000geneAnnots.rds')
saveRDS(inst,file = 'analysis/00.cmapRanks/instances.rds')

print('data saved')

# data.table::frankv




# # nathaniel's version
natMat = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/fc.matrix.RDS.XZ")
natInstances = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/final.assign.mapping.RDS.XZ")
nat = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/sample.mapping.RDS.XZ")
natGeneAnnots = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Gene_Info.txt')
natGeneAnnots = natGeneAnnots[match(rownames(natMat),natGeneAnnots$pr_gene_id),]
gc()


goodGenes = natGeneAnnots$pr_is_bing==1
natGeneAnnots = natGeneAnnots[goodGenes,]
natMat = natMat[goodGenes,]

natInstanceInfo = natInstances$trt_full %>% stringr::str_split(' ')

chem = natInstances$trt_full %>% stringr::str_extract('(?<= |^)[^ ]*(?= for)')
doses =  natInstances$trt_full %>% stringr::str_extract('[^ ]*(?= .*? .*? .*? .*? )')
cellLine =  natInstances$trt_full %>% stringr::str_extract('(?<= in ).*')
time = natInstances$trt_full %>% stringr::str_extract('(?<= for ).*(?= in )')

natInstances$chem = chem
natInstances$dose = doses
natInstances$cellLine = cellLine
natInstances$time = time

saveRDS(natInstances, 'analysis/00.cmapRanks/NatInstances.rds')
saveRDS(natGeneAnnots, 'analysis/00.cmapRanks/NatGeneAnnots.rds')


natRankMatrix = natMat %>% apply(2,frankv,order = -1)
rownames(natRankMatrix) = rownames(natMat)

saveRDS(natRankMatrix, 'analysis/00.cmapRanks/NatRankMatrix.rds')


# fwd version
fwdCMAP = cmapR::parse.gctx('data-raw/CD_signatures_full_42809x22268.gctx')
gc()
fwdMetadata = readr::read_csv('data-raw/CD_signature_metadata.csv')
fwdMetadata %<>% mutate(pert_iname =  instances$pert_iname[match(pert_id,instances$pert_id)])
fwdGeneAnnots = readr::read_csv('data-raw/Probes_full_metadata.csv')
fwdMatrix = fwdCMAP@mat
instances = readr::read_csv('data-raw/Drugs_metadata.csv')

L1000geneAnnots = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Gene_Info.txt')
L1000geneAnnots %<>% filter(pr_is_bing == 1)

goodGenes = fwdGeneAnnots$pr_gene_symbol %in% L1000geneAnnots$pr_gene_symbol
fwdGeneAnnots = fwdGeneAnnots[goodGenes,]
fwdMatrix = fwdMatrix[,goodGenes]


all(colnames(fwdMatrix) == fwdGeneAnnots$pr_id)
all(rownames(fwdMatrix) == fwdMetadata$sig_id)
fwdRanks = fwdMatrix %>% apply(1,frankv, order = -1)
rownames(fwdRanks) = colnames(fwdMatrix)
all(colnames(fwdRanks) == fwdMetadata$sig_id)


saveRDS(fwdMatrix,'analysis/00.cmapRanks/FWDranks.rds')
saveRDS(fwdGeneAnnots,'analysis/00.cmapRanks/FWDgeneAnnots.rds')
saveRDS(fwdMetadata,'analysis/00.cmapRanks/FWDinstances.rds')

