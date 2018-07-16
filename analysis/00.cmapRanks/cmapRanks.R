library(rhdf5)
library(data.table)
library(cmapR)
cmapData = cmapR::parse.gctx('/space/scratch/nlim/LINCS/GSE92742-Level5/Level5_Data.gctx')
print('read DE data')
gc()
print('garbage collection complete')
inst = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Sig_Info.txt')
L1000geneAnnots = readr::read_tsv('/space/scratch/nlim/LINCS/GSE92742-Level5/Gene_Info.txt')
print('read metadata')

# measured or nicely imputed


inst %<>% rename(cmap_name = pert_iname)
rownames(inst) = inst$sig_id

print('fiddling')


goodGenes = L1000geneAnnots$pr_is_bing==1
L1000geneAnnots = L1000geneAnnots[goodGenes,]
deMatrix = cmapData@mat[goodGenes,]

print('good Genes filtered')


rankMatrix = deMatrix %>% apply(2,frank)

print('ranks calcualted')

saveRDS(rankMatrix,file = 'analysis/00.cmapRanks/rankMatrix.rds')
saveRDS(L1000geneAnnots,file = 'analysis/00.cmapRanks/rankMatrix.rds')
saveRDS(inst,file = 'analysis/00.cmapRanks/instances.rds')

print('data saved')

# data.table::frankv