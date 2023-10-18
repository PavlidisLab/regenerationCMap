library(R.utils)
library(rhdf5)
library(data.table)
library(cmapR)
library(magrittr)
library(dplyr)


dir.create('data-raw/lincs1000_data',recursive = TRUE)



download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz",
              destfile = 'data-raw/lincs1000_data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz',method = 'wget')
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz",
              destfile = 'data-raw/lincs1000_data/GSE92742_Broad_LINCS_sig_info.txt.gz',method=  'wget')
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz",
              destfile = "data-raw/lincs1000_data/Gene_Info.txt",method = 'wget')

gunzip("data-raw/lincs1000_data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz")
gunzip("data-raw/lincs1000_data/GSE92742_Broad_LINCS_sig_info.txt.gz")


cmapData = cmapR::parse.gctx("data-raw/lincs1000_data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")

gc()

inst = readr::read_tsv('data-raw/lincs1000_data/GSE92742_Broad_LINCS_sig_info.txt.gz')
L1000geneAnnots = readr::read_tsv('data-raw/lincs1000_data/Gene_Info.txt')


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

saveRDS(rankMatrix,file = 'data-raw/lincs1000_data/rankMatrix.rds',compress= 'bzip2')
saveRDS(L1000geneAnnots,file = 'data-raw/lincs1000_data/L1000geneAnnots.rds')
saveRDS(inst,file = 'data-raw/lincs1000_data/instances.rds')