library(R.utils)
library(rhdf5)
library(data.table)
library(cmapR)
library(magrittr)
library(dplyr)

dir.create('data-raw/FWD_data',recursive = TRUE)



download.file("https://maayanlab.cloud/l1000fwd/download/CD_signatures_full_42809x22268.gctx",
              destfile = 'data-raw/FWD_data/CD_signatures_full_42809x22268.gctx',method = 'wget')
download.file("https://maayanlab.cloud/l1000fwd/download/Drugs_metadata.csv",
              destfile = 'data-raw/FWD_data/Drugs_metadata.csv',method=  'wget')
download.file("https://maayanlab.cloud/l1000fwd/download/Probes_full_metadata.csv",
              destfile = 'data-raw/FWD_data/Probes_full_metadata.csv',method=  'wget')
download.file("https://maayanlab.cloud/l1000fwd/download/CD_signature_metadata.csv",
              destfile = 'data-raw/FWD_data/CD_signature_metadata.csv',method=  'wget')

fwdCMAP = cmapR::parse.gctx('data-raw/FWD_data/CD_signatures_full_42809x22268.gctx')
instances = readr::read_csv('data-raw/FWD_data/Drugs_metadata.csv')
fwdMetadata = readr::read_csv('data-raw/FWD_data/CD_signature_metadata.csv')
fwdGeneAnnots = readr::read_csv('data-raw/FWD_data/Probes_full_metadata.csv')
L1000geneAnnots = readr::read_tsv('data-raw/lincs1000_data/Gene_Info.txt')


fwdMetadata %<>% mutate(pert_iname =  instances$pert_iname[match(pert_id,instances$pert_id)])
fwdMatrix = fwdCMAP@mat
L1000geneAnnots %<>% filter(pr_is_bing == 1)


goodGenes = fwdGeneAnnots$pr_gene_symbol %in% L1000geneAnnots$pr_gene_symbol
fwdGeneAnnots = fwdGeneAnnots[goodGenes,]
fwdMatrix = fwdMatrix[,goodGenes]


all(colnames(fwdMatrix) == fwdGeneAnnots$pr_id)
all(rownames(fwdMatrix) == fwdMetadata$sig_id)
fwdRanks = fwdMatrix %>% apply(1,frankv, order = -1)
rownames(fwdRanks) = colnames(fwdMatrix)
all(colnames(fwdRanks) == fwdMetadata$sig_id)


saveRDS(fwdRanks,'data-raw/FWD_data/FWDranks.rds')
saveRDS(fwdGeneAnnots,'data-raw/FWD_data/FWDgeneAnnots.rds')
saveRDS(fwdMetadata,'data-raw/FWD_data/FWDinstances.rds')

