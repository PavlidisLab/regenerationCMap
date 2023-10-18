library(R.utils)
library(rhdf5)
library(data.table)
library(cmapR)
library(magrittr)
library(dplyr)

# these will be replaced with download links
dir.create('data-raw/pav_data')

natMat = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/fc.matrix.RDS.XZ")
natInstances = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/final.assign.mapping.RDS.XZ")
nat = readRDS("/home/nlim/MDE/RScripts/DataFreeze/Packaged/LINCS/GSE92742/Interim/sample.mapping.RDS.XZ")
natGeneAnnots = readr::read_tsv('data-raw/lincs1000_data/Gene_Info.txt')
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

saveRDS(natInstances, 'data-raw/pav_data/pavInstances.rds')
saveRDS(natGeneAnnots, 'data-raw/pav_data/pavGeneAnnots.rds')


natRankMatrix = natMat %>% apply(2,frankv,order = -1)
rownames(natRankMatrix) = rownames(natMat)

saveRDS(natRankMatrix, 'data-raw/pav_data/pavRankMatrix.rds')