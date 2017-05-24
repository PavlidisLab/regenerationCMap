library(data.table)
library(ogbox)
library(readxl)
library(dplyr)
library(magrittr)
library(devtools)
library(purrr)
# "complete genelist" is taken from erna directly
options(java.parameters = "-Xmx1024m")

unzip('data-raw/2013-092_CST_IPonly_2017_pipeline.zip',exdir = 'data-raw/ErnaData')

genesVoomLimma = read_excel('data-raw/ErnaData/DEA/Complete_geneList_voom_limma.xlsx')
genesLimma = read_excel('data-raw/ErnaData/DEA/Complete_geneList_limma.xlsx')
genesEdger = read_excel('data-raw/ErnaData/DEA/Complete_geneList_EdgeR.xlsx')

genesVoomLimmaRUVrk1 = read_excel('data-raw/ErnaData/DEA_RUVrk1/Complete_geneList_voom_limma.xlsx')
genesLimmaRUVrk1 = read_excel('data-raw/ErnaData/DEA_RUVrk1/Complete_geneList_limma.xlsx')
genesEdgerRUVrk1 = read_excel('data-raw/ErnaData/DEA_RUVrk1/Complete_geneList_EdgeR.xlsx')

genesVoomLimmaNoOutlier = read_excel('data-raw/ErnaData/DEA_wo_3day_naive_IP1//Complete_geneList_voom_limma.xlsx')
genesLimmaNoOutlier = read_excel('data-raw/ErnaData/DEA_wo_3day_naive_IP1/Complete_geneList_limma.xlsx')
genesEdgerNoOutlier = read_excel('data-raw/ErnaData/DEA_wo_3day_naive_IP1/Complete_geneList_EdgeR.xlsx')

use_data(genesVoomLimma)
use_data(genesLimma)
use_data(genesEdger)

use_data(genesVoomLimmaRUVrk1)
use_data(genesLimmaRUVrk1)
use_data(genesEdgerRUVrk1)

use_data(genesVoomLimmaNoOutlier)
use_data(genesLimmaNoOutlier)
use_data(genesEdgerNoOutlier)

# sigGenesVoomLimma =  read_excel('data-raw/ErnaData/DEA/Significant_geneList_voom_limma_FDR0.1.xlsx')


# hitlist = fread('data-raw/Complete_geneList_voom_limma.csv',data.table = FALSE)
# use_data(hitlist)

biocLite("ConnectivityMap")

download.file('https://portals.broadinstitute.org/cmap/msigdb_gene_sets.zip',destfile = 'data-raw/msigdb_gene_sets.zip')
unzip('data-raw/msigdb_gene_sets.zip',exdir = 'data-raw/')
MsigUp = readLines('data-raw/msigdb_up_mapped_to_HG_U133A.gmt') %>% strsplit('\t')
names = MsigUp %>% map_chr(1) %>% gsub(pattern = '_UP',replacement = '',x=.)

MsigDown = readLines('data-raw/msigdb_dn_mapped_to_HG_U133A.gmt')%>% strsplit('\t')

names2 = MsigDown %>% map_chr(1) %>% gsub(pattern = '_DN',replacement = '',x=.)

assertthat::see_if(all(names==names2))

MSigDB = lapply(1:length(MsigUp), function(i){
    upTags = MsigUp[[i]][c(-1,-2)]
    downTags = MsigDown[[i]][c(-1,-2)]
    return(list(upTags = upTags,downTags = downTags))
})
names(MSigDB) = names
devtools::use_data(MSigDB,overwrite = TRUE)

# # cmap downloads do not actually require login
# 
# # the first file is the ranked list of DE genes
# download.file('ftp://ftp.broad.mit.edu/pub/cmap/rankMatrix.txt.zip',destfile = 'data-raw/rankMatrix.txt.zip')
# unzip('data-raw/rankMatrix.txt.zip',exdir = 'data-raw/')
# rankMatrix = fread('data-raw/rankMatrix.txt',data.table=FALSE)
# # first line causes problems when reading. not sure why
# colnames = ogbox::checkLines('data-raw/rankMatrix.txt',lines =1) %>% {.[[1]]} %>% strsplit('\t') %>% {.[[1]]}
# cn(rankMatrix) = colnames
# use_data(rankMatrix)
# 
# # this file has the description of the datasets
# downoad.file('https://portals.broadinstitute.org/cmap/cmap_instances_02.xls', destfile = 'data-raw/cmap_instances_02.xls')
# 
# instances = loadWorkbook('data-raw/cmap_instances_02.xls')
# instances = readWorksheet(instances, sheet = 1, header = TRUE)
# 
# # last few lines are footnotes
# instances = instances[1:(nrow(instances)-9),]
# 
# geoDictionary = c("HG-U133A"= "GPL96",
#                   'HT_HG-U133A' = "GPL3921",
#                   'HT_HG-U133A_EA' = 'GPL3921')# they use an early access version of the chip with a 676 extra probes (https://portals.broadinstitute.org/cmap/help_topics_linkified.jsp). the probe data for this chip is not available
# 
# instances %<>% mutate(GPL = replaceElement(array3,geoDictionary) %>% {.$newVector})
# use_data(instances)


dir.create('data-raw/GemmaAnnots')
getGemmaAnnot('GPL96',chipFile = 'data-raw/GemmaAnnots/GPL96')
getGemmaAnnot('GPL3921',chipFile = 'data-raw/GemmaAnnots/GPL3921')



# save the randomized matrices for quick calculation
library(memoise)
library(ConnectivityMap)
load_all()
data("instances")
data("rankMatrix")
d = 100000
randomV = function(length,d){
    replicate(d,sample(x=1:ncol(rankMatrix),
                       size = length,replace = FALSE) %>% sort)
}
memoRandomV = memoise(randomV)
memoKsCalc = memoise(ksCalc)

chems = instances$cmap_name %>% unique
instanceLengths = chems %>% sapply(function(chem){
    chemInstances = rownames(instances)[instances$cmap_name %in% chem]
    length(chemInstances)
})

allVRandoms = instanceLengths %>% unique %>% sapply(function(x){
    print(x)
    return(memoRandomV(x,d))
})

allVRandoms %>% sapply(function(x){
    memoKsCalc(x,nrow(instances))
})

use_data(memoRandomV,overwrite = TRUE)
use_data(memoKsCalc,overwrite = TRUE)

ksPerm = ksCalc(Vrandoms,nrow(instances))
