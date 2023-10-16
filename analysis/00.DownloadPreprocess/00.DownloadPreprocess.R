library(data.table)
library(ogbox)
library(readxl)
library(dplyr)
library(magrittr)
library(devtools)
library(purrr)
library(usethis)
library(rhdf5)
library(cmapR)
# "complete genelist" is taken from erna directly
options(java.parameters = "-Xmx1024m")

# processing differential expression raw data. ----
# unzip('data-raw/2013-092_CST_IPonly_2017_pipeline.zip',exdir = 'data-raw/dif_exp_data')

dif_exp_data = read_excel('data-raw/ErnaData/DEA_wo_3day_naive_IP1/Complete_geneList_EdgeR.xlsx')
usethis::use_data(dif_exp_data,overwrite = TRUE)

# historic download link for the MSigDB gene sets used by the original connectivity map
# this link is now removed but was accessed May 10th 2017.
# download.file('https://portals.broadinstitute.org/cmap/msigdb_gene_sets.zip',destfile = 'data-raw/msigdb_gene_sets2.zip')
unzip('data-raw/msigdb_gene_sets.zip',exdir = 'data-raw/')
MsigUp = readLines('data-raw/msigdb_up_mapped_to_HG_U133A.gmt') %>% strsplit('\t')
names = MsigUp %>% map_chr(1) %>% gsub(pattern = '_UP',replacement = '',x=.)

MsigDown = readLines('data-raw/msigdb_dn_mapped_to_HG_U133A.gmt')%>% strsplit('\t')


MSigDB = lapply(1:length(MsigUp), function(i){
    upTags = MsigUp[[i]][c(-1,-2)]
    downTags = MsigDown[[i]][c(-1,-2)]
    return(list(upTags = upTags,downTags = downTags))
})
names(MSigDB) = names
usethis::use_data(MSigDB,overwrite = TRUE)



dir.create('data-raw/GemmaAnnots')
gemmaAPI::getAnnotation('GPL96',file = 'data-raw/GemmaAnnots/GPL96',overwrite = TRUE)
# gemmaAPI::getAnnotation('GPL3921',file = 'data-raw/GemmaAnnots/GPL3921_2',overwrite = TRUE)




# save the randomized matrices for quick calculation --------
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

instanceLenghts = table(instances$cmap_name)


allVRandoms = instanceLengths %>% unique %>% sapply(function(x){
    print(x)
    return(memoRandomV(x,d))
})

allVRandoms %>% sapply(function(x){
    memoKsCalc(x,nrow(instances))
})

usethis::use_data(memoRandomV,overwrite = TRUE)
usethis::use_data(memoKsCalc,overwrite = TRUE)

ksPerm = ksCalc(Vrandoms,nrow(instances))

# precalculate msgibdb enrichment --------------

MSigDB_enrich = MSigDB %>% mclapply(function(signature){
    print(signature$name)
    upTags = signature$upTags 
    
    downTags = signature$downTags
    
    out = connectivityMapEnrichment(upTags,downTags, rankMatrix,instances,d = 100000)
    return(out$chemScores)
},mc.cores=16)
usethis::use_data(MSigDB_enrich)






