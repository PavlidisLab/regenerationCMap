library(glue)
library(magrittr)
library(dplyr)
library(ogbox)
library(lazyeval)
library(parallel)
library(ConnectivityMap)
library(memoise)
library(MSigDB)
data("rankMatrix")
data("instances")
devtools::load_all()

for(i in 1:length(MSigDB)){
    MSigDB[[i]]$name = names(MSigDB)[i]
}

MSigDB_enrich = MSigDB %>% mclapply(function(signature){
    print(signature$name)
    upTags = signature$upTags 
    
    downTags = signature$downTags
    
    out = connectivityMapEnrichment(upTags,downTags, rankMatrix,instances,d = 100000)
    return(out$chemScores)
},mc.cores=16)
MSigDB_enrichFast = MSigDB_enrich
devtools::use_data(MSigDB_enrichFast,overwrite = TRUE)
