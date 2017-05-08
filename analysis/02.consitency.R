library(ogbox)
library(purrr)
library(gplots)
library(glue)
datasets = list.files('analysis/results/',recursive = FALSE)

results = datasets %>% lapply(function(x){
    files = list.files(glue('analysis/results/{x}'),full.names = TRUE,pattern = 'enrichment')
    groups =  list.files(glue('analysis/results/{x}'),pattern = 'enrichment')
    files %<>% lapply(ogbox::read.design)
    names(files) = groups
    return(files)
}) 

names(results) = datasets


dir.create('analysis/consistency')
names(results[[1]]) %>% lapply(function(x){
    print(x)
    groupResults = purrr::map(results,x)
    enrichment = purrr::map(groupResults,'enrichment') %>% {.[sapply(.,length)>0]} %>% as.df
    pVals = purrr::map(groupResults,'p') %>% {.[sapply(.,length)>0]} %>% as.df
    # keep = pVals %>% apply(1,function(y){
    #     any(y<0.1)
    # })
    # enrichment = enrichment[keep,]
    
    png(glue('analysis/consistency/{x}.png'),width = 1200,height = 800)
    enrichment %>%as.matrix %>%cor(method = 'spearman') %>%  heatmap.2(trace='none',col = viridis::viridis(20),margins= c(14,14),main = x)
    dev.off()
})
