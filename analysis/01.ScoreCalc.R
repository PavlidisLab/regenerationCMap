library(glue)
library(magrittr)
library(dplyr)
library(ogbox)
library(lazyeval)
library(homologene)
library(parallel)
library(ConnectivityMap)
library(memoise)
library(VennDiagram)
library(gplots)
library(pheatmap)
library(purrr)
data("rankMatrix")
data("instances")
library(UpSetR)
devtools::load_all()

FDRLimit = 0.05
# temporary for testing purposes
# geneCountLimit = 250


groups = c("E12_1_week_IP_vs_naive_adult_3_IP",
           "E12_2_week_IP_vs_naive_adult_3_IP",
           "E12_3_day_IP_vs_naive_adult_3_IP",
           "naive_1_week_IP_vs_naive_adult_3_IP",
           "naive_2_weeks_IP_vs_naive_adult_3_IP",
           "naive_3_days_IP_vs_naive_adult_3_IP")

groupShorthands = c("E12_1_week_IP_vs_naive_adult_3_IP" = 'regen 1 week',
                    "E12_2_week_IP_vs_naive_adult_3_IP" = 'regen 2 week',
                    "E12_3_day_IP_vs_naive_adult_3_IP" = 'regen 3 days',
                    "naive_1_week_IP_vs_naive_adult_3_IP" = 'naive 1 week',
                    "naive_2_weeks_IP_vs_naive_adult_3_IP" = 'naive 2 weeks',
                    "naive_3_days_IP_vs_naive_adult_3_IP" = 'naive 3 days')

# datasets = list(genesVoomLimma= genesVoomLimma,
#                 genesLimma = genesLimma,
#                 genesEdger = genesEdger,
#                 genesVoomLimmaRUVrk1= genesVoomLimmaRUVrk1,
#                 genesLimmaRUVrk1 = genesLimmaRUVrk1,
#                 genesEdgerRUVrk1 = genesEdgerRUVrk1,
#                 genesVoomLimmaNoOutlier= genesVoomLimmaNoOutlier,
#                 genesLimmaNoOutlier = genesLimmaNoOutlier,
#                 genesEdgerNoOutlier = genesEdgerNoOutlier)
datasets = c('genesVoomLimma',
             'genesVoomLimmaRUVrk1',
             'genesVoomLimmaNoOutlier',
             'genesLimma' ,
             'genesEdger' ,
             'genesLimmaRUVrk1',
             'genesEdgerRUVrk1',
             'genesLimmaNoOutlier' ,
             'genesEdgerNoOutlier' )

datasetShorthands = c('genesVoomLimma' = "VoLimma",
                      'genesVoomLimmaRUVrk1' = "VoLimmaaRuv",
                      'genesVoomLimmaNoOutlier' = "VoLimmanoOut",
                      'genesLimma' = "Limma",
                      'genesEdger' = "Edge",
                      'genesLimmaRUVrk1' = "LimmaRUV",
                      'genesEdgerRUVrk1' = 'EdgeRUV',
                      'genesLimmaNoOutlier' = 'LimmanoOut',
                      'genesEdgerNoOutlier' = 'EdgenoOut' )

# groups = datasets %>% lapply(function(x){
#     colnames(teval(x))[grepl('FDR',colnames(teval(x)))]
# }
# ) %>% unlist %>% unique


# dataset  ='genesEdgerNoOutlier'
dataset = 'genesVoomLimma'
#group = 'naive_3_days_IP_vs_naive_adult_3_IP'
#group = 'naive_2_weeks_IP_vs_naive_adult_3_IP'
# group = 'naive_2_weeks_IP_vs_naive_adult_3_IP'
group = 'E12_1_week_IP_vs_naive_adult_3_IP'
dir.create('data-raw/tags', showWarnings = FALSE)

dir.create(glue('analysis/results/enrichmentMonolith/'),recursive= TRUE,showWarnings = FALSE)
dir.create(glue('analysis/results/heatmaps/'),recursive= TRUE,showWarnings = FALSE)

allResults = datasets %>% mclapply(function(dataset){
    print(dataset)
    dir.create(glue('analysis/results/enrichment/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/results/perturbagenHitlist/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/results/perturbagenSpecificHitlist/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/results/perturbagenHitlist01/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/results/perturbagenHitlistPos/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/results/perturbagenHitlistPos01/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    
    
    dir.create(glue('analysis/results/instanceScores/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    
    dir.create(glue('analysis/probes/{dataset}'),recursive= TRUE, showWarnings = FALSE)
    
    results = groups %>% lapply(function(group){
        print(group)
        if(any(grepl('FDR_pVal',colnames(ogbox::teval(dataset))))){
            pVal = 'pVal_'
        } else{
            pVal = ''
        }
        
        filter_criteriaUp = lazyeval::interp(~FC > 0 & Pval < FDRLimit, 
                                             FC = as.name(glue('logFC_{group}')),
                                             Pval =  as.name(glue('FDR_{pVal}{group}')))
        upTags = ogbox::teval(dataset) %>%
            dplyr::filter_(filter_criteriaUp) %>% 
            dplyr::arrange_(.dots = c(glue('desc(logFC_{group})'))) %>% # this line won't be necessary later on
            dplyr::select(gene) %>% 
            unlist %>%
            mouse2human %>% {.$humanGene} %>% unique %>% 
            ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96')  #%>% 
        #  {.[1:geneCountLimit,]} # this line is temporary and should be removed later
        
        cat(upTags$Probe,file = glue('analysis/probes/{dataset}/{group}_upTags.grp'),sep='\n')
        
        
        filter_criteriaDown = lazyeval::interp(~FC < 0 & Pval < FDRLimit, 
                                               FC = as.name(glue('logFC_{group}')),
                                               Pval =  as.name(glue('FDR_{pVal}{group}')))
        
        downTags = ogbox::teval(dataset) %>%
            dplyr::filter_(filter_criteriaDown) %>% 
            dplyr::arrange_(.dots = c(glue('logFC_{group}'))) %>% # this line won't be necessary later on
            dplyr::select(gene) %>% 
            unlist %>%
            mouse2human %>% {.$humanGene} %>% unique %>% 
            ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96')  #%>% 
        # {.[1:min(geneCountLimit,nrow(.)),]} # this line is temporary and should be removed later
        
        cat(downTags$Probe,file = glue('analysis/probes/{dataset}/{group}_downTags.grp'),sep='\n')
        
        n = rankMatrix %>% nrow # "Let n be the number of probesets in the featureset (22,283)"
        if(nrow(upTags)==0 & nrow(downTags)==0){
            return(NULL)
        }
        # instance = 'inst_7402'
        
        out = connectivityMapEnrichment(upTags$Probe,downTags$Probe, rankMatrix,instances,d = 100000)
        
        out$chemScores$specificity = 1:nrow(out$chemScores) %>%  sapply(function(i){
            chem = rownames(out$chemScores)[i]
            enrichments = MSigDB_enrich %>% sapply(function(x){
                x[chem,'enrichment']
            })
            sign = sign(out$chemScores[chem,'enrichment'])
            if(sign==-1){
                enrichments = enrichments[enrichments<=0]
            } else if(sign ==1){
                enrichments = enrichments[enrichments>=0]
                
            }
            sum(abs(out$chemScores[chem,'enrichment']) < abs(enrichments))/length(enrichments)
        })
        
        out$chemScores$reliable = out$chemScores$nonNull>0.5 & out$chemScores$instanceCount > 1 
        
        out$hits =  out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.05,] %>% {.[order(.$FDR),]} %>% rownames
        out$hits01 =  out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.1,] %>% {.[order(.$FDR),]} %>% rownames
        
        out$hitsPos = out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.05 & out$chemScores$enrichment>0,] %>% {.[order(.$FDR),]} %>% rownames
        out$hitsPos01 =  out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.1 & out$chemScores$enrichment>0,] %>% {.[order(.$FDR),]} %>% rownames
        
        out$specificHits = out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.05 & out$chemScores$specificity<0.1,] %>% rownames
        
        write.table(out$chemScores,
                    file = glue('analysis/results/enrichment/{dataset}/{group}_enrichment.tsv'),
                    sep='\t',quote = FALSE,col.names=NA)
        write.table(out$instanceScores,
                    file = glue('analysis/results/instanceScores/{dataset}/{group}_instanceScores.tsv'),
                    sep='\t',quote = FALSE,col.names=NA)
        write.table(out$hits,
                    file = glue('analysis/results/perturbagenHitlist/{dataset}/{group}_hitlists.tsv'),
                    sep='\t',quote = FALSE,col.names=FALSE,row.names=FALSE)
        

        write.table(out$hits01,
                    file = glue('analysis/results/perturbagenHitlist01/{dataset}/{group}_hitlists.tsv'),
                    sep='\t',quote = FALSE,col.names=FALSE,row.names=FALSE)
        
        write.table(out$hitsPos,
                    file = glue('analysis/results/perturbagenHitlistPos/{dataset}/{group}_hitlists.tsv'),
                    sep='\t',quote = FALSE,col.names=FALSE,row.names=FALSE)
        write.table(out$hitsPos01,
                    file = glue('analysis/results/perturbagenHitlistPos01/{dataset}/{group}_hitlists.tsv'),
                    sep='\t',quote = FALSE,col.names=FALSE,row.names=FALSE)
        
        return(out)
    }# ,mc.cores= 6
    )
    names(results) = groups
    
    results =  results[(results %>% sapply(class))!='NULL']
    
    scores = results %>% purrr::map('chemScores') %>% purrr::map('enrichment') %>% as.df %>% as.matrix %>% cor(method = 'spearman') %>%
    {diag(.) = NA;.}
    
    colnames(scores) %<>% replaceElement(dictionary =  groupShorthands) %$% newVector
    rn(scores) = cn(scores)
    
    png(glue('analysis/results/heatmaps/{dataset}.png'),width = 1200,height = 800)
    scores %>%
        heatmap.2(trace='none',col = viridis::viridis(20),margins= c(14,14),main = dataset,symbreaks = FALSE)
    dev.off()
    
    
    results = results[(results %>% sapply(class))!='NULL']
    monolith = do.call(cbind, purrr::map(results,'chemScores'))
    
    write.table(monolith,
                file = glue::glue('analysis/results/enrichmentMonolith/{dataset}_monolith.tsv'),
                sep='\t',quote = FALSE,col.names=NA)
    return(results)
    
    
    
},mc.cores = 9)

names(allResults) = datasets
saveRDS(allResults,'analysis/allResults.rds')


files = list.files('analysis/results/enrichment/', recursive = TRUE,full.names = TRUE)
names = list.files('analysis/results/enrichment/', recursive = TRUE)

names %<>% sapply(function(x){
    x %>% gsub(pattern = '_enrichment.tsv',replacement = '', x = .) %>% strsplit('/') %>% {.[[1]]} %>% 
        replaceElement(groupShorthands) %$% newVector %>% replaceElement(datasetShorthands) %$% newVector %>%
        paste(collapse =' ')
})

allEnrich = files %>% sapply(function(x){
    data.table::fread(x,data.table = FALSE)$enrichment
}) 

colnames(allEnrich) = names

# monolith heatmap -----------
png(glue('analysis/results/heatmaps/monolithHeatmap.png'),width = 1400,height = 1200)

annotCol = data.frame(Normalization = stringr::str_extract(names,'^.*?(?=\\s)'),
                      group = stringr::str_extract(names,ogbox::regexMerge(groupShorthands,exact= TRUE)),stringsAsFactors = FALSE)
rownames(annotCol) = colnames(allEnrich)

annotColors = list(Normalization = annotCol$Normalization %>% replaceElement(c('Edge' = 'cyan',
                                                                               'EdgeRUV' = 'blue1',
                                                                               'EdgenoOut' = 'cadetblue4',
                                                                               'Limma' = 'chartreuse',
                                                                               'LimmaRUV' = 'chartreuse4',
                                                                               'LimmanoOut' = 'aquamarine',
                                                                               'VoLimma' = 'brown4',
                                                                               'VoLimmanoOut' = 'coral',
                                                                               'VoLimmaaRuv' = 'brown1')) %$% dictionary,
                   group = annotCol$group %>% toColor %$% palette)

allEnrich %>% cor(method = 'spearman') %>% {diag(.) = NA;.} %>% 
    pheatmap(annotation_col =annotCol ,fontsize = 15,border_color = NA,
             annotation_colors = annotColors,color = viridis::viridis(20))
#heatmap.2(trace='none',col = viridis::viridis(20),margins= c(14,14),symbreaks = FALSE,cexCol = 1.1,cexRow = 1.1)
dev.off()

filesEnrichment =   list.files('analysis/results/enrichment/', recursive = FALSE,full.names = TRUE)
datasetNames = basename(filesEnrichment)


# this gets the definitive hitlist
datasetEnrichments = filesEnrichment %>% 
    sapply(function(x){
        datasetFiles = list.files(x,full.names = TRUE)
        out = datasetFiles %>% lapply(function(y){
            read.design(y) %>% arrange(desc(reliable),FDR)
        })
        names(out) = basename(datasetFiles) %>% stringr::str_replace('_enrichment.tsv','') %>% replaceElement(groupShorthands) %$% newVector
        
        ineligableRegenPos = out[grepl(pattern = 'naive',names(out))] %>% sapply(function(y){
            y %>% filter(reliable == TRUE & enrichment > 0 & FDR<0.05) %$% X
        })
        
        ineligableNaivePos = out[grepl(pattern = 'regen',names(out))] %>% sapply(function(y){
            y %>% filter(reliable == TRUE & enrichment > 0 & FDR<0.05) %$% X
        })
        
        ineligableRegenNeg = out[grepl(pattern = 'regen',names(out))] %>% sapply(function(y){
            y %>% filter(reliable == TRUE & enrichment < 0 & FDR<0.05) %$% X
        })
        
        ineligableNaiveNeg = out[grepl(pattern = 'naive',names(out))] %>% sapply(function(y){
            y %>% filter(reliable == TRUE & enrichment < 0 & FDR<0.05) %$% X
        })
        
        topRegenMarkersPos = out[grepl(pattern = 'regen',names(out))] %>% sapply(function(y){
            y %>% filter(!X %in% unlist(ineligableRegenPos) & reliable==TRUE & enrichment > 0 & FDR <0.05) %$% X
        },simplify = FALSE)
        
        topRegenMarkersNeg = out[grepl(pattern = 'regen',names(out))] %>% sapply(function(y){
            y %>% filter(!X %in% unlist(ineligableRegenNeg) & reliable==TRUE & enrichment < 0 & FDR <0.05) %$% X
        },simplify = FALSE)
        
        return(list(topRegenMarkersPos = topRegenMarkersPos,
                    topRegenMarkersNeg = topRegenMarkersNeg))
        
    },simplify=FALSE)
names(datasetEnrichments) = basename(filesEnrichment)

intersectList(datasetEnrichments$genesEdgerNoOutlier$topRegenMarkersPos[c('regen 1 week','regen 2 week')])
intersectList(datasetEnrichments$genesLimmaNoOutlier$topRegenMarkersPos[c('regen 1 week','regen 2 week')])
intersectList(datasetEnrichments$genesVoomLimmaNoOutlier$topRegenMarkersPos[c('regen 1 week','regen 2 week')])
intersectList(datasetEnrichments$genesEdger$topRegenMarkersPos[c('regen 1 week','regen 2 week')])


# upsets ----------------

files = list.files('analysis/results/perturbagenHitlist/', recursive = FALSE,full.names = TRUE)
datasetHitlists = files %>% sapply(function(x){
    datasetFiles = list.files(x,full.names = TRUE)
    out = datasetFiles %>% sapply(function(y){
        readLines(y)
    })
    names(out) = basename(datasetFiles) %>% stringr::str_replace('_hitlists.tsv','') %>% replaceElement(groupShorthands) %$% newVector
    return(out)
})

names(datasetHitlists) = basename(files)

dir.create('analysis/results/upsetWithinGroups')
names(datasetHitlists) %>% sapply(function(x){
    png(glue::glue('analysis/results/upsetWithinGroups/{x}.png'),width = 600,height = 600)
    upset(UpSetR::fromList(datasetHitlists[[x]]),nsets = 6, text.scale = 2)
    dev.off()
})


files = list.files('analysis/results/perturbagenHitlist01//', recursive = FALSE,full.names = TRUE)
datasetHitlists = files %>% sapply(function(x){
    datasetFiles = list.files(x,full.names = TRUE)
    out = datasetFiles %>% sapply(function(y){
        readLines(y)
    })
    names(out) = basename(datasetFiles) %>% stringr::str_replace('_hitlists.tsv','') %>% replaceElement(groupShorthands) %$% newVector
    return(out)
})

names(datasetHitlists) = basename(files)


dir.create('analysis/results/upsetWithinGroups01')
names(datasetHitlists) %>% sapply(function(x){
    png(glue::glue('analysis/results/upsetWithinGroups01/{x}.png'),width = 600,height = 600)
    upset(UpSetR::fromList(datasetHitlists[[x]]),nsets = 6, text.scale = 2)
    dev.off()
})




files = list.files('analysis/results/perturbagenHitlistPos//', recursive = FALSE,full.names = TRUE)
datasetHitlists = files %>% sapply(function(x){
    datasetFiles = list.files(x,full.names = TRUE)
    out = datasetFiles %>% sapply(function(y){
        readLines(y)
    })
    names(out) = basename(datasetFiles) %>% stringr::str_replace('_hitlists.tsv','') %>% replaceElement(groupShorthands) %$% newVector
    return(out)
})

names(datasetHitlists) = basename(files)


dir.create('analysis/results/upsetWithinGroupsPos')
names(datasetHitlists) %>% sapply(function(x){
    png(glue::glue('analysis/results/upsetWithinGroupsPos/{x}.png'),width = 600,height = 600)
    upset(UpSetR::fromList(datasetHitlists[[x]]),nsets = 6, text.scale = 2)
    dev.off()
})



files = list.files('analysis/results/perturbagenHitlistPos01//', recursive = FALSE,full.names = TRUE)
datasetHitlists = files %>% sapply(function(x){
    datasetFiles = list.files(x,full.names = TRUE)
    out = datasetFiles %>% sapply(function(y){
        readLines(y)
    })
    names(out) = basename(datasetFiles) %>% stringr::str_replace('_hitlists.tsv','') %>% replaceElement(groupShorthands) %$% newVector
    return(out)
})

names(datasetHitlists) = basename(files)


dir.create('analysis/results/upsetWithinGroupsPos01')
names(datasetHitlists) %>% sapply(function(x){
    png(glue::glue('analysis/results/upsetWithinGroupsPos01/{x}.png'),width = 600,height = 600)
    upset(UpSetR::fromList(datasetHitlists[[x]]),nsets = 6, text.scale = 2)
    dev.off()
})
# groupContrasts = names(datasetHitlists) %>% sapply(function(x){
#     day3ExclusiveRegen = datasetHitlists[[x]]$`regen 3 days`[!datasetHitlists[[x]]$`regen 3 days` %in% datasetHitlists[[x]]$`naive 3 days`]
#     week1ExclusiveRegen  = datasetHitlists[[x]]$`regen 1 week`[!datasetHitlists[[x]]$`regen 1 week` %in% datasetHitlists[[x]]$`naive 1 week`]
#     
#     
#     regenLists = datasetHitlists[[x]][grepl(pattern = 'regen', names(datasetHitlists[[x]]))] %>% unlist %>% unique
#     naiveList = datasetHitlists[[x]][grepl(pattern = 'naive', names(datasetHitlists[[x]]))] %>% unlist %>% unique
#     commonList = intersect(regenLists,naiveList)
#     exclusiveRegen = regenLists[!regenLists %in% commonList]
#     exclusiveNaive = naiveList[!naiveList %in% commonList]
#     return(list(exclusiveNaive = exclusiveNaive,
#                 exclusiveRegen = exclusiveRegen,
#                 commonList = commonList))
# },simplify = FALSE)
# 
# groupContrasts %>% purrr::map('commonList') %>% fromList %>% upset(nsets = 9)

# hitlistNames = list.files('analysis/results/perturbagenHitlist/', recursive = TRUE)
# 
# hitlistNames %<>% sapply(function(x){
#     x %>% gsub(pattern = '_hitlists.tsv',replacement = '', x = .,perl = TRUE) %>% strsplit('/') %>% {.[[1]]} %>% 
#         replaceElement(groupShorthands) %$% newVector %>% replaceElement(datasetShorthands) %$% newVector %>%
#         paste(collapse =' ')
# })
# 
# hitlists = files %>% sapply(function(x){
#     print(x)
#     readLines(x)
# }) 
# names(hitlists) = hitlistNames
# 
# hitlistNames
# 
# all(names == hitlistNames)
