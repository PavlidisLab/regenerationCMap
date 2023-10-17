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




dir.create('analysis/01.CmapEnrichment/results/tags', showWarnings = FALSE,recursive = TRUE)
dir.create('analysis/01.CmapEnrichment/results/enrichment/',recursive= TRUE,showWarnings = FALSE)
dir.create('analysis/01.CmapEnrichment/results/instanceScores/',recursive= TRUE,showWarnings = FALSE)

results = groups %>% lapply(function(group){
    print(group)

    
    filter_criteriaUp = lazyeval::interp(~FC > 0 & Pval < FDRLimit, 
                                         FC = as.name(glue('logFC_{group}')),
                                         Pval =  as.name(glue('FDR_{group}')))
    
    
    upTags = dif_exp_data %>%
        dplyr::filter_(filter_criteriaUp) %>% 
        dplyr::arrange_(.dots = c(glue('desc(logFC_{group})'))) %>% # this line won't be necessary later on
        dplyr::select(gene) %>% 
        unlist %>%
        mouse2human %>% {.$humanGene} %>% unique %>% 
        ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96') 
    
    cat(upTags$Probe,file = glue('analysis/01.CmapEnrichment/results/tags/{group}_upTags.grp'),sep='\n')
    
    
    filter_criteriaDown = lazyeval::interp(~FC < 0 & Pval < FDRLimit, 
                                           FC = as.name(glue('logFC_{group}')),
                                           Pval =  as.name(glue('FDR_{group}')))
    
    downTags = dif_exp_data %>%
        dplyr::filter_(filter_criteriaDown) %>% 
        dplyr::arrange_(.dots = c(glue('logFC_{group}'))) %>% # this line won't be necessary later on
        dplyr::select(gene) %>% 
        unlist %>%
        mouse2human %>% {.$humanGene} %>% unique %>% 
        ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96')
    
    cat(downTags$Probe,file = glue('analysis/01.CmapEnrichment/results/tags/{group}_downTags.grp'),sep='\n')
    
    if(nrow(upTags)==0 & nrow(downTags)==0){
        return(NULL)
    }
    
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
    
    
    write.table(out$chemScores,
                file = glue('analysis/01.CmapEnrichment/results/enrichment/{group}_enrichment.tsv'),
                sep='\t',quote = FALSE,col.names=NA)
    
    write.table(out$instanceScores,
                file = glue('analysis/01.CmapEnrichment/results/instanceScores/{group}_instanceScores.tsv'),
                sep='\t',quote = FALSE,col.names=NA)
    
    return(out)
})

names(results) = groups

scores = results %>% purrr::map('chemScores') %>% purrr::map('enrichment') %>% as.df %>% as.matrix %>% cor(method = 'spearman') %>%
    {diag(.) = NA;.}

colnames(scores) %<>% replaceElement(dictionary =  groupShorthands) %$% newVector
rn(scores) = cn(scores)


png(glue('analysis/01.CmapEnrichment/results/heatmap.png'),width = 1200,height = 800)
scores %>%
    heatmap.2(trace='none',col = viridis::viridis(20),margins= c(14,14),symbreaks = FALSE)
dev.off()

combined = do.call(cbind, purrr::map(results,'chemScores'))

write.table(combined,
            file = glue::glue('analysis/01.CmapEnrichment/results/combined.tsv'),
            sep='\t',quote = FALSE,col.names=NA)

