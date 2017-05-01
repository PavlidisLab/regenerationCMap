library(glue)
library(magrittr)
library(dplyr)
library(ogbox)
library(lazyeval)
library(homologene)
library(purrr)
library(parallel)
library(ConnectivityMap)
data("rankMatrix")
data("instances")
devtools::load_all()

FDRLimit = .01
# temporary for testing purposes
geneCountLimit = 250


groups = c("E12_1_week_IP_vs_naive_adult_3_IP",
           "E12_2_week_IP_vs_naive_adult_3_IP",
           "E12_3_day_IP_vs_naive_adult_3_IP",
           "naive_1_week_IP_vs_naive_adult_3_IP",
           "naive_2_weeks_IP_vs_naive_adult_3_IP",
           "naive_3_days_IP_vs_naive_adult_3_IP")

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

dir.create('data-raw/tags', showWarnings = FALSE)
datasets %>% sapply(function(dataset){
    print(dataset)
    dir.create(glue('analysis/results/{dataset}'),recursive= TRUE,showWarnings = FALSE)
    dir.create(glue('analysis/probes/{dataset}'),recursive= TRUE, showWarnings = FALSE)
    
    groups %>% sapply(function(group){
        print(group)
        filter_criteriaUp = lazyeval::interp(~FC > 0 & Pval < FDRLimit, 
                                             FC = as.name(glue('logFC_{group}')),
                                             Pval =  as.name(glue('FDR_pVal_{group}')))
        
        upTags = ogbox::teval(dataset) %>%
            dplyr::filter_(filter_criteriaUp) %>% 
            dplyr::arrange_(.dots = c(glue('desc(logFC_{group})'))) %>% # this line won't be necessary later on
            dplyr::select(gene) %>% 
            unlist %>%
            mouse2human %>% {.$humanGene} %>% unique %>% 
            ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96')  %>% 
             {.[1:geneCountLimit,]} # this line is temporary and should be removed later
        
        cat(upTags$Probe,file = glue('analysis/probes/{dataset}/{group}_upTags.grp'),sep='\n')
        
        
        filter_criteriaDown = lazyeval::interp(~FC < 0 & Pval < FDRLimit, 
                                               FC = as.name(glue('logFC_{group}')),
                                               Pval =  as.name(glue('FDR_pVal_{group}')))
        
        downTags = ogbox::teval(dataset) %>%
            dplyr::filter_(filter_criteriaDown) %>% 
            dplyr::arrange_(.dots = c(glue('logFC_{group}'))) %>% # this line won't be necessary later on
            dplyr::select(gene) %>% 
            unlist %>%
            mouse2human %>% {.$humanGene} %>% unique %>% 
            ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96')  %>% 
            {.[1:min(geneCountLimit,nrow(.)),]} # this line is temporary and should be removed later
        
        cat(downTags$Probe,file = glue('analysis/probes/{dataset}/{group}_downTags.grp'),sep='\n')
        
        n = rankMatrix %>% nrow # "Let n be the number of probesets in the featureset (22,283)"
        
        
        scores = mclapply(colnames(rankMatrix),function(instance){
            print(instance)
            orderedProbes = rownames(rankMatrix)[rankMatrix[instance] %>% order] 
            
            Vup = match(upTags$Probe,orderedProbes) %>% sort
            Vdown = match(downTags$Probe,orderedProbes) %>% sort
            return(scoreCalc(Vup,Vdown,n))
            
        },
        mc.cores=5)
        names(scores) = colnames(rankMatrix)
        # "Order all n probe sets by the extent of their differential expression for the current instance i"
        # the data itself is in ranks so it is already ordered but lets do what they say
        
        # "p to be max( si) and q to be min( si) across all instances in the collection c"
        p = max(map_dbl(scores,'score'))
        q = min(map_dbl(scores,'score'))
        
        scores = scores %>% lapply(function(x){
            s = x[['score']]
            ConScore = (s>0)*(s/p) + (s<0)*(-s/q)
            x = c(x, ConScore = ConScore)
        })
        scores %<>% as.df %>% t %>% as.df
        scores$instance = colnames(rankMatrix)
        
        # "The Kolmogorov-Smirnov statistic is computed for the set of t instances
        # in the list of all n instances in a result ordered in descending order of
        # connectivity score and up score (see how connectivity score is
        # calculated), giving an enrichment score ks0."
        scores = scores %>% arrange(desc(ConScore),desc(kUp))
        
        chems = instances$cmap_name %>% unique
        
        d = 10000
        
        set.seed(1)
        confidence = chems %>% sapply(function(chem){
            print(chem)
            chemInstances = rownames(instances)[instances$cmap_name %in% chem]
            V = match(chemInstances,scores$instance) %>% sort
            ks0 = compiler::cmpfun(ksCalc)(V,nrow(instances))
            hede = compiler::cmpfun(ksCalc)
            
            permute = 1:d %>% mclapply(function(i){

                Vrandom = sample(x=1:length(scores$instance),
                                 size = length(chemInstances),replace = FALSE) %>% sort
                ksm = ksCalc(Vrandom,nrow(instances))
            },mc.cores=5) %>% unlist
            proc.time() - ptm
            q = sum(abs(permute) >= abs(ks0))
            p = q/d
            
            return(c(enrichment = ks0, p = p))
        }) %>% t
        
        write.table(confidence,
                    file = glue('analysis/results/{dataset}/{group}_enrichment.tsv'),
                    sep='\t',quote = FALSE,col.names=NA)
        write.table(scores,
                    file = glue('analysis/results/{dataset}/{group}_instanceScores.tsv'),
                    sep='\t',quote = FALSE,col.names=NA)
    },simplify= FALSE)
},simplify = FALSE)