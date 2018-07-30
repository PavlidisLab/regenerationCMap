print(Sys.getpid())
library(dplyr)
library(cmapQuery)
library(magrittr)
library(glue)
library(homologene)

dir.create('analysis/01.L1000Analysis/fwdWebRun/similar',showWarnings = FALSE, recursive = TRUE)
dir.create('analysis/01.L1000Analysis/fwdWebRun/opposite',showWarnings = FALSE, recursive = TRUE)


print("loading data")
load('data/genesEdgerNoOutlier.rda')


dataset = genesEdgerNoOutlier

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

FDRLimit = 0.05


groups %>% lapply(function(group){
    print(group)
    if(any(grepl('FDR_pVal',colnames(dataset)))){
        pVal = 'pVal_'
    } else{
        pVal = ''
    }
    filter_criteriaUp = lazyeval::interp(~FC > 0 & Pval < FDRLimit, 
                                         FC = as.name(glue('logFC_{group}')),
                                         Pval =  as.name(glue('FDR_{pVal}{group}')))
    
    upGenes = dataset %>%
        dplyr::filter_(filter_criteriaUp) %>% 
        dplyr::arrange_(.dots = c(glue('desc(logFC_{group})'))) %>% # this line won't be necessary later on
        dplyr::select(gene) %>% 
        unlist %>%
        mouse2human %>% {.$humanGene} %>% unique
    
    filter_criteriaDown = lazyeval::interp(~FC < 0 & Pval < FDRLimit, 
                                           FC = as.name(glue('logFC_{group}')),
                                           Pval =  as.name(glue('FDR_{pVal}{group}')))
    
    downGenes = dataset %>%
        dplyr::filter_(filter_criteriaDown) %>% 
        dplyr::arrange_(.dots = c(glue('logFC_{group}'))) %>% # this line won't be necessary later on
        dplyr::select(gene) %>% 
        unlist %>%
        mouse2human %>% {.$humanGene} %>% unique 
    
    results = cmapQuery::L1000FWD(upGenes,downGenes)
    write.table(results$similar,
                file = glue('analysis/01.L1000Analysis/fwdWebRun/similar/{groupShorthands[group]}'))
    write.table(results$opposite,
                file = glue('analysis/01.L1000Analysis/fwdWebRun/opposite/{groupShorthands[group]}'))
})