library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(glue)
library(ConnectivityMap)
library(rlang)

dir.create('analysis/02.L1000Analysis/results')

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


results = c('analysis/01.CmapEnrichment/results/enrichment/',
            'analysis/02.L1000Analysis/L1000Results/chemScores/',
            'analysis/02.L1000Analysis/pavResults/chemScores/',
            'analysis/02.L1000Analysis/fwdResults/chemScores/'
            #'analysis/01.L1000Analysis/fwdWebRun/similar/',
            #'analysis/01.L1000Analysis/cdsWebRun/'
            )
resultFiles = results %>% lapply(list.files,full.names=TRUE)
names(resultFiles)= c('CMAP','L1000','Pavlab','FWDdata')#,'FWDweb','CDSweb')
nameCol = c('...1','rowname','rowname','rowname','pert_desc','pert_desc')
names(nameCol) = names(resultFiles)
# nameScore = c('sign(enrichment)*(1-p)','sign(enrichment)*p','sign(enrichment)*(1-p)','sign(enrichment)*p','scores','score')
nameScore = c('enrichment','enrichment','enrichment','enrichment','scores','score')

# nameScore = c('p','p','p','p','scores','score')

reliable = c('reliable','reliable','reliable','reliable',NA,NA)
names(reliable) = names(resultFiles)
names(nameScore) = names(resultFiles)
readFun = c(readr::read_tsv,read.table,read.table,read.table,read.table,read.table)

allResults = 1:length(resultFiles) %>% lapply(function(i){
    results = resultFiles[[i]] %>% lapply(readFun[[i]]) %>% lapply(tibble::rownames_to_column)
    names(results)= basename(resultFiles[[i]]) %>% gsub('_enrichment.tsv','',.) %>%
        ogbox::replaceElement(dictionary = groupShorthands) %$% newVector
    results = results[sort(names(results))]
    
    return(results)
})
names(allResults) = names(resultFiles)

corCombinations = combn(names(allResults),2)

corMatrix = matrix(0,nrow = length(allResults),ncol =length(allResults))
colnames(corMatrix) = names(allResults)
rownames(corMatrix) = names(allResults)
groupCors = rep(list(corMatrix),length(groupShorthands))
groupCounts = rep(list(corMatrix),length(groupShorthands))
names(groupCors) = groupShorthands
names(groupCounts) = groupShorthands


reliableFilter = FALSE
for(x in groupShorthands){
    diag(groupCors[[x]]) = NA
    diag(groupCounts[[x]])= NA
    for(i in 1:ncol(corCombinations)){
        group1Name = corCombinations[1,i]
        group2Name =  corCombinations[2,i]
        scoreTable1 = allResults[[group1Name]][[x]]
        if(!is.na(reliable[[group1Name]]) & reliableFilter){
            scoreTable1 = scoreTable1[scoreTable1[[reliable[[group1Name]]]],]
        }
        
        scoreTable2 = allResults[[group2Name]][[x]]
        if(!is.na(reliable[[group2Name]]) & reliableFilter){
            scoreTable2 = scoreTable2[scoreTable2[[reliable[[group2Name]]]],]
        }
        
        commonChems = intersect(scoreTable1[[nameCol[group1Name]]],scoreTable2[[nameCol[group2Name]]])
        
        subTable1 = scoreTable1[match(tolower(commonChems),tolower(scoreTable1[[nameCol[group1Name]]])),]
        subTable2 = scoreTable2[match(tolower(commonChems),tolower(scoreTable2[[nameCol[group2Name]]])),]
        scores1 = ogbox::teval(glue("with(subTable1,{nameScore[group1Name]})"))
        scores2 = ogbox::teval(glue("with(subTable2,{nameScore[group2Name]})"))
        
        groupCors[[x]][group1Name,group2Name] = cor(scores1,scores2,method = 'spearman')
        groupCors[[x]][group2Name,group1Name] = cor(scores1,scores2,method = 'spearman')
        groupCounts[[x]][group2Name,group1Name] = length(scores1)
        groupCounts[[x]][group1Name,group2Name] = length(scores1)
        
        
    }
}

plots = list
names(groupCors)[1:3] %>% lapply(function(name){
    x = groupCors[[name]][1:4,1:4] %>% reshape2::melt()
    y = groupCounts[[name]][1:4,1:4] %>% reshape2::melt()
    
    x = data.frame(Var1 = x$Var1, Var2 = x$Var2, `Spearman's ρ` = x$value, text = paste0(round(x$value,2),'\n(',y$value,')'),check.names=FALSE)
    
    x %>% ggplot(aes(x= Var1,y  = Var2,fill =`Spearman's ρ`, label = text)) + 
        geom_tile() +  scale_fill_continuous(low = 'white',high = '#46A948',na.value = 'black') + 
        geom_text() + ggtitle(name) + xlab('') + ylab('') -> p
   return(p)
}) ->plots

to_save = plots[[1]] + plots[[2]] + plots[[3]]


ggsave('analysis/02.L1000Analysis/results/cor_plot.png',to_save,width = 16, height = 6)