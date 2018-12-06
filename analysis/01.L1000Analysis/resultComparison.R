library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(glue)
library(ConnectivityMap)
library(rlang)

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


results = c('analysis/results/enrichment/genesEdgerNoOutlier/',
            'analysis/01.L1000Analysis/L1000Results/chemScores/',
            'analysis/01.L1000Analysis/natResults/chemScores/',
            'analysis/01.L1000Analysis/fwdResults/chemScores/',
            'analysis/01.L1000Analysis/fwdWebRun/similar/',
            'analysis/01.L1000Analysis/cdsWebRun/')
resultFiles = results %>% lapply(list.files,full.names=TRUE)
names(resultFiles)= c('CMAP','L1000','Pavlab','FWDdata','FWDweb','CDSweb')
nameCol = c('X1','rowname','rowname','rowname','pert_desc','pert_desc')
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
    
    
chems = tolower(c('DL-thiorphan','milrinone','triflusal','thiorphan','adiphenine'))

names(allResults) %>% lapply(function(resultName){
    names(allResults[[resultName]])[4:6] %>% lapply(function(groupName){
        frame = allResults[[resultName]][[groupName]]
        as_tibble(frame[tolower(frame[[nameCol[resultName]]]) %in% chems,])
    }) ->out
    
    names(out) = names(allResults[[resultName]][4:6]) 
    return(out)
}) -> usedChems

names(usedChems) = names(allResults)


instances = list(L1000 =  readRDS('analysis/00.cmapRanks/instances.rds'),
     Pavlab = readRDS('analysis/00.cmapRanks/NatInstances.rds'),
     FWDData =  readRDS('analysis/00.cmapRanks/FWDinstances.rds'))

# FWD seems to have changed some inames. not sure why but they are recognizable by their pert_ids

validChems = instances$L1000  %>% filter(pert_type %in% c('trt_cp','trt_lig')) %>% select(pert_iname,pert_id) %>% unlist %>% unique
validChems = c(validChems,instances$FWDData %>% filter(pert_id %in% validChems) %$% pert_iname) %>% unique

hitlists = names(allResults)[1:4] %>%  lapply(function(resultName){
    names(allResults[[resultName]])[4:6] %>% lapply(function(groupName){
        nm = unname(nameCol[resultName])
        nmcol = rlang::sym(nm)
        frame = allResults[[resultName]][[groupName]]
        # frame %>% filter(enrichment > sort(frame$enrichment,decreasing = TRUE)[100] & reliable)
        frame %<>% filter(enrichment > 0 & reliable & enrichment>0.5 & p<0.1)
        if(resultName!='CMAP'){
            frame %<>% filter(!!nmcol %in% validChems) #  %>% filter(!grepl('^(BRD)',!!nmcol))
        }
        return(frame)
    }) %>% purrr::map(nameCol[resultName]) %>% ogbox::intersectList()
})

names(hitlists) = c('cmap','L1000','nat','fwd')


instL1000 = readRDS('analysis/00.cmapRanks/instances.rds')
data('instances')
instPavlab = readRDS('analysis/00.cmapRanks/NatInstances.rds')
instFWD =  readRDS('analysis/00.cmapRanks/FWDinstances.rds')


cmapSubset =instances %>% tibble::rownames_to_column() %>% filter(tolower(cmap_name) %in% chems) %>% arrange(`cmap_name`)
L1000Subset = instL1000 %>% filter(tolower(pert_iname) %in% chems & cell_id %in% unique(cmapSubset$cell2)) %>% arrange(`pert_iname`)
pavlabSubset = instPavlab %>% filter(tolower(chem) %in% chems & cellLine %in% unique(cmapSubset$cell2)) %>% arrange(`chem`)
fwdSubset = instFWD %>% tibble::rownames_to_column() %>% filter(tolower(pert_iname) %in% chems & cell_id %in% unique(cmapSubset$cell2)) %>% arrange(`pert_iname`)

data("rankMatrix")
cmapRankSubset = rankMatrix[,cmapSubset$rowname]
gpl96 = gemmaAPI::getAnnotation('GPL96')
rankMedian = cmapRankSubset %>%unlist %>%median
# get the probesets at the highest extremes in general
extremeProbes = abs(cmapRankSubset- rankMedian) %>% apply(1,sum) %>% sort(decreasing = TRUE)
useProbes = gemmaAPI::annotationGeneMatch(names(extremeProbes),gpl96) %>% duplicated %>% not %>% {extremeProbes[.]} %>% names
useProbes = gemmaAPI::annotationGeneMatch(useProbes,gpl96) %>% {.[!is.na(.)]} %>% names
cmapRankSubset = cmapRankSubset[useProbes,]
rownames(cmapRankSubset) = gemmaAPI::annotationGeneMatch(rownames(cmapRankSubset),gpl96)


rankMatrixL1000 = readRDS('analysis/00.cmapRanks/rankMatrix.rds')
L1000geneAnnots = readRDS('analysis/00.cmapRanks/L1000geneAnnots.rds')
L1000RankSubset = rankMatrixL1000[,L1000Subset$sig_id]
# genes = L1000geneAnnots[match(rownames(L1000RankSubset),L1000geneAnnots$pr_gene_id),]$pr_gene_symbol
genes = L1000geneAnnots$pr_gene_symbol
rownames(L1000RankSubset) = genes
L1000RankSubset = L1000RankSubset[L1000geneAnnots$pr_is_lm ==1,]

commonGenes = ogbox::intersectMult(rownames(L1000RankSubset),rownames(cmapRankSubset))

cmapRankSubsetCommon = cmapRankSubset[commonGenes,]
L1000RankSubsetCommon = L1000RankSubset[commonGenes,]

colnames(cmapRankSubsetCommon) = paste('cmap',cmapSubset$cmap_name,cmapSubset$cell2,cmapSubset$duration..h.,'h')
colnames(L1000RankSubsetCommon) = paste('l1000',L1000Subset$pert_iname,L1000Subset$cell_id,L1000Subset$pert_itime)

cbind(cmapRankSubsetCommon,L1000RankSubsetCommon) %>% cor(method = 'spearman') %>% {diag(.) = NA;.}%>% pheatmap::pheatmap()

cor(L1000RankSubsetCommon,cmapRankSubsetCommon,method = 'spearman') %>% pheatmap::pheatmap(cluster_rows = F,cluster_cols = F)

L1000RankSubset[,]

allResults$CMAP$`regen 1 week` %>% filter(tolower(X1) %in% chems)
allResults$CMAP$`regen 2 week` %>% filter(X1 %in% chems)
allResults$CMAP$`regen 3 days` %>% filter(X1 %in% chems)


allResults$L1000$`regen 1 week` %>% filter(rowname %in% chems)
allResults$L1000$`regen 2 week` %>% filter(rowname %in% chems)
allResults$L1000$`regen 3 days` %>% filter(rowname %in% chems)



allResults$Pavlab$`regen 1 week` %>% filter(rowname %in% chems)
allResults$Pavlab$`regen 2 week` %>% filter(rowname %in% chems)
allResults$Pavlab$`regen 3 days` %>% filter(rowname %in% chems)



allResults$FWDdata$`regen 1 week` %>% filter(rowname %in% chems)
allResults$FWDdata$`regen 2 week` %>% filter(rowname %in% chems)
allResults$FWDdata$`regen 3 days` %>% filter(rowname %in% chems)




allResults$FWDweb$`regen 1 week` %>% filter(rowname %in% tolower(allResults$FWDweb$`regen 1 week`$pert_desc))
fwdScores$`regen 2 week` %>% filter(rowname %in% fwdWeb$`regen 2 week`$pert_desc)
fwdScores$`regen 3 days` %>% filter(rowname %in% fwdWeb$`regen 3 days`$pert_desc)

commonChems = ogbox::intersectMult(oldResults$`regen 2 week`$X1,newResults$`regen 2 week`$rowname,natScores$`regen 2 week`$rowname,fwdScores$`regen 2 week`$rowname)
commonChemsOldNew = ogbox::intersectMult(oldResults$`regen 2 week`$X1,newResults$`regen 2 week`$rowname)
commonChemsoldNat = ogbox::intersectMult(oldResults$`regen 2 week`$X1,natScores$`regen 2 week`$rowname)
commonChemsoldFWD = ogbox::intersectMult(oldResults$`regen 2 week`$X1,fwdScores$`regen 2 week`$rowname)


oldCommon = oldResults$`regen 2 week`[match(commonChems,oldResults$`regen 2 week`$X1),]$enrichment
newCommon = newResults$`regen 2 week`[match(commonChems,newResults$`regen 2 week`$rowname),]$enrichment
natCommon = natScores$`regen 2 week`[match(commonChems,natScores$`regen 2 week`$rowname),]$enrichment
fwdCommmon = fwdScores$`regen 2 week`[match(commonChems,fwdScores$`regen 2 week`$rowname),]$enrichment

data.frame(-oldCommon,newCommon,natCommon,fwdCommmon) %>% cor %>% reshape2::melt() %>% ggplot(aes(x= Var1,y  = Var2,fill =value, label = value)) + 
    geom_tile() +  scale_fill_continuous(low = 'white',high = '#46A948',na.value = 'white') + 
    geom_text(aes(label = value %>% round(2) %>% format(nsmall=2)))






oldCommon = oldResults$`regen 1 week`[match(commonChems,oldResults$`regen 2 week`$X1),]$enrichment
newCommon = newResults$`regen 1 week`[match(commonChems,newResults$`regen 2 week`$rowname),]$enrichment
natCommon = natScores$`regen 1 week`[match(commonChems,natScores$`regen 2 week`$rowname),]$enrichment
fwdCommmon = fwdScores$`regen 1 week`[match(commonChems,fwdScores$`regen 2 week`$rowname),]$enrichment

data.frame(oldCommon,newCommon,natCommon,fwdCommmon) %>% cor %>% pheatmap::pheatmap()


oldCommon = oldResults$`regen 3 days`[match(commonChems,oldResults$`regen 2 week`$X1),]$enrichment
newCommon = newResults$`regen 3 days`[match(commonChems,newResults$`regen 2 week`$rowname),]$enrichment
natCommon = natScores$`regen 3 days`[match(commonChems,natScores$`regen 2 week`$rowname),]$enrichment
fwdCommmon = fwdScores$`regen 3 days`[match(commonChems,fwdScores$`regen 2 week`$rowname),]$enrichment

data.frame(oldCommon,newCommon,natCommon,fwdCommmon) %>% cor %>% pheatmap::pheatmap()



goodOldResults = oldResults$`regen 2 week` %>% arrange(desc(enrichment)) %>% filter(reliable) %>% filter(X1 %in% commonChems)

newResults$`regen 1 week`[match(goodOldResults$X1[1:10],newResults$`regen 2 week`$rowname),]
natScores$`regen 1 week`[match(goodOldResults$X1[1:10],natScores$`regen 2 week`$rowname),]




higsScoreCommon = oldResults$`regen 2 week`[match(commonChemsoldNat,oldResults$`regen 2 week`$X1),] %>% arrange(desc(enrichment)) %>% 
    filter(reliable) %>% {.$X1[1:100]}


oldCommon = oldResults$`regen 2 week`[match(higsScoreCommon,oldResults$`regen 2 week`$X1),]$enrichment
newCommon = newResults$`regen 2 week`[match(higsScoreCommon,newResults$`regen 2 week`$rowname),]$enrichment
natCommon = natScores$`regen 2 week`[match(higsScoreCommon,natScores$`regen 2 week`$rowname),]$enrichment
data.frame(oldCommon,newCommon,natCommon) %>% cor %>% reshape2::melt() %>% ggplot(aes(x= Var1,y  = Var2,fill =value, label = value)) + 
    geom_tile() +  scale_fill_continuous(low = 'white',high = '#46A948',na.value = 'white') + 
    geom_text(aes(label = value %>% round(2) %>% format(nsmall=2)))

