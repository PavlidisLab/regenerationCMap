# install_github("oganm/MSigDB")
library(MSigDB)
# install_github("oganm/ogbox")
library(ogbox)
# install_github("oganm/homologene")
library(homologene)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ConnectivityMap")
library(ConnectivityMap)

library(dplyr)
library(magrittr)

data(instances)
data(rankMatrix)
# replication trial
devtools::load_all()

# calculation attempt for instance 868 for logFC_E12_3_day_IP_vs_naive_3_days_IP
# web interface results up=.133, down=-.188, score=.856
# tests can be run at https://portals.broadinstitute.org/cmap/index.jsp
# OganM username. Ones tagged 'real'

# getting the top 250 genes from my data as the upregulated list. order should not matter but it is ordered from most DE to least DE
# some irrelevant magic about turning genes into probeset occurs here
# ultimately the input is the probe section of this dataframe
# up genes --------------
posList = hitlist %>%
    arrange(desc(logFC_E12_3_day_IP_vs_naive_3_days_IP)) %>% 
    filter(logFC_E12_3_day_IP_vs_naive_3_days_IP>0) %>% 
    dplyr::select(V1) %>%unlist %>% mouse2human %>% {.$humanGene} %>% unique %>% ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96') %>% {.[1:250,]}

# tests for the web interface
# reverted list
posListRev = posList %>%
    mutate(Probe = rev(Probe),Gene.Symbol = rev(Gene.Symbol))

# shuffled list
posListShuffle = posList %>% 
    mutate(Probe = sample(Probe),Gene.Symbol = rev(Gene.Symbol))

cat(posList$Probe,file = 'posList.grp',sep='\n')
cat(posListRev$Probe,file = 'posListRev.grp',sep='\n')
cat(posListShuffle$Probe,file = 'posListRandom.grp',sep='\n')

# down genes -----------------
negList = hitlist %>%
    arrange((logFC_E12_3_day_IP_vs_naive_3_days_IP)) %>% 
    filter(logFC_E12_3_day_IP_vs_naive_3_days_IP<0) %>% 
    dplyr::select(V1) %>%unlist %>% mouse2human %>% {.$humanGene} %>% unique %>% ogbox::gemmaProbesetMatch('data-raw/GemmaAnnots/GPL96') %>% {.[1:250,]}

negListRev = negList %>% 
    mutate(Probe = rev(Probe),Gene.Symbol = rev(Gene.Symbol))
negListShuffle = negList %>% 
    mutate(Probe = sample(Probe),Gene.Symbol = rev(Gene.Symbol))

cat(negList$Probe,file = 'negList.grp',sep='\n')
cat(negListRev$Probe,file = 'negListRev.grp',sep='\n')
cat(negListShuffle$Probe,file = 'negListRandom.grp',sep='\n')


# all these files return the same output for all instances


# calculation attempt ---------------------------------------------------

# "Let n be the number of probesets in the featureset (22,283)"
n = rankMatrix %>% nrow

# "t be the number of probe sets (or tags) in the tag list." Both 250 but let's assume it can change
tPos = nrow(posList)
tNeg = nrow(negList)

# "Order all n probe sets by the extent of their differential expression for the current instance i"
# the data itself is in ranks so it is already ordered but lets do what they say
orderedProbes = rownames(rankMatrix)[rankMatrix$inst_868 %>% order]

# "Construct a vector V of the position (1... n) of each probe set in the tag
# list in the ordered list of all probe sets and sort these components in
# ascending order such that V(j) is the position of tag j, where j = 1, 2, ..., t." 
# This is the tricky part. Position of tag j itselt is embedded in the formula
# which causes order of probesets to matter
Vpos = match(posList$Probe,orderedProbes) %>% sort
Vneg = match(negList$Probe,orderedProbes) %>% sort

# "compute the following 2 values" I don't really need the t here but wanted to keep the notation consistent
a = function(V,t){
    1:t %>% sapply(function(j){
        g=j/t - V[j]/n
    }) %>% max
}


b = function(V,t){
    1:t %>% sapply(function(j){
        V[j]/n - (j-1)/t
    }) %>% max
}

aPos = a(Vpos,tPos)
bPos = b(Vpos,tPos)

aNeg = a(Vneg,tNeg)
bNeg = b(Vneg,tNeg)

# "set ks= a if (a>b), set ks = -b if b > a"
ks = function(a,b){
    if(a>b){
        return(a)
    } else if (b>a){
        return(-b)
    }
}
kPos = ks(aPos,bPos) # this is supposed to be UP (.133)
kNeg = ks(aNeg,bNeg) # this is supposed to be down (-.188) at least the signs are right I guess

# "The connectivity score Si is set to zero where ksiup and ksidown have the same sign. Otherwise, set si to be ksiup - ksidown,"
score = function(kPos,kNeg){
    if(sign(kPos) == sign(kNeg)){
        return(0)
    } else {
        kPos - kNeg
    }
}

score(kPos,kNeg) # this is supposed to be the score (.856)

