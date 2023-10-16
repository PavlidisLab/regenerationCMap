
# L1000 precalcs ------------------
print(Sys.getpid())
print('reading data')

L1000geneAnnots = readRDS('data-raw/lincs1000_data/L1000geneAnnots.rds')
inst = readRDS('data-raw/lincs1000_data/instances.rds')

devtools::load_all()
library(cmapQuery)
library(dplyr)
library(magrittr)
calculateKs = TRUE
print('pre-calcing random Ks')
L1000PreCalc = preCalcRandomKs(inst$pert_iname)

saveRDS(L1000PreCalc,'data-raw/lincs1000_data/L1000PreCalc.rds')




# FWD pre-calcs --------------
L1000geneAnnots = readRDS('data-raw/FWD_data/FWDgeneAnnots.rds')
inst = readRDS('data-raw/FWD_data/FWDinstances.rds')

L1000PreCalc = preCalcRandomKs(inst$pert_iname)

saveRDS(L1000PreCalc,'data-raw/FWD_data/FWDPreCalc.rds')



# Pavlab pre-calcs --------
L1000geneAnnots = readRDS('data-raw/pav_data/pavGeneAnnots.rds')
inst = readRDS('data-raw/pav_data/pavInstances.rds')
    
L1000PreCalc = preCalcRandomKs(inst$pert_iname)
saveRDS(L1000PreCalc,'data-raw/pav_data/pavPreCalc.rds')

