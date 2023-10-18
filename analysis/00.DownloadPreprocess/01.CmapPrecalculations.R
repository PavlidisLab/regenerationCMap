# save the randomized matrices for quick calculation
library(memoise)
library(ConnectivityMap)
devtools::load_all()
data("instances")
data("rankMatrix")
d = 100000
randomV = function(length,d){
    replicate(d,sample(x=1:ncol(rankMatrix),
                       size = length,replace = FALSE) %>% sort)
}
memoRandomV = memoise(randomV)
memoKsCalc = memoise(ksCalc)

instanceLengths = table(instances$cmap_name)


allVRandoms = instanceLengths %>% unique %>% sapply(function(x){
    print(x)
    return(memoRandomV(x,d))
})

allVRandoms %>% sapply(function(x){
    memoKsCalc(x,nrow(instances))
})

use_data(memoRandomV,overwrite = TRUE)
use_data(memoKsCalc,overwrite = TRUE)