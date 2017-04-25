library(data.table)
library(ogbox)
library(XLConnect)
library(GEOquery)
library(dplyr)
library(magrittr)
library(devtools)
# "complete genelist" is taken from erna directly
hitlist = fread('data-raw/Complete_geneList_voom_limma.csv',data.table = FALSE)
use_data(hitlist)


# cmap downloads do not actually require login

# the first file is the ranked list of DE genes
download.file('ftp://ftp.broad.mit.edu/pub/cmap/rankMatrix.txt.zip',destfile = 'data-raw/rankMatrix.txt.zip')
unzip('data-raw/rankMatrix.txt.zip',exdir = 'data-raw/')
rankMatrix = fread('data-raw/rankMatrix.txt',data.table=FALSE)
# first line causes problems when reading. not sure why
colnames = ogbox::checkLines('data-raw/rankMatrix.txt',lines =1) %>% {.[[1]]} %>% strsplit('\t') %>% {.[[1]]}
cn(rankMatrix) = colnames
use_data(rankMatrix)

# this file has the description of the datasets
downoad.file('https://portals.broadinstitute.org/cmap/cmap_instances_02.xls', destfile = 'data-raw/cmap_instances_02.xls')

instances = loadWorkbook('data-raw/cmap_instances_02.xls')
instances = readWorksheet(instances, sheet = 1, header = TRUE)

# last few lines are footnotes
instances = instances[1:(nrow(instances)-9),]

geoDictionary = c("HG-U133A"= "GPL96",
                  'HT_HG-U133A' = "GPL3921",
                  'HT_HG-U133A_EA' = 'GPL3921')# they use an early access version of the chip with a 676 extra probes (https://portals.broadinstitute.org/cmap/help_topics_linkified.jsp). the probe data for this chip is not available

instances %<>% mutate(GPL = replaceElement(array3,geoDictionary) %>% {.$newVector})
use_data(instances)


dir.create('data-raw/GemmaAnnots')
getGemmaAnnot('GPL96',chipFile = 'data-raw/GemmaAnnots/GPL96')
getGemmaAnnot('GPL3921',chipFile = 'data-raw/GemmaAnnots/GPL3921')


