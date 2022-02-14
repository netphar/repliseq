setwd('/Users/zagidull/Documents/elli_projects/repli_seq/raw_data_repli_seq')
library('tidyverse')

# chrX gene names. To be excluded
chrX_ins <- c("Dmd", "Nlgn3", "Il13ra2", "Gpc4")
chrX_del <- c("Huwe1", "Phka1", "Mcf2", "Zbtb33")

# reading in raw data + preparing for joining by gene name
pol2_ins <- read_tsv('insHitsRNApol2statusNEW.txt',
                     col_types = cols(.default = col_character(), 
                                      rnapol2 = col_factor())
                     ) %>%
    filter(!genename %in% chrX_ins) %>%
    as.data.frame()

pol2_del <- read_tsv('delHitsRNApol2statusNEW.txt',
                     col_types = cols(.default = col_character(), 
                                      rnapol2 = col_factor())
                     ) %>%
    filter(!genename %in% chrX_del) %>%
    as.data.frame()

s50_ins <- read_tsv('insertions_top50_NEW.csv',
                    col_types = cols(.default = col_character(), 
                                     s50 = col_double())
                    ) %>%
    select(genes, s50) %>%
    filter(!genes %in% chrX_del) %>%
    as.data.frame()

s50_del <- read_tsv('deletions_top50_NEW.csv',
                    col_types = cols(.default = col_character(), 
                                     s50 = col_double())
                    ) %>%
    select(genes, s50) %>%
    filter(!genes %in% chrX_ins) %>%
    as.data.frame()

len_ins <- read_tsv('hotspotInsertionsGeneCoordinatesNEW_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     start = col_double(), 
                                     end = col_double())
                    ) %>%
    select(genename, start, end) %>%
    filter(!genename %in% chrX_ins) %>%
    mutate(diff = end-start) %>%
    select(-c(start,end)) %>%
    as.data.frame()

len_del <- read_tsv('hotspotDeletionsGeneCoordinatesNEW_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     start = col_double(), 
                                     end = col_double())
                    ) %>%
    select(genename, start, end) %>%
    filter(!genename %in% chrX_del) %>%
    mutate(diff = end-start) %>%
    select(-c(start,end)) %>%
    as.data.frame()
    
msi_ins <- read_tsv('MSI_ins_Hotspot_for_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     mutatedSamplesFrc = col_double())
                    ) %>%
    filter(!genename %in% chrX_ins) %>%
    as.data.frame()

msi_del <- read_tsv('MSI_del_Hotspot_for_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     mutatedSamplesFrc = col_double())
                    ) %>%
    filter(!genename %in% chrX_del) %>%
    as.data.frame()

# join by gene name
ins <- inner_join(len_ins, pol2_ins, by = 'genename') %>%
    inner_join(., msi_ins, by = 'genename') %>%
    inner_join(., s50_ins, by = c('genename' = 'genes')) %>%
    as.data.frame

del <- inner_join(len_del, pol2_del, by = 'genename') %>%
    inner_join(., msi_del, by = 'genename') %>%
    inner_join(., s50_del, by = c('genename' = 'genes')) %>%
    as.data.frame

# data clean
rm(chrX_del, chrX_ins, len_del, len_ins, msi_del, msi_ins, pol2_del, pol2_ins, s50_del, s50_ins)

# correlations
## rather weak correlation between s50 and length in insertions (absent in deletions. 
## There are virtually no correlations between three variables in deletions. 
cor(ins[,c('diff','mutatedSamplesFrc','s50')])
cor(del[,c('diff','mutatedSamplesFrc','s50')])

# for correlation of binary to continuous 
# let's fit log reg using Pol2 status as response
# no need to standardize the data, since non-regularized log reg
## neither in insertions, nor in deletions pol2 status can be significantly associated with mutational frequency
summary(glm(rnapol2~mutatedSamplesFrc, data = ins, family = 'binomial'))
summary(glm(rnapol2~mutatedSamplesFrc, data = del, family = 'binomial'))


