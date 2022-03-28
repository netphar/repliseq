library('tidyverse')
library('gridExtra')
library('splines')
library('broom')

# chrX gene names. To be excluded
chrX_ins <- c("Dmd", "Nlgn3", "Il13ra2", "Gpc4")
chrX_del <- c("Huwe1", "Phka1", "Mcf2", "Zbtb33")

# reading in raw data + preparing for joining by gene name
pol2_ins <- read_tsv('KOinsHitsRNApol2statusNEW_NEW.txt',
                     col_types = cols(.default = col_character(), 
                                      rnapol2 = col_factor())
                     ) %>%
    filter(!genename %in% chrX_ins) %>%
    as.data.frame()

pol2_del <- read_tsv('KOdelHitsRNApol2statusNEW_NEW.txt',
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
## rather weak, but significant correlation between s50 and length in insertions. p=0.0312 with 95%CI (0.02434869, 0.55261412). 
## There are no stat signif correlations between three variables in deletions. 
cor(ins[,c('diff','mutatedSamplesFrc','s50')])
cor.test(ins$diff, ins$s50) 
cor(del[,c('diff','mutatedSamplesFrc','s50')])

# re: replication timing and mutatedSamplesFrc
## there is a weak (0.1) positive association between s50 and instability in deletions
## there is a weak negative (-0.17) negative association between s50 and instability in insertions
## but neither are significant
cor.test(del$s50, del$mutatedSamplesFrc)
cor.test(ins$s50, ins$mutatedSamplesFrc)

# USING t.test to find association between rnapol2 vs length
## suggested by Esa

ttestPolII <- function(genomic_feature){
    # genomic feature is the colname, so one of c('diff','s50','mutatedSamplesFrc')
    print('pvals')
    p_del <- round(t.test(del[del$rnapol2=='yes', genomic_feature], del[del$rnapol2=='no', genomic_feature])$p.value, 3)
    p_ins <- round(t.test(ins[ins$rnapol2=='yes', genomic_feature], ins[ins$rnapol2=='no', genomic_feature])$p.value, 3)
    return(c(p_del, p_ins))
}

meanPolII <- function(){
    print('actual scores')
    df_del <- del %>%
        group_by(rnapol2) %>%
        summarise(
            m_diff=median(diff),
            m_s50=mean(s50),
            m_mutatedSamplesFrc=mean(mutatedSamplesFrc),
            .groups = 'drop_last'
        ) %>%
        mutate(label='del', .before='rnapol2') %>%
        as.data.frame
    
    df_ins <- ins %>%
        group_by(rnapol2) %>%
        summarise(
            m_diff=median(diff),
            m_s50=mean(s50),
            m_mutatedSamplesFrc=mean(mutatedSamplesFrc),
            .groups = 'drop_last'
        ) %>%
        mutate(label='ins', .before='rnapol2') %>%
        as.data.frame
    return(rbind(df_del, df_ins))
}

## there is a difference between gene lengths between Pol II positive and negative genes with p-val=0.002 for MSI-del. PolII positive are shorter (mean length=62549), PolII negative are longer (mean length=188772) 
## p-val (0.65) is not significant for gene length differences in the MSI-ins group 

f_names <- c('diff','s50','mutatedSamplesFrc')
meanPolII()

# t-test p-vals
## multiple hypothesis testing is not needed since we have different data. To check again
sapply(f_names, ttestPolII) %>%
    as.data.frame() %>% 
    mutate(label=c('del','ins'), .before='diff') 