setwd('/Users/zagidull/Documents/elli_projects/repli_seq/raw_data_repli_seq')
library('tidyverse')
library('gridExtra')
library('splines')
library('broom')

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
## rather weak, but significant correlation between s50 and length in insertions. p=0.0345 with 95%CI (0.02434869, 0.55261412). 
## There are no stat signif correlations between three variables in deletions. 
cor(ins[,c('diff','mutatedSamplesFrc','s50')])
cor.test(ins$diff, ins$s50) 
cor(del[,c('diff','mutatedSamplesFrc','s50')])

# re: replication timing and mutatedSamplesFrc
## there is a weak (0.09) positive association between s50 and instability in deletions
## there is a weak negative (-0.17) negative association between s50 and instability in insertions
## but neither are significant
cor.test(del$s50, del$mutatedSamplesFrc)
cor.test(ins$s50, ins$mutatedSamplesFrc)

# for correlation of binary to continuous 
# ToDo stratify by Pol II status, and then do a t-test of genomic features.
# let's fit log reg using Pol2 status as response
# no need to standardize the data, since non-regularized log reg
## neither in insertions, nor in deletions pol2 status can be significantly associated with mutational frequency
logreg_ins_mutatedSamplesFrc <- glm(rnapol2~mutatedSamplesFrc, data = ins, family = 'binomial')
logreg_del_mutatedSamplesFrc <- glm(rnapol2~mutatedSamplesFrc, data = del, family = 'binomial')
summary(logreg_ins_mutatedSamplesFrc)
summary(logreg_del_mutatedSamplesFrc)

ins %>% group_by(rnapol2) %>%
    summarise(m_mutated = mean(mutatedSamplesFrc),
              m_diff = mean(diff),
              m_s50 = mean(s50))


# USING Log.Reg to find association between rnapol2 vs length
## significant association in del for pol2 status vs gene length
## for every extra bp (of gene length) we get 0.00086% higher chance of it having positive Pol II status with p-val of <0.001
logreg_ins_diff <- glm(rnapol2~diff, data = ins, family = 'binomial')
summary(logreg_ins_diff)
t.test(ins[ins$rnapol2=='yes','diff'], ins[ins$rnapol2=='no','diff'])


logreg_del_diff <- glm(rnapol2~diff, data = del, family = 'binomial')
summary(logreg_del_diff)
t.test(del[del$rnapol2=='yes','diff'], del[del$rnapol2=='no','diff'])

summary(logreg_ins_diff)
summary(logreg_del_diff)
# odds from log-odds
exp(cbind(OR = coef(logreg_del_diff), confint(logreg_del_diff)))

# rnapol2 vs s50
## borderline significance of s50 for both ins / del. But let's assume that due to limited samples size it is not. at least for now
logreg_ins_s50 <-glm(rnapol2~s50, data = ins, family = 'binomial')
logreg_del_s50 <-glm(rnapol2~s50, data = del, family = 'binomial')
summary(logreg_ins_s50)
summary(logreg_del_s50)

# actually s50+diff makes models discernibly better for del, vs just s50
logreg_ins_s50diff <-glm(rnapol2~s50+diff, data = ins, family = 'binomial')
logreg_del_s50diff <-glm(rnapol2~s50+diff, data = del, family = 'binomial')
summary(logreg_ins_s50diff)
summary(logreg_del_s50diff)
anova(logreg_del_diff, logreg_del_s50diff, test = 'Chisq') 

# USING t.test to find association between rnapol2 vs length
## suggested by Esa
ttestPolII <- function(genomic_feature){
    p_del <- round(t.test(del[del$rnapol2=='yes', genomic_feature], del[del$rnapol2=='no', genomic_feature])$p.value, 3)
    p_ins <- round(t.test(ins[ins$rnapol2=='yes', genomic_feature], ins[ins$rnapol2=='no', genomic_feature])$p.value, 3)
    return(c(p_del, p_ins))
}

meanPolII <- function(){
    df_del <- del %>%
        group_by(rnapol2) %>%
        summarise(
            m_diff=mean(diff),
            m_s50=mean(s50),
            m_mutatedSamplesFrc=mean(mutatedSamplesFrc),
            .groups = 'drop_last'
        ) %>%
        mutate(label='del', .before='rnapol2') %>%
        as.data.frame
    
    df_ins <- ins %>%
        group_by(rnapol2) %>%
        summarise(
            m_diff=mean(diff),
            m_s50=mean(s50),
            m_mutatedSamplesFrc=mean(mutatedSamplesFrc),
            .groups = 'drop_last'
        ) %>%
        mutate(label='ins', .before='rnapol2') %>%
        as.data.frame
    return(rbind(df_del, df_ins))
}

## there is a difference between gene lengths between Pol II positive and negative genes with p-val=0.004 for MSI-del. 
## PolII positive are shorter (mean length=84075), PolII negative are longer (mean length=215305) 
## p-val (0.65) is not significant for gene length differences in the MSI-ins group 

f_names <- c('diff','s50','mutatedSamplesFrc')
meanPolII()

# t-test p-vals
## multiple hypothesis testing is not needed since we have different data. To check again
sapply(f_names, ttestPolII) %>%
    as.data.frame() %>% 
    mutate(label=c('del','ins'), .before='diff') 



# hist of s50. However, spread is better shown using mean/sd or median/iqr
grid.arrange(ggplot(ins, aes(x=s50)) + geom_histogram(), 
             ggplot(del, aes(x=s50)) + geom_histogram())



# two-sided t-test with H0 of means being equal
## rejected with p-val of 0.036
t.test(ins$s50, del$s50)

