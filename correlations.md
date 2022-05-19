correlations
================

This notebook contains correlations analysis.

``` r
rm(list=ls(all=TRUE))
library('tidyverse')
library('broom')
```

Let’s find ChrX genes. They are excluded from analysis

``` r
# chrX gene names. To be excluded
#chrX_ins <- c("Dmd", "Nlgn3", "Il13ra2", "Gpc4") <- old, for top50 names
#chrX_del <- c("Huwe1", "Phka1", "Mcf2", "Zbtb33") <- old, for top50 names

chrX_del <- read_tsv('MSIdelGeneCoordinatesNEW_bulat.txt',
                     col_types = cols(.default = col_character())
                     ) %>% 
    filter(chr=='chrX') %>%
    select(genename) %>%
    flatten_chr

chrX_ins <- read_tsv('MSIinsGeneCoordinatesNEW_bulat.txt',
                     col_types = cols(.default = col_character())
                     ) %>% 
    filter(chr=='chrX') %>%
    select(genename) %>%
    flatten_chr   
```

load data

``` r
# reading in raw data + preparing for joining by gene name
pol2_ins <- read_tsv('MSIinskoGenesRNApol2statusNEW_NEW.txt',
                     col_types = cols(.default = col_character(), 
                                      rnapol2 = col_factor())
                     ) %>%
    filter(!genename %in% chrX_ins) %>%
    as.data.frame

pol2_del <- read_tsv('MSIdelkoGenesRNApol2statusNEW_NEW.txt',
                     col_types = cols(.default = col_character(), 
                                      rnapol2 = col_factor())
                     ) %>%
    filter(!genename %in% chrX_del) %>%
    as.data.frame

s50_ins <- read_tsv('insertions_all_090522.csv',
                    col_types = cols(.default = col_character(), 
                                     s50 = col_double())
                    ) %>%
    select(genes, s50) %>%
    filter(!genes %in% chrX_del) %>%
    as.data.frame

s50_del <- read_tsv('deletions_all_090522.csv',
                    col_types = cols(.default = col_character(), 
                                     s50 = col_double())
                    ) %>%
    select(genes, s50) %>%
    filter(!genes %in% chrX_ins) %>%
    as.data.frame

len_ins <- read_tsv('MSIinsGeneCoordinatesNEW_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     start = col_double(), 
                                     end = col_double())
                    ) %>%
    select(genename, start, end) %>%
    filter(!genename %in% chrX_ins) %>%
    mutate(diff = end-start) %>%
    select(-c(start,end)) %>%
    as.data.frame

len_del <- read_tsv('MSIdelGeneCoordinatesNEW_bulat.txt',
                    col_types = cols(.default = col_character(), 
                                     start = col_double(), 
                                     end = col_double())
                    ) %>%
    select(genename, start, end) %>%
    filter(!genename %in% chrX_del) %>%
    mutate(diff = end-start) %>%
    select(-c(start,end)) %>%
    as.data.frame
    
msi_ins <- read_tsv('MSIinskoTargetGenesNEW.txt',
                    col_types = cols(.default = col_character(), 
                                     mutatedSamplesFrc = col_double())
                    ) %>%
    select(c(genename,mutatedSamplesFrc)) %>%
    filter(!genename %in% chrX_ins) %>%
    as.data.frame

msi_del <- read_tsv('MSIdelkoTargetGenesNEW.txt',
                    col_types = cols(.default = col_character(), 
                                     mutatedSamplesFrc = col_double())
                    ) %>%
    select(c(genename,mutatedSamplesFrc)) %>%
    filter(!genename %in% chrX_del) %>%
    as.data.frame
```

join all by gene names. results in two files - ins and del

``` r
ins <- inner_join(len_ins, pol2_ins, by = 'genename') %>%
    inner_join(., msi_ins, by = 'genename') %>%
    inner_join(., s50_ins, by = c('genename' = 'genes')) %>%
    as.data.frame

del <- inner_join(len_del, pol2_del, by = 'genename') %>%
    inner_join(., msi_del, by = 'genename') %>%
    inner_join(., s50_del, by = c('genename' = 'genes')) %>%
    as.data.frame

rm(chrX_del, chrX_ins, len_del, len_ins, msi_del, msi_ins, pol2_del, pol2_ins, s50_del, s50_ins)
```

### Let’s see correlations in ins and dels without any stratification

``` r
tidy_cor <- function(object){
    l <- vector('list', length = 3)
    l[[1]] <- tidy(cor.test(object$diff, object$mutatedSamplesFrc))
    l[[2]] <- tidy(cor.test(object$diff, object$s50))
    l[[3]] <- tidy(cor.test(object$s50, object$mutatedSamplesFrc))
    return(do.call(rbind, l))
}
tb <- tidy_cor(ins)
tb <- add_column(tb, .before = 1, var1=c('diff','diff','s50'))
tb <- add_column(tb, .before = 2, var2=c('mutatedSamplesFrc','s50','mutatedSamplesFrc'))
tb <- add_column(tb, .before = 1, type='ins')

tb1 <- tidy_cor(del)
tb1 <- add_column(tb1, .before = 1, var1=c('diff','diff','s50'))
tb1 <- add_column(tb1, .before = 2, var2=c('mutatedSamplesFrc','s50','mutatedSamplesFrc'))
tb1 <- add_column(tb1, .before = 1, type='del')
bind_rows(tb, tb1)
```

    ## # A tibble: 6 × 11
    ##   type  var1  var2      estimate statistic  p.value parameter conf.low conf.high
    ##   <chr> <chr> <chr>        <dbl>     <dbl>    <dbl>     <int>    <dbl>     <dbl>
    ## 1 ins   diff  mutatedS…  -0.0925    -0.841 4.03e- 1        82  -0.301      0.124
    ## 2 ins   diff  s50         0.513      5.41  6.01e- 7        82   0.336      0.655
    ## 3 ins   s50   mutatedS…  -0.116     -1.06  2.94e- 1        82  -0.322      0.101
    ## 4 del   diff  mutatedS…   0.0357     0.890 3.74e- 1       622  -0.0429     0.114
    ## 5 del   diff  s50         0.287      7.47  2.65e-13       622   0.213      0.358
    ## 6 del   s50   mutatedS…   0.0449     1.12  2.63e- 1       622  -0.0337     0.123
    ## # … with 2 more variables: method <chr>, alternative <chr>

-   the only positively & significantly correlated pair in insertions is
    diff vs s50 (0.513)
-   the only positively & significantly correlated pair in deletions is
    diff vs s50 (0.287)

### Let’s see what happens when we stratify by PolII status

``` r
ttestPolII <- function(genomic_feature){
    # genomic feature is the colname, so one of c('diff','s50','mutatedSamplesFrc')
    print('pvals')
    p_del <- t.test(del[del$rnapol2=='yes', genomic_feature], del[del$rnapol2=='no', genomic_feature])$p.value
    p_ins <- t.test(ins[ins$rnapol2=='yes', genomic_feature], ins[ins$rnapol2=='no', genomic_feature])$p.value
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

## OLD for top50 genes
## there is a difference between gene lengths between Pol II positive and negative genes with p-val=0.002 for MSI-del. PolII positive are shorter (mean length=62549), PolII negative are longer (mean length=188772) 
## p-val (0.65) is not significant for gene length differences in the MSI-ins group 

f_names <- c('diff','s50','mutatedSamplesFrc')
meanPolII()
```

    ## [1] "actual scores"

    ##   label rnapol2  m_diff     m_s50 m_mutatedSamplesFrc
    ## 1   del     yes 51689.0 0.3208069          0.05475504
    ## 2   del      no 78777.0 0.4106137          0.05677716
    ## 3   ins     yes 56564.5 0.2867391          0.05039526
    ## 4   ins      no 72367.0 0.3634211          0.04665072

``` r
# t-test p-vals
## multiple hypothesis testing is not needed since we have different data
sapply(f_names, ttestPolII) %>%
    as.data.frame() %>% 
    mutate(label=c('del','ins'), .before='diff') 
```

    ## [1] "pvals"
    ## [1] "pvals"
    ## [1] "pvals"

    ##   label         diff          s50 mutatedSamplesFrc
    ## 1   del 8.966018e-06 6.700679e-10         0.2954364
    ## 2   ins 2.271116e-02 4.740439e-02         0.1270231

### Pol II vs gene length / s50 / num mutated samples. For all genes

#### FOR MSI-del

-   there is a difference between gene lengths in Pol II positive and
    negative genes,
-   PolII positive are shorter (length=51689), PolII negative are longer
    (length=78777), with p-val=8.966018e-06
-   there is a difference between s50 in Pol II positive and negative
    genes,
-   PolII positive have lower s50 (s50=0.321), PolII negative are higher
    s50 (s50=0.411), with p-val=6.700679e-10

#### FOR MSI-ins.

-   there is a difference between gene lengths in Pol II positive and
    negative genes,
-   PolII positive are shorter (length=56565), PolII negative are longer
    (length=72367), with p-val=0.02271116
-   there is a difference between s50 in Pol II positive and negative
    genes,
-   PolII positive have lower s50 (s50=0.287), PolII negative have
    higher s50 (s50=0.363), with p-val=0.04740439
-   No significant differences in mutated Samples Frc for either
