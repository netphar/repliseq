correlations
================

This notebook contains correlations analysis.

``` r
rm(list=ls(all=TRUE))
library('tidyverse')
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library('broom')
```

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

``` r
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

``` r
ins <- inner_join(len_ins, pol2_ins, by = 'genename') %>%
    inner_join(., msi_ins, by = 'genename') %>%
    inner_join(., s50_ins, by = c('genename' = 'genes')) %>%
    as.data.frame

del <- inner_join(len_del, pol2_del, by = 'genename') %>%
    inner_join(., msi_del, by = 'genename') %>%
    inner_join(., s50_del, by = c('genename' = 'genes')) %>%
    as.data.frame
```

``` r
rm(chrX_del, chrX_ins, len_del, len_ins, msi_del, msi_ins, pol2_del, pol2_ins, s50_del, s50_ins)
```

``` r
# correlations
cor(ins[,c('diff','mutatedSamplesFrc','s50')])
```

    ##                          diff mutatedSamplesFrc        s50
    ## diff               1.00000000       -0.09250722  0.5131294
    ## mutatedSamplesFrc -0.09250722        1.00000000 -0.1159092
    ## s50                0.51312939       -0.11590918  1.0000000

``` r
cor.test(ins$diff, ins$s50) 
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  ins$diff and ins$s50
    ## t = 5.4136, df = 82, p-value = 6.011e-07
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.3356611 0.6554202
    ## sample estimates:
    ##       cor 
    ## 0.5131294

``` r
cor(del[,c('diff','mutatedSamplesFrc','s50')])
```

    ##                         diff mutatedSamplesFrc        s50
    ## diff              1.00000000        0.03568217 0.28707573
    ## mutatedSamplesFrc 0.03568217        1.00000000 0.04491564
    ## s50               0.28707573        0.04491564 1.00000000

``` r
cor.test(del$diff, del$s50) 
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  del$diff and del$s50
    ## t = 7.4743, df = 622, p-value = 2.647e-13
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2133951 0.3575091
    ## sample estimates:
    ##       cor 
    ## 0.2870757

``` r
## OLD result for top50 genes
## rather weak, but significant correlation between s50 and length in insertions. p=0.0312 with 95%CI (0.02434869, 0.55261412). 
## There are no stat signif correlations between three variables in deletions. 
```

### Correlation result, all together. For all unstable MSI genes:

-   a correlation of 0.5131294 between s50 and length in insertions with
    p=6.011e-07 and 95%CI (0.3356611, 0.6554202)
-   a correlation of 0.2870757 between s50 and length in deletions with
    p=2.647e-13 and 95%CI (0.2133951, 0.3575091)

``` r
cor.test(del$s50, del$mutatedSamplesFrc)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  del$s50 and del$mutatedSamplesFrc
    ## t = 1.1213, df = 622, p-value = 0.2626
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.033692  0.122971
    ## sample estimates:
    ##        cor 
    ## 0.04491564

``` r
cor.test(ins$s50, ins$mutatedSamplesFrc)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  ins$s50 and ins$mutatedSamplesFrc
    ## t = -1.0567, df = 82, p-value = 0.2937
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3222952  0.1009958
    ## sample estimates:
    ##        cor 
    ## -0.1159092

``` r
## for top50
## there is a weak (0.1) positive association between s50 and instability in deletions
## there is a weak negative (-0.17) negative association between s50 and instability in insertions
## but neither are significant
```

### S50 vs num mutated samples for all genes

-   there is a weak (0.044) positive association between s50 and
    mutatedSamplesFrc in deletions
-   there is a weak negative (-0.11) negative association between s50
    and mutatedSamplesFrc in insertions but neither are significant,
    p=0.2626 and p=0.2937 respectively

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
