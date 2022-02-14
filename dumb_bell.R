setwd('/Users/zagidull/Documents/elli_projects/repli_seq/raw_data_repli_seq')
library('tidyverse')
library('ggalt')
library('gridExtra')

df = read.csv('dataForFigureBulat.txt', sep = '\t')
xlim_del=c(6, 20)
xlim_ins=c(6, 20)
cutoff=0.05
num_breaks=10

plot_WT_del <- df %>%
    filter(gt=='wt' & type=='del') %>%
    select(ms_length, new_ms_length, count) %>%
    mutate(perc=ifelse(test=count/sum(count)>cutoff, yes = paste0(round(100*count/sum(count),1),'%'),no='')) %>%
   #arrange(desc(perc)) %>%
    #filter(perc > 10) %>%
    ggplot(aes(x=ms_length, xend=new_ms_length, y=1:nrow(.))) + 
    theme_minimal() + 
    theme(text = element_text(size=20),
          plot.title = element_text(colour = "#FAAB18"),
          axis.title.x = element_text(colour = '#222222')) +
    geom_dumbbell(colour = "#dddddd",dot_guide = FALSE,
                  size = 1.5,
                  colour_x = "#FAAB18",
                  colour_xend = "#1380A1",
                  size_x = 3,
                  size_xend = 1.5) + 
    geom_text(aes(label=perc, hjust=-.5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=num_breaks)) +
    scale_y_discrete(labels = NULL, breaks = NULL) + 
    coord_cartesian(xlim = xlim_del, clip = 'off') +
    labs(x='',y='',title = 'WT/deletions')

plot_KO_del <- df %>%
    filter(gt=='ko' & type=='del') %>%
    select(ms_length, new_ms_length, count) %>%
    mutate(perc=ifelse(test=count/sum(count)>cutoff, yes = paste0(round(100*count/sum(count),1),'%'),no='')) %>%
    ggplot(aes(x=ms_length, xend=new_ms_length, y=1:nrow(.))) + 
    theme_minimal() + 
    theme(text = element_text(size=20),
          plot.title = element_text(colour = "#FAAB18"),
          axis.title.x = element_text(colour = '#222222')) +
    geom_dumbbell(colour = "#dddddd",dot_guide = FALSE,
                  size = 1.5,
                  colour_x = "#FAAB18",
                  colour_xend = "#1380A1",
                  size_x = 3,
                  size_xend = 1.5) + 
    geom_text(aes(label=perc, hjust=-.5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=num_breaks)) +
    scale_y_discrete(labels = NULL, breaks = NULL) + 
    coord_cartesian(xlim = xlim_del, clip = 'off') +
    labs(x='',y='',title = 'KO/deletions')

plot_WT_ins <- df %>%
    filter(gt=='wt' & type=='ins') %>%
    select(ms_length, new_ms_length, count) %>%
    mutate(perc=ifelse(test=count/sum(count)>cutoff, yes = paste0(round(100*count/sum(count),1),'%'),no='')) %>%
    ggplot(aes(x=ms_length, xend=new_ms_length, y=1:nrow(.) )) + 
    theme_minimal() + 
    theme(text = element_text(size=20),
          plot.title = element_text(colour = '#1380A1'),
          axis.title.x = element_text(colour = '#222222')) +
    geom_dumbbell(colour = "#dddddd",
                  size = 1.5,
                  colour_x = "#1380A1",
                  colour_xend = "#FAAB18",
                  size_x = 3,
                  size_xend = 1.5) + 
    geom_text(aes(label=perc, hjust=1.5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=num_breaks)) +
    scale_y_discrete(labels = NULL, breaks = NULL) + 
    coord_cartesian(xlim = xlim_ins, clip = 'off') +
    labs(x='',y='',title = 'WT/insertions')

plot_KO_ins <- df %>%
    filter(gt=='ko' & type=='ins') %>%
    select(ms_length, new_ms_length, count) %>%
    mutate(perc=ifelse(test=count/sum(count)>cutoff, yes = paste0(100*round(count/sum(count),1),'%'),no='')) %>%
    ggplot(aes(x=ms_length, xend=new_ms_length, y=1:nrow(.) )) + 
    theme_minimal() + 
    theme(text = element_text(size=20),
          plot.title = element_text(colour = "#1380A1"),
          axis.title.x = element_text(colour = '#222222')) +
    geom_dumbbell(colour = "#dddddd",
                  size = 1.5,
                  colour_x = "#1380A1",
                  colour_xend = "#FAAB18",
                  size_x = 3,
                  size_xend = 1.5) + 
    geom_text(aes(label=perc, hjust=1.5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=num_breaks)) +
    scale_y_discrete(labels = NULL, breaks = NULL) +
    coord_cartesian(xlim = xlim_ins, clip = 'off') +
    labs(x='',y='', title = 'KO/insertions')

grid.arrange(plot_WT_del, plot_KO_del, plot_WT_ins, plot_KO_ins, nrow = 2, ncol = 2)

