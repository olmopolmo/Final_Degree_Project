library(dplyr)
library(seqinr)
setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/")
fig2<-read.csv("fig2e_fastqs-info.csv")
fig4<-read.csv("fig4_fastqs-info.csv")
fig10<-read.csv("fig10_fastqs-info.csv")


fig10[165,]

#FIGURE 10

fig10 %>% filter(grepl('HDR', Experiment.Title)) #163 HDR
fig10 %>% filter(grepl('HDR', Experiment.Title)) %>% filter(grepl('to',pegRNA)) #138 substitutions
fig10 %>% filter(grepl('HDR', Experiment.Title)) %>% filter(grepl('ins',pegRNA)) #25 insertions


fig10 %>% filter(grepl('PE3', Experiment.Title)) #139 PE
fig10 %>% filter(grepl('PE3', Experiment.Title)) %>% filter(grepl('to',pegRNA)) #108 substitutions
fig10 %>% filter(grepl('PE3', Experiment.Title)) %>% filter(grepl('ins',pegRNA)) #31 insertions

#Figure 4
fig4 %>% filter(grepl('', Experiment.Title)) #294 PE
fig4 %>%  filter(grepl('to',pegRNA)) #210 substitutions
fig4 %>% filter(grepl('ins',pegRNA)) #42 insertions
fig4 %>% filter(grepl('del',pegRNA)) #42 deletions


#TOTAL:
#HDR: 187, 138 subs, 25 ins, 24 dels 
#PE:  433, 318 subs, 42 ins, 42 dels

unique(fig10$gRNA)
"GGCCCAGACTGAGCACGTGA"  "GAGTCCGAGCAGAAGAAGAA"  "GGAATCCCTTCTGCAGCACC"  "GTCATCTTAGTCATTACCTG"  "GCATGGTGCACCTGACTCCTG" "GCAGTGGTGGGGGGCCTTGG" 