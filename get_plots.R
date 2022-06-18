#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(plyr)
library(dplyr)
library(seqinr)
library(ShortRead)
library(ggpubr)
library(ggmsa)
library(seqmagick)
library(stringr)
library(tidyr)
library(ggseqlogo)
library(plotly)
library(cowplot)

######################
##### CRISPR-GA-1 like plot, but with ggplot
######################

plotIndels_gg <- function(indels_data, cut){
  # Split insertions from deletions
  l <- length(indels_data$Modification) + indels_data$wt_reads[1] + data$t_reads[1]
  dels <- indels_data %>% filter(Modification == "del") 
  ins <- indels_data %>% filter(Modification == "ins") 
  
  #Metadata deletion
  if (length(dels$X)>0){
    metadels <- as.data.frame(table(dels$Start,dels$Length))
    metadels <- metadels[which(metadels$Freq>0),]
    colnames(metadels) <- c("Start","length","Freq")
    metadels$End <- 0
    metadels$changes <- sum((as.integer(as.character(metadels$Freq)) * as.integer(as.character(metadels$length)))) #Now the changes column corrspond to a constant that will be used in the creation of accumulative. This variable stores the total number of changes per position
    
    for (i in 1:length(metadels$Start)){
      metadels$End[i] = as.integer(as.character(metadels$Start[i])) + as.integer(as.character(metadels$length[i])) - 1
    }
  }else{
    metadels <- list()
  }
  
  exportJson <- jsonlite::toJSON(metadels)
  write(exportJson, paste(sample_name,"_metadels.json",sep=""))
  
  #Metadata insertion  
  if(length(ins$X)>0){
    metains <- as.data.frame(table(ins$Start,ins$Length))
    metains <- metains[which(metains$Freq>0),]
    colnames(metains) <- c("Start","length","Freq")
    metains$End <- 0
    metains$changes <- sum((as.integer(as.character(metains$Freq)) * as.integer(as.character(metains$length)))) #Now the changes column corrspond to a constant that will be used in the creation of accumulative. This variable stores the total number of changes per position
    
    for (i in 1:length(metains$Start)){
      metains$End[i] = as.integer(as.character(metains$Start[i])) + as.integer(as.character(metains$length[i])) - 1
    }  
  }else{
    metains <- list()
  }
  exportJson <- jsonlite::toJSON(metains)
  write(exportJson, paste(sample_name,"_metains.json",sep=""))
  
  
# # Count accumulated insertions and deletions by position
#   accum.dels=c()
#   accum.dels=unlist(apply(dels,1,function(X) return( seq(as.integer(X[3]), as.integer(X[3]) + as.integer(X[4])) )))
#   accum.ins=c()
#   accum.ins=unlist(apply(ins,1,function(X) return( seq(as.integer(X[3]), as.integer(X[3]) + as.integer(X[4])) )))
#   
#   # Two options: dejarlo como esta o relentizar el proceso con tal de tener un indicador que diga si es start end o intermediate
#   
#   # Plot
#   as.data.frame(table(accum.dels))
#   a_dels <- data.frame(location = c(accum.dels))
#   a_ins <- data.frame(location = c(accum.ins)) 
#   
# 
#   ###############################################
#   #         ACUMMULATIVE LOCATIONS              #
#   ###############################################
#   if ( dim(dels)[1] != 0 ){
#     acc_dels <- as.data.frame(table(a_dels))
#     acc_dels <- acc_dels[acc_dels$Freq>0,]
#     acc_dels$percent <- (acc_dels$Freq/sum(acc_dels$Freq))*100
# #    del_acumm <- plot_ly(y = ~acc_dels$percent, x= ~as.numeric(as.character(acc_dels$a_dels)),type = "bar",hovertemplate = paste('Percentage: %{y}','Counts: ',acc_dels$Freq)) %>%
# #      layout(xaxis = list(title="Accumulative Deletion Location",rangeslider = list()), yaxis = list(fixedrange = FALSE,title = "Reads Percentage(%)"), barmode = 'stack')%>% 
# #      layout(shapes = list(list(
# #        type = "line",
# #        y0 = 0,
# #        y1 = 100,
# #        yref = "paper",
# #        x0 = cut_site,
# #        x1 = cut_site,
# #        line = list(color = "red", dash="dot")
# #      )))%>%
# #      add_text(showlegend = FALSE, x=cut_site, y=max(acc_dels$percent),text = "Cut site", hovertemplate = paste("Cut site ", cut_site),
# #               textfont = list(color= c("red")))
#   }else {
#     acc_dels <- data.frame()
#   }
# 
# 
#   if ( dim(ins)[1] != 0 ){
#     acc_ins <- as.data.frame(table(a_ins))
#     acc_ins <- acc_ins[acc_ins$Freq>0,]
#     acc_ins$percent <- (acc_ins$Freq/sum(acc_ins$Freq))*100
#     
# #    ins_acumm <- plot_ly(y = ~acc_ins$percent, x= ~as.numeric(as.character(acc_ins$a_ins)),type = "bar",hovertemplate = paste('Percentage: %{y}','Counts: ',acc_ins$Freq)) %>%
# #      layout(xaxis = list(title="Accumulative Insertion Location",rangeslider = list()), yaxis = list(fixedrange = FALSE,title = "Reads Percentage(%)"), barmode = 'stack')%>% 
# #      layout(shapes = list(list(
# #        type = "line",
# #        y0 = 0,
# #        y1 = 100,
# #        yref = "paper",
# #        x0 = cut_site,
# #        x1 = cut_site,
# #        line = list(color = "red", dash="dot")
# #      )))%>%
#  #     add_text(showlegend = FALSE, x=cut_site, y=max(acc_ins$percent),text = "Cut site", hovertemplate = paste("Cut Site ", cut_site),
#  #             textfont = list(color= c("red")))
#   }else {
#     acc_ins <- data.frame()
#   }
  
  # if( (dim(ins)[1] != 0) &&  (dim(dels)[1] != 0) ){
  #   accummulative <- subplot(del_acumm, ins_acumm, shareY = T,shareX = F,titleX = T,titleY = T)
  # } else if( dim(ins)[1] != 0 ){
  #   accummulative <- subplot(ins_acumm, shareY = T,shareX = F,titleX = T,titleY = T)
  # } else {
  #   accummulative <- subplot(del_acumm, shareY = T,shareX = F,titleX = T,titleY = T)
  # }
 
  #Elements for dynamic table
  sequece_length <- max(indels_data$Start+indels_data$Length)
  exportJson <- jsonlite::toJSON(sequece_length)
  write(exportJson, paste(sample_name,"_length.json",sep=""))
  
  # modif_end <- indels_data$Start + indels_data$Length
  #
  # exportJson <- jsonlite::toJSON(modif_end)
  # write(exportJson, paste(sample_name,"_end_pos_reads.json",sep=""))
  
  exportJson <- jsonlite::toJSON(l)
  write(exportJson, paste(sample_name,"_Total_reads.json",sep=""))
  
  # exportJson <- jsonlite::toJSON(acc_ins)
  # write(exportJson, paste(sample_name,"_Ins_dynamic.json",sep=""))
  # 
  # exportJson <- jsonlite::toJSON(acc_dels)
  # write(exportJson, paste(sample_name,"_Del_dynamic.json",sep=""))
    
  
  #INSERTIONS
  ###############################################
  #         LOCATION                            #
  ###############################################
  
  pattern = ins
  
  #defensive programming against file with no insertions
  if(length(pattern$Modification)>0){
    
    LOC <- as.data.frame(table(pattern$Start,pattern$Length,pattern$ins_nt))
    LOC <- LOC[which(LOC$Freq>0),]
    LOC$percent <- (LOC$Freq/l)*100
    
    
    #AXIS -> PERCENTAGE, BAR -> COUNT
    Ins_location <- plot_ly() %>%
    layout(xaxis = list(title="Insertion Location",rangeslider = list()), yaxis = list(fixedrange = FALSE,title = "Reads Percentage(%)"), barmode = 'stack')
    
    Ins_location <- Ins_location %>% add_trace(x = ~as.numeric(as.character(LOC$Var1)), 
                                               y= ~LOC$percent, color = ~LOC$Var3, type = 'bar',
                                               marker = list(~LOC$Freq), legendgroup = ~LOC$Var3,
                                               hovertemplate = paste('Percentage: %{y}','Counts: ',LOC$Freq))
    
  
    Ins_location <- Ins_location %>%   layout(shapes = list(list(
      type = "line",
      y0 = 0,
      y1 = 100,
      yref = "paper",
      x0 = cut_site,
      x1 = cut_site,
      line = list(color = "red", dash="dot")
    ))) %>% add_text(showlegend = FALSE, x=cut_site,y=max(LOC$percent),text = "Cut site", hovertemplate = paste("Cut site ", cut_site),
                     textfont = list(color= c("red")))
    
    
    ###############################################
    #                SIZES                        #
    ###############################################
    # SIZE <- as.data.frame(table(pattern$Length,pattern$patterns))
    # SIZE <- SIZE[which(SIZE$Freq>0),]
    # SIZE$percent <- (SIZE$Freq/l)*100
    
    #AXIS -> PERCENTAGE, BAR -> COUNT
    Ins_sizes <- plot_ly() %>%
      layout(xaxis = list(title="Insertion Sizes",rangeslider = list()), yaxis = list(title = "Reads Percentage(%)"), barmode = 'stack')
    Ins_sizes <- Ins_sizes %>% add_trace(x = ~as.numeric(as.character(LOC$Var2)), 
                                         y= ~LOC$percent,  color = ~LOC$Var3, type = 'bar',
                                         marker = list(~LOC$Freq), showlegend = F, legendgroup = ~LOC$Var3,
                                         hovertemplate = paste(' Percentage: %{y}','Count: ',LOC$Freq))
    
    ###############################################
    #                MERGE PLOTS                  #
    ###############################################
    
    insertions <- subplot(Ins_location,Ins_sizes, titleY = F, titleX = T) %>% layout(title = "Insertions",yaxis = list(title = "Reads Percentage (%)"))
  }else{
    insertions<-empty_plot("No Insertions were detected.")
  }
  #DELETIONS
  #pattern filter (A,C,T,G) -> NHEJ
  pattern <- dels
  
  #defensive programming against file with no insertions
  if(length(pattern$Modification)){
      
    for(i in 1:length(pattern$patterns)){
      if(length(strsplit(pattern$patterns[i],split = "")[[1]])==1){
        pattern$patterns[i] <- "NHEJ"
        
      }
    }
    
    ###############################################
    #         LOCATION                            #
    ###############################################
    
    
    LOC <- as.data.frame(table(pattern$Start,pattern$patterns,pattern$Length))
    LOC <- LOC[which(LOC$Freq>0),]
    LOC$percent <- (LOC$Freq/l)*100
    
    
    #AXIS -> PERCENTAGE, BAR -> COUNT
    Del_location <- plot_ly() %>%
    layout(xaxis = list(title="Deletion Location",rangeslider = list()), yaxis = list(fixedrange = FALSE,title = "Reads Percentage(%)"), barmode = 'stack')
    Del_location <- Del_location %>% add_trace(x = ~as.numeric(as.character(LOC$Var1)), 
                                               y= ~LOC$percent, color = ~LOC$Var2, type = 'bar',
                                               marker = list(~LOC$Freq),  legendgroup = ~LOC$Var2,
                                               hovertemplate = paste('Percentage: %{y}','Counts: ',LOC$Freq))
    
    Del_location <- Del_location %>%   layout(shapes = list(list(
      type = "line",
      y0 = 0,
      y1 = 100,
      yref = "paper",
      x0 = cut_site,
      x1 = cut_site,
      line = list(color = "red", dash="dot")
    )))%>%
      add_text(showlegend = FALSE, x=cut_site,y=max(LOC$percent),text = "Cut site", hovertemplate = paste("Cut site ", cut_site),
               textfont = list(color= c("red")))
    
    
    ###############################################
    #                SIZES                        #
    ###############################################
    # SIZE <- as.data.frame(table(pattern$Length,pattern$patterns))
    # SIZE <- SIZE[which(SIZE$Freq>0),]
    # SIZE$percent <- (SIZE$Freq/l)*100
  
    #AXIS -> PERCENTAGE, BAR -> COUNT
    Del_sizes <- plot_ly() %>%
      layout(xaxis = list(title="Deletion Sizes",rangeslider = list()), yaxis = list(title = "Reads Percentage(%)"), barmode = 'stack')
    Del_sizes <- Del_sizes %>% add_trace(x = ~as.numeric(as.character(LOC$Var3)), 
                                         y= ~LOC$percent, color = ~LOC$Var2, type = 'bar',
                                         marker = list(~LOC$Freq), showlegend=F, legendgroup = ~LOC$Var2,
                                         hovertemplate = paste('Size: ',LOC$Var3,' Percentage: %{y}'))
    
    
    ###############################################
    #                MERGE PLOTS                  #
    ###############################################
    
    deletions <- subplot(Del_location,Del_sizes, titleY = F, titleX = T) %>% layout(title = "Deletions",yaxis = list(title = "Reads Percentage (%)"))
  }else{
    deletions<-empty_plot("No deletions were detected.")
  }
  plots <- list(deletions,insertions)
  return(plots)
}

#plotIndels_gg <- function(indels_data, cut){
  # Split insertions from deletions
  #dels <- indels_data %>% filter(Modification == "del") 
  #ins <- indels_data %>% filter(Modification == "ins") 
  # Count accumulated insertions and deletions by position
  #accum.dels=c()
  #accum.dels=unlist(apply(dels,1,function(X) return( seq(as.integer(X[3]), as.integer(X[3]) + as.integer(X[4])) )))
  #accum.ins=c()
  #accum.ins=unlist(apply(ins,1,function(X) return( seq(as.integer(X[3]), as.integer(X[3]) + as.integer(X[4])) )))
  
  # Plot
  #a <-  ggplot(dels, aes(x=as.numeric(Length))) + 
    #geom_bar(color="black", fill="purple") + 
    #theme_classic() + xlab("Deletions sizes") + ylab("Counts") +
    #scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
  #b <-  ggplot(ins, aes(x=as.numeric(Length))) + 
    #geom_bar(color="black", fill="purple") + 
    #theme_classic() + xlab("Insertions sizes") + ylab("") +
    #scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
  # c <- ggplot(s_sizes, aes(x=as.numeric(as.character(Var1)), y=as.numeric(Freq))) + 
  #   geom_bar(stat="identity", color="black", fill="purple") + 
  #   theme_classic() + xlab("Substitutions sizes") + ylab("") +
  #   scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

  #a_dels <- data.frame(location = c(accum.dels))
  #d <- ggplot(a_dels, aes(x=location)) + 
    #geom_bar(color="black", fill="purple") + 
    #geom_vline(xintercept = cut, color="red") +
    #theme_classic() + xlab("Deletions location") + ylab("Counts") +
    #scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
  #a_ins <- data.frame(location = c(accum.ins)) 
  #e <- ggplot(a_ins, aes(x=location)) + 
    #geom_bar(color="black", fill="purple") + 
    #geom_vline(xintercept = cut, color="red") +
    #theme_classic() + xlab("Insertions location") + ylab("") +
    #scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
  # f <- ggplot(s_location, aes(x=as.numeric(as.character(Var1)), y=as.numeric(Freq))) + 
  #   geom_bar(stat="identity", color="black", fill="purple") + 
  #   theme_classic() + xlab("Substitutions location") + ylab("") +
  #   geom_vline(xintercept = cut_site, color="red") +
  #   scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

  #figure <- ggarrange(ggarrange(a, b, d, e,  ncol = 2, nrow = 2), subs_plot, ncol = 2)
  #figure <- ggarrange(a, b, d, e, ncol = 2, nrow = 2)
  #return(figure)
#}

#######
#### Reverse complement
#######
strComp=function(X){
  return(c2s(rev(comp(s2c(X)))))
}

######
### Cut site
######
get_cutSite <- function(gRNA_seq, reference, rel_cut){
  ### gRNA seq to upper case
  gRNA_seq <- toupper(gRNA_seq)
  ### Check orientation of gRNA in reference sequence and cut site
  rvComp_seq <- strComp(as.character(sread(reference)[[1]]))
  align <- pairwiseAlignment(toupper(gRNA_seq), toupper(as.character(sread(reference)[[1]])), type="global-local")
  alignRevComp <- pairwiseAlignment(toupper(gRNA_seq), toupper(rvComp_seq), type="global-local")
  
  if (score(align) > score(alignRevComp)){
    if (rel_cut > 0){
      cut_site <- start(subject(align))+rel_cut-1
    } else {
      cut_site <- start(subject(align)) + ( nchar(gRNA_seq) + rel_cut ) -1
    }
  } else {
    if (rel_cut > 0){
      cut_site <- nchar(as.character(sread(reference)[[1]])) - end(subject(alignRevComp)) + (nchar(gRNA_seq)-rel_cut)
    } else {
      cut_site <- nchar(as.character(sread(reference)[[1]])) - end(subject(alignRevComp)) - rel_cut
    }
  }
  return(cut_site)
}


#######
#### Empty Plots
#######
empty_plot <- function(title = NULL){
  p <- plotly_empty(type = "scatter", mode = "markers") %>%
    config(
      displayModeBar = FALSE
    ) %>%
    layout(
      title = list(
        text = title,
        yref = "paper",
        y = 0.5
      )
    )
  return(p)
} 

######################
##### Alleles percentages plot
######################

indel_count_seq <- function(indel_data){
  # Get table with percentages of indels among all genotypes
  indel_count <- plyr::count(indel_data, vars = c("Start", "Length", "Modification"))
  indel_count <- indel_count[order(indel_count$freq, decreasing = TRUE),]
  indel_count$Perc <- (indel_count$freq/sum(indel_count$freq)) * 100
  return(indel_count)
}

get_sequences <- function(indels, ref_seq){
  edited_sequences <- c()
  seqs <- lapply(c(1:dim(indels)[1]), function(num) {
    ref_seq_test <- as.character(sread(ref_seq))
    if (indels$Modification[[num]] == "del"){
      substr(ref_seq_test, indels$Start[[num]], indels$Start[[num]]+indels$Length[[num]]-1) <- paste0(rep("-", indels$Length[[num]]), collapse = "")
      edited_sequences <- c(edited_sequences, ref_seq_test)
    } else {
      seq <- paste(substring(ref_seq_test, 1, indels$Start[[num]]), paste0(rep("N", indels$Length[[num]]), collapse = ""), substring(ref_seq_test, indels$Start[[num]]+1, nchar(ref_seq_test)), sep = "")
      edited_sequences <- c(edited_sequences, seq)
    }
  })
}

######################
##### Substitutions at -+25 of the gRNA plot
######################
subs_plot <- function(subsperc, gRNA_seq, cut_site){
  ### Get the gRNA region
  pre_cut_site <- cut_site - (nchar(gRNA_seq)-3)  ## insted of 17 to allow different gRNA lengths instead of only 20
  post_cut_site <- cut_site + 3 + 1
  ### gRNA nucleotides to use them in the axis
  ref_nt <- stringr::str_split(gRNA_seq, "")[[1]]
  ### Get the plot
  plot <- ggplot(data=subsperc %>% filter(pos > pre_cut_site - 25) %>% filter(pos < post_cut_site + 25), aes(x=ordered(pos), y=percentage, fill=nucleotide, alpha=ifelse(percentage<95,1,0))) +
      geom_bar(stat="identity") + theme_classic() + scale_fill_manual("Nuclotides", values = c("A" = "#109648", "C" = "#86b7ed", "G" = "#f7b32b", "T" = "#d62839", "-" = "#f2f2f2")) +
      scale_x_discrete(breaks = (pre_cut_site+1):(post_cut_site-1), labels = ref_nt) + xlab(NULL) + scale_alpha(guide = 'none') + ylab('nt (%)') +
      geom_text(aes(label = ifelse(percentage<50 & percentage>5,percentage,"")),
                hjust = 0, vjust = 1.5, angle = 90, nudge_x = -.5,
                size = 2.5, family = "Public Sans") + 
      theme(text = element_text(size = 11, family = "Public Sans"), legend.position = "bottom")
  return(plot)
}

############
### Substitutions logo plot
############
subs_logo_plot <- function(subsperc, gRNA_seq, cut_site){
  ### Positions related to the cut site
  pre_cut_site <- cut_site - (nchar(gRNA_seq)-3)
  post_cut_site <- cut_site + 6 + 1
  #
  ### From percentages table to matrix
  fper_filtered <- subsperc %>% filter(pos > pre_cut_site) %>% filter(pos < post_cut_site) %>% filter(nucleotide != "-")
  nuc_pos <- pivot_wider(fper_filtered, names_from = pos, values_from = percentage)
  nuc_pos[is.na(nuc_pos)] <- 0
  nuc_pos[nuc_pos == 100] <- 0
  nuc_pos_matrix <- data.matrix(nuc_pos[,-1])
  nuc_pos_matrix <- nuc_pos_matrix+0.00001 ### This is done to aboid problemes in case of a matrix full of 0s
  rownames(nuc_pos_matrix) <- nuc_pos$nucleotide
  # nuc_pos_matrix <- nuc_pos_matrix+0.00001
  #
  ### gRNA nucleotides to use them in the axis
  ref_nt <- c(str_split(gRNA_seq, "")[[1]], "N", "G", "G")
  logo <- ggseqlogo(nuc_pos_matrix, method='custom', seq_type='dna') + ylab('nt (%)')
  logo$scales$scales[[1]] <- scale_x_continuous(breaks= 1:length(ref_nt),labels=ref_nt)
  logo_plot <- logo + annotate('rect', xmin = length(ref_nt)-2.5, xmax = length(ref_nt)+0.5, ymin = -5, ymax = 105, alpha = .1, col='black', fill='yellow')
  return(logo_plot)
}

##################
### Top variants plots
##################

### Logo plot of top indels 
get_logo_top_vars <- function(selected_reference, change_len, pattern_indel, modification, str_pos, cut_site){
  # This function is used to generate the logo of the alleles that apperar most times in data.
  all_pre_nt <- list(A = "", C = "", T = "", G = "")
  for (j in c("A", "C", "G", "T")){
    per_nt <- lapply(1:length(selected_reference), { 
      function(i){
        if(i < 11 || i > 10+l || modification == "ins"){
          if(selected_reference[i] == j){
            return(1-0.00001)
          } else {
            return(0.00001)
          }
        } else {
          return(0.00001)
        }
      }
    }
    )
    all_pre_nt[[j]] <- per_nt
  }
  ## Data frame with the probability of each position
  df<-NULL
  df <- rbind(df, unlist(all_pre_nt$A))
  df <- rbind(df, unlist(all_pre_nt$C))
  df <- rbind(df, unlist(all_pre_nt$G))
  df <- rbind(df, unlist(all_pre_nt$T))
  pos_matrix <- data.matrix(df)
  rownames(pos_matrix) <- c("A", "C", "G", "T")
  if(modification == "ins"){
    cut_correction <- change_len
  } else {
    cut_correction <- 0
  }
  
  logo <- ggseqlogo(pos_matrix, method='custom', seq_type='dna') + ylab('') +
    theme(plot.title = element_text(size = 12 , face = "bold" , hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle(paste0(l, "nt ", mod)) +
    annotate('segment', x = cut_site-str_pos+11+0.5+cut_correction , xend=cut_site-str_pos+11+0.5+cut_correction, y=-0.08, yend=1.08, size=1, color="red")
  
  if (is.na(pattern_indel) && modification != "ins"){
    return(logo)
  } else if (modification == "ins"){
    logo_plot <- logo + annotate('rect', xmin = 11-0.5, xmax = 10+change_len+0.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow')
    return(logo_plot)
  } else if (pattern_indel == "NHEJ"){
    return(logo)
  } else {
    mh_len <- nchar(pattern_indel)
    logo_plot <- logo + annotate('rect', xmin = 11+change_len-0.5, xmax = 10+change_len+mh_len+0.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow')
    return(logo_plot)
  }
}



##################
### Get all plots
##################

### Data
sample_name <- args[1]
data <- read.csv(args[2])
ref_seq <- readFasta(args[3])
gR <- args[4]
#subs_info <- args[5]
subs_plup <- read.csv(args[5], row.names = 1)
rel_cut_site <- as.numeric(args[6])

#load(subs_info)
#empty data troubleshooter, notice that	3 is choose to avoid counting small empty df product of	read.csv function
a<-as.character(data)
sampleName<-gsub('.{11}$', '', a)
checkFaulty<-grep('Faulty', sampleName)
a<-as.character(args[2])
a<-gsub('.{11}$', '', a)
checkEmpty<-grep('Empt', a)
checkFaulty<-grep('Badl', a)
if (dim(data)[2]>3 && length(checkFaulty) == 0 && length(checkEmpty) == 0){ ### Substitutions percentages plot
  cut_site <- get_cutSite(gR, ref_seq, rel_cut_site)
  if (dim(subs_plup)[1] == 0){
    system(paste0("touch ", sample_name, "_subs-perc_plot.png"))
  } else {
    sp <- subs_plot(subs_plup, gR, cut_site)
    ggsave(paste0(sample_name, "_subs-perc_plot.png"), width = 12, height = 1.5)
    ### Substitutions logo plot
    sp_logo <- subs_logo_plot(subs_plup, gR, cut_site)
    ggsave(paste0(sample_name, "_subs-perc_plot_LOGO.png"),  width = 12, height = 1.5)
  }
  
  ### Counts plot
  #GET THE HTMLS
  figure_counts <- plotIndels_gg(indels_data = data, cut =  cut_site)
  htmlwidgets::saveWidget(as_widget(figure_counts[[1]]), paste0(sample_name,"_Deletions.html"))
  htmlwidgets::saveWidget(as_widget(figure_counts[[2]]), paste0(sample_name,"_Insertions.html"))
#  htmlwidgets::saveWidget(as_widget(figure_counts[[3]]), paste0(sample_name,"_accumulative.html"))
    
  ### Alleles plot (deletions)
  # indel_count <- indel_count_seq(data)
  # dels_count <- indel_count %>% filter(Modification == "del") ## Just deletions
  # if (dim(dels_count)[1] != 0){
  #   dels_count$Start <- as.numeric(as.character(dels_count$Start))
  #   dels_count$Length <- as.numeric(as.character(dels_count$Length))
  #   seqs <- get_sequences(dels_count, ref_seq)
  #   
  #   # Write fasta with the sequences and use it to get the alleles plot
  #   if (length(seqs) < 10){
  #     rep_seqs <- seqs[1:length(seqs)]
  #   } else {
  #     rep_seqs <- seqs[1:10]
  #   }
  #   write.fasta(sequences = rep_seqs, names = paste0(round(dels_count$Perc, 2), "_", dels_count$Modification, dels_count$Length, "_", seq(1, dim(dels_count)[1]))[1:10], file.out = paste0(sample_name, "_top10.fa"))
  #   nt_sequences <- paste0(sample_name, "_top10.fa")
  #   figure_alleles <- ggmsa(nt_sequences, cut_site-40, cut_site+40, color = "Chemistry_NT", seq_name = TRUE)
  #   ggsave(paste0(sample_name, "_delAlleles_plot.png"))
  #   
  # } else {
  #   system(paste0("touch ",	sample_name, "_delAlleles_plot.png"))
  # }
  
  ######### Low variability or top variants plot
  wt <- data$wt_reads[1] ### 7299
  templata_based <- data$t_reads[1] ### 0 
  total_char <- wt + templata_based + dim(data)[1]
  
  delCols_indels <- data %>% select(c(Modification, Start, Length, freq))
  uniq_indels <- delCols_indels %>% distinct()
  uniq_indels_sorted <- uniq_indels[order(uniq_indels$freq, decreasing = TRUE),]
  
  # Check if there are enogh indels to have 5 top
  if(dim(delCols_indels)[1] < 4 ){
    num_top <- dim(delCols_indels)[1] + 1
  } else { num_top <- 5 }
  
  candidates <- rbind(uniq_indels_sorted[1:num_top-1,], c("wt", 0, 0, wt), c("template-based", 0, 0, templata_based))
  candidates_sorted <- candidates[order(candidates$freq, decreasing = TRUE),]
  top_5 <- candidates_sorted[1:num_top,]
  
  top5_names <- lapply(c(1:dim(top_5)[1]), 
                       function(i){
                         if( top_5[i,]$Modification %in% c("del", "ins") ){
                           return(paste0(top_5[i,]$Start, "_", top_5[i,]$Modification, top_5[i,]$Length))
                         } else {
                           return(top_5[i,]$Modification)
                         }
                       }
  )
  
  reads_classes <- c("Other alleles", "Top alleles", unlist(top5_names))
  reads_counts <- c(total_char - sum(as.numeric(top_5$freq)), sum(as.numeric(top_5$freq)), as.numeric(top_5$freq))
  reads_summary <- data.frame(classes = unlist(reads_classes), counts = unlist(reads_counts))
  reads_summary$parents =  c("", "", rep("Top alleles", num_top)) 
  fig <- plot_ly(reads_summary,
                 labels = ~classes,
                 parents = ~parents,
                 values = ~counts,
                 type = 'sunburst',
                 branchvalues = 'total',
                 textinfo = "label+percent entry",
                 marker = list(colors = c("#f2f2f2", "#1cf453", "#ff632c", "#9394f7", "#f71e71", "#19f9ce", "#fff933")))
  
  htmlwidgets::saveWidget(as_widget(fig), paste0(sample_name,"_top.html"))
  
  ### Logo plot 
  all_each_logo <- list()
  list_num <- 1
  sel_top <- top_5 %>% filter(Modification %in% c("del", "ins"))
  if (dim(sel_top)[1] > 0){
    for (i in c(1:dim(sel_top)[1])){
      selfil_1 <- data %>% filter(Modification ==  sel_top[i,]$Modification) %>% filter(Start ==  sel_top[i,]$Start) %>% filter(Length ==  sel_top[i,]$Length)
      s <- selfil_1[1,]$Start
      l <- selfil_1[1,]$Length
      if( s > 12  && (nchar(as.character(sread(ref_seq)[[1]])) > (s+l+9) )){
        ref_splited <- c(str_split(toupper(as.character(sread(ref_seq)[[1]])), "")[[1]])
        sel_ref <- ref_splited[(s-10):(s+l+9)]
        mod <- selfil_1[1,]$Modification
        if (mod == "ins"){
          sel_ref <- c(sel_ref[1:10], c(str_split(toupper(selfil_1[1,]$ins_nt), "")[[1]]), sel_ref[11:length(sel_ref)])
        }
        p <- selfil_1[1,]$patterns
        each_logo <- get_logo_top_vars(sel_ref, l, p, mod, s, cut_site)
        all_each_logo[[list_num]] <- each_logo
        list_num <- list_num + 1
      }
    }
  }
  
  if (length(all_each_logo) > 0){
    plot_grid(plotlist =  all_each_logo, ncol = 1)
    ggsave(paste0(sample_name, "_top-alleles_LOGO.png"))
  } else {
    system(paste0("touch ",	sample_name, "_top-alleles_LOGO.png"))
  }
  
} else {
  fig<-empty_plot("No alignments were produced.
  Please check your files and references")
  htmlwidgets::saveWidget(as_widget(fig), paste0(sample_name,"_Deletions.html"))
  htmlwidgets::saveWidget(as_widget(fig), paste0(sample_name,"_Insertions.html"))
  #htmlwidgets::saveWidget(as_widget(fig), paste0(sample_name,"_accumulative.html"))
  # system(paste0("touch ",	sample_name, "_delAlleles_plot.png"))
  
  #Elements for dynamic table IF THERE ARE NO ALIGNMENTS
  empty_list = list()
  exportJson <- jsonlite::toJSON(empty_list)
  write(exportJson, paste(sample_name,"_length.json",sep=""))
  write(exportJson, paste(sample_name,"_Total_reads.json",sep=""))
  write(exportJson, paste(sample_name,"_metadels.json",sep=""))
  write(exportJson, paste(sample_name,"_metains.json",sep=""))
  
  system(paste0("touch ",	sample_name, "_top.html"))
  system(paste0("touch ",	sample_name, "_top-alleles_LOGO.png"))
  system(paste0("touch ",	sample_name, "_counts_plot.png"))
  system(paste0("touch ",	sample_name, "_subs-perc_plot_LOGO.png"))
  system(paste0("touch ",	sample_name, "_subs-perc_plot.png"))
}
