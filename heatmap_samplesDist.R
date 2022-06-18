#!/usr/bin/env Rscript

library(ggplotify)
library(pheatmap)
library(philentropy)
library(heatmaply)

df_lenghts <- data.frame()
df_dist <- data.frame()
cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

distances_samples <- function(csv1, sample_names){
  file1 <- read.csv(csv1)
  splitted <- split(file1, file1$sample)
  # names(splitted) <- names(splitted)[names(splitted) != "sample"]
  samples_to_rm <- c()
  for (i in 1:length(splitted)){
    if (names(splitted)[i] != "sample"){
      file1 <- splitted[[i]]
      if (nrow(file1) > 1){
        file1<-subset(file1,file1$Modification=="del")
        
        ####positions transformation
        file1$Dist <- as.numeric(file1$Start) - as.numeric(file1$cut_site)
        file1$Dist <- sort(abs(file1$Dist))
        ############################
        
        file1$Length <- sort(file1$Length)
        
        file1_l=as.data.frame((table(file1$Length)))
        file1_s=as.data.frame((table(file1$Dist)))
        if (nrow(file1_l) > 0){
          colnames(file1_l) <- c("x", "Freq")
          colnames(file1_s) <- c("x", "Freq")
          df_lenghts <- data.frame(cbind.fill(df_lenghts, file1_l$Freq))
          df_dist <- data.frame(cbind.fill(file1_s$Freq, df_dist))
          
          for (j in 1:ncol(df_lenghts)){
            colnames(df_lenghts)[j] <- names(splitted)[j+1]
            colnames(df_dist)[j] <- names(splitted) [j+1]
          }      
        }
      } else {
        samples_to_rm <-c(samples_to_rm, i)
      } 
    } else {
      samples_to_rm <-c(samples_to_rm, i)
    }
  }

  #### Select final sample names
  pre_sel_names <- names(splitted)[-samples_to_rm]
  pre_sel_names_2 <- unlist(lapply(1:length(pre_sel_names), function(x){strsplit(pre_sel_names[x], "sample_")[[1]][2]}))
  pre_sel_names_3 <- unlist(lapply(1:length(pre_sel_names_2), 
    function(x){
      sname <- strsplit(pre_sel_names_2[x], "[.]")[[1]]
      if(length(sname) == 2){
        return(sname[2])
      } else {
        return(pre_sel_names_2[x])
      }
    }))

  
  df <- df_lenghts
  df[is.na(df)] <- 0
  df <- t(as.matrix(df))
  rownames(df) <- sample_names[(as.numeric(pre_sel_names_3)+1)]
  
  df_s <- df_dist
  df_s[is.na(df_s)] <- 0
  df_s <- t(as.matrix(df_s))
  rownames(df_s) <-  sample_names[(as.numeric(pre_sel_names_3)+1)]
  
  ### Kullback-leibler Length
  # d1_kl <- philentropy::distance(prop.table(df), method = "kullback-leibler")
  # colnames(d1_kl) <- c(names(df_lenghts))
  # rownames(d1_kl) <- c(names(df_lenghts))
  # d1_kll_plot <- heatmaply(d1_kl, colors = RdBu, main = 'Kullback-leibler Length', file = "kl-length.png")
  # 
  # ### Kullback-leibler Position
  # d1_kl_s <- philentropy::distance(prop.table(df_s), method = "kullback-leibler")
  # colnames(d1_kl_s) <-  c(names(df_dist))
  # rownames(d1_kl_s) <-  c(names(df_dist))
  # d1_klp_plot <- heatmaply(d1_kl_s, colors = RdBu, main = 'Kullback-leibler Position', file = "kl-position.png")
  
  ### Jensen Length
  d1_jen <- philentropy::distance(prop.table(df),method = "jensen_difference")
  colnames(d1_jen) <- sample_names[(as.numeric(pre_sel_names_3)+1)]
  rownames(d1_jen) <- sample_names[(as.numeric(pre_sel_names_3)+1)]
  d1_kl_plot_len <- heatmaply(d1_jen, colors = RdBu)
  
  ### Jensen Position
  d1_jen_s <- philentropy::distance(prop.table(df_s),method = "jensen_difference")
  colnames(d1_jen_s) <- sample_names[(as.numeric(pre_sel_names_3)+1)]
  rownames(d1_jen_s) <- sample_names[(as.numeric(pre_sel_names_3)+1)]
  d1_kl_plot <- heatmaply(d1_jen_s, colors = RdBu)
  
  #### Put both plots together
  jensen_plot <- subplot(d1_kl_plot, d1_kl_plot_len, margin = .05)
  htmlwidgets::saveWidget(as_widget(jensen_plot), "jensen_pos-len.html")
  
  return()
}

args = commandArgs(trailingOnly=TRUE)

csv1 <- args[1]
a <- read.csv(args[2], sep = "'", header = F)

sample_names <- c()
for (i in c(2:(dim(a)[1]-1)) ){
  sample_names <- c(sample_names, strsplit(a[i,2], "_")[[1]][1])
}

# sample_names <- c("AV-TnpB-1","AV-TnpB-5A","AV-TnpB-5B","AV-TnpB-6","AV-TnpB-7","AV-TnpB-8","AV-TnpB-9","AV-TnpB-10","AV-TnpB-11")
distances_samples(csv1, sample_names)
