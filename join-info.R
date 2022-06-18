library(dplyr)
library(seqinr)

##################################################################################
################### PE DELETIONS + INSERTIONS + SUBSTITUTIONS ####################
##################################################################################


setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/")
fastq_info = read.csv("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/fastqs-info.csv")
pegs_info = read.csv("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/pegs-info.csv")

### Data from sra run selector
samples_ids = read.csv("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/srr-srx.csv")

fastq_info$pegRNA <- rep(1, dim(fastq_info)[1])       
fastq_info$spacer <- rep(1, dim(fastq_info)[1])                        
fastq_info$X3extension <- rep(1, dim(fastq_info)[1]) 
fastq_info$PBS.length <- rep(1, dim(fastq_info)[1]) 
fastq_info$RT.template.length <- rep(1, dim(fastq_info)[1]) 

names(samples_ids) <- c("Experiment.Accession", "Run")
fastq_info <- merge(fastq_info, samples_ids, by = "Experiment.Accession")

#### Data from figure 2
# Pegs info --> 470
# Fastqs --> 2957
# 
# fig4 --> 107 --> 453 fastqs
final_data <- data.frame()

fig4 <- pegs_info %>% filter(grepl("^4", Figure))



# FIG 4 A
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse = "")
  }
  
  if (peginfo[2] %in% c("4a")){
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title))  %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension 
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_a <- rbind(final_data_a, data_related)
    } else {
      final_data_a <- data_related
      first_row <- first_row + 1
    }
  }
}

# FIG 4 B
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o" || strsplit(peginfo[3],"")[[1]][5] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse = "")
  }
  
  if (peginfo[2] %in% c("4b")){
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title))  %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_b <- rbind(final_data_b, data_related)
    } else {
      final_data_b <- data_related
      first_row <- first_row + 1
    }
  }
}

# FIG 4 C
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o" || strsplit(peginfo[3],"")[[1]][5] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse =  "")
  }
  
  if (peginfo[2] %in% c("4c")){
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title))  %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_c <- rbind(final_data_c, data_related)
    } else {
      final_data_c <- data_related
      first_row <- first_row + 1
    }
  }
}

# FIG 4 H
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o" || strsplit(peginfo[3],"")[[1]][5] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse =  "")
  }
  
  if (peginfo[2] %in% c("4_h")){
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4_", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title)) %>% filter(grepl(peginfo[1], Experiment.Title)) %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_h <- rbind(final_data_h, data_related)
    } else {
      final_data_h <- data_related
      first_row <- first_row + 1
    }
  }
}

# FIG 4 G
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o" || strsplit(peginfo[3],"")[[1]][5] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse =  "")
  }
  
  if (peginfo[2] %in% c("4_g")){
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4_", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title)) %>% filter(grepl(peginfo[1], Experiment.Title)) %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_g <- rbind(final_data_g, data_related)
    } else {
      final_data_g <- data_related
      first_row <- first_row + 1
    }
  }
}


# FIG 4 D E F
first_row <- 0
pegs_not_found_in_fastqs <- c()
for (i in c(1:dim(fig4)[1])){
  peginfo <- strsplit(fig4[i,]$pegRNA, "_")[[1]]
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse =  "")
  }
  
  if (peginfo[2] %in% c("4a", "4b", "4c","4h","4g")){ #perguntale a maarta pq faltan tantos
    pegs_not_found_in_fastqs <- c(pegs_not_found_in_fastqs, fig4[i,])
  } else {
    data_related <- fastq_info %>% filter(grepl(paste0("Fig4_", strsplit(peginfo[2],"")[[1]][2]), Experiment.Title)) %>% filter(grepl(peginfo[1], Experiment.Title)) %>% filter(grepl(peginfo[3], Experiment.Title))
    data_related$pegRNA <- fig4[i,]$pegRNA
    data_related$gRNA <- fig4[i,]$spacer
    data_related$template <- fig4[i,]$X3â...extension
    data_related$PBS.length <- fig4[i,]$PBS.length
    data_related$RT.template.length <- fig4[i,]$RT.template.length
    if (first_row != 0){
      final_data_def <- rbind(final_data_def, data_related)
    } else {
      final_data_def <- data_related
      first_row <- first_row + 1
    }
  }
}

final_data <- rbind(final_data_a,final_data_b,final_data_c,final_data_def)
write.csv(final_data, "fig4_fastqs-info.csv")                  
setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/FIG4_PE_RUN/")

for (i in 1:length(final_data$Experiment.Accession)){
  write.fasta(as.string = T,sequences = final_data$template[i], open = "w",names = final_data$Run[i], file.out = paste(final_data$Run[i],".fasta",sep = ""))  
}


###################
# FIG 10 HDR + PE #
###################

fig10 <- pegs_info %>% filter(grepl("^EDF10", Figure))
fig10 <- fig10[-4,]

# FIG 10
first_row <- 0
pegs_not_found_in_fastqs <- c()

for (i in c(1:dim(fig10)[1])){
  peginfo <- strsplit(fig10[i,]$pegRNA, "_")[[1]]
  peginfo[3] <- gsub(" ","",peginfo[3])
  
  if(strsplit(peginfo[3],"")[[1]][4] %in% "o"){
    peginfo[3] <- paste(strsplit(peginfo[3],"to")[[1]],collapse ="")
    
  }
  
  if (peginfo[2] %in% c("ED10")){
    data_related <- fastq_info %>% filter(grepl("FigED10", fastq_info$Experiment.Title))  %>% filter(grepl(peginfo[3], Experiment.Title))  %>% filter(grepl(peginfo[1], Experiment.Title)) 
    data_related$pegRNA <- fig10[i,]$pegRNA
    data_related$gRNA <- fig10[i,]$spacer
    data_related$template <- fig10[i,]$X3â...extension
    data_related$PBS.length <- fig10[i,]$PBS.length
    data_related$RT.template.length <- fig10[i,]$RT.template.length
    if (first_row != 0){
      final_data_10 <- rbind(final_data_10, data_related)
    } else {
      final_data_10 <- data_related
      first_row <- first_row + 1
    }
  }
}

length(unique(final_data_10$Run))
length(final_data_10$Run)
write.csv(final_data_10, "fig10_fastqs-info.csv")                  
setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/FIG10_HDRPE_RUN/")

for (i in 1:length(final_data_10$Experiment.Accession)){
  write.fasta(as.string = T,sequences = final_data_10$template[i], open = "w",names = final_data_10$Run[i], file.out = paste(final_data_10$Run[i],".fasta",sep = ""))  
}


##################################################################################
################### HDR DELETIONS + INSERTIONS + SUBSTITUTIONS ###################
##################################################################################

info = read.csv("/Users/ASUS/Downloads/SraRunInfo (1).csv")
pegs_info = read.csv("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/fig2e_HDR.csv",header = T,sep = ";")

fastq_info = as.data.frame(cbind(info$Run,info$SampleName))
colnames(fastq_info) <- c("Run", "Sample")


fig2e <- fastq_info %>% filter(grepl("2e", Sample))


final_data <- data.frame()


fig2e$Key <- NA
for (i in 1:length(fig2e$Run)){
  Sample_elements <- strsplit(fig2e$Sample[i], "-")[[1]]
  if(Sample_elements[3]=="EvoPreQ1"){
    key <- paste(Sample_elements[2],Sample_elements[3],Sample_elements[5],substr(Sample_elements[6],1,5),sep = "-")  
  }else{
    key <- paste(Sample_elements[2],Sample_elements[5],substr(Sample_elements[6],1,5),sep = "-")
  }
  
  fig2e$Key[i] <- key
}


pegs_info$Door <- NA
for (i in 1:length(pegs_info$Spacer)){
  Sample_matches <- strsplit(pegs_info$ï..Sample_ID[i], "_")[[1]]
  if(length(Sample_matches)==6){
    door <- paste(Sample_matches[1],Sample_matches[6],Sample_matches[4],substr(Sample_matches[5],1,5),sep = "-")  
  }else{
    door <- paste(Sample_matches[1],Sample_matches[4],substr(Sample_matches[5],1,5),sep = "-")
  }
  
  pegs_info$Door[i] <- door
}

fig2e$X3_extension <- ""
fig2e$spacer <- ""
fig2e$template <- ""

for (i in 1:length(fig2e$Run)){
  fig2e$spacer[i] <- pegs_info$Spacer[which(fig2e$Key[i]==pegs_info$Door)[1]]
  fig2e$X3_extension[i] <- pegs_info$X3_extension[which(fig2e$Key[i]==pegs_info$Door)[1]]
  fig2e$template[i] <- fig2e$X3_extension[i]
  
}

setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/")
write.csv(fig2e, "fig2e_fastqs-info.csv")                  

for (i in 1:length(fig2e$Run)){
  write.fasta(as.string = T,sequences = fig2e$template[i], open = "w",names = fig2e$Run[i], file.out = paste(fig2e$Run[i],".fasta",sep = ""))  
}
