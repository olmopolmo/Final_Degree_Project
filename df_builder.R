setwd("/Users/Asus/Desktop/BIOINFO/PRBB/Training_Dataset/RandomForest_HDR/HDR_results/HEK4_PE/")

samples <- read.csv(file = "Results/all-samples_indels.csv")
s <- unique(samples$sample)
s <- s[-2]
gRNA <- "GGCACTGCGGCTGGAGGTGG" #1
#gRNA <- "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGGACCGAGTCGGTCC" #2 

#################
### HDR       ###
#################
final_df <- data.frame()
for(chr in 1:length(s)){
  #################
  # Load the data #
  #################
  sample_id <- s[chr]
  sequence_name <- sample_id
  gRNA_sequence = toupper(gRNA) 
  template_ids <- read.csv(paste("template-based_reads-id/",sample_id,"_template-reads.csv",sep=""))
  indels_ids <-  read.csv(paste("Results/",sample_id,"_indels.csv", sep = ""))
  ref_fasta = paste("Data/",sample_id,"_reference.fasta", sep = "")
  #build template
  template = paste("Data/",sample_id,"_template.fasta", sep = "")
  
  cut_site = get_cutSite(gRNA_seq = gRNA_sequence,reference = readFasta(ref_fasta))
  
  # Bam of the Reference vs reads
  bamPath <- paste("Data/",sample_id,".sorted.bam", sep = "")
  bamFile <- BamFile(bamPath)
  bam <- scanBam(bamFile)
  RefvsRead_bam <- bam[[1]]
  temp_seq <- sread(readFasta(template))[[1]]
  ref_seq <- sread(readFasta(ref_fasta))[[1]]
  ref <- toupper(read.fasta(ref_fasta,as.string = T))
  ### Generate new reference (reference with change done by the template)
  newRef_Results <- newReference(template_sequence = temp_seq, reference_sequence = ref_seq)
  NewRef <- newRef_Results[1]
  
  a1 <-  paste("bash -c 'minimap2 -d reference.mmi Results/",sample_id,"_reference-correctOrient.fasta'", sep = "")
  a2 <- paste("bash -c '", "minimap2 -a -A 29 -B 17 -O 25 -E 2 reference.mmi NewRef.fasta ", "> template.sam'", sep = "")
  b2 <- paste("bash -c 'samtools view -bS template.sam -o template.bam'")
  system(a1)
  system(a2)
  system(b2)
  
  # Bam of Reference vs NewReference
  bamPath <- "template.bam" #NEED TO BE CHANGED DEPENDING ON THE DATA
  bamFile <- BamFile(bamPath)
  bam <- scanBam(bamFile)
  RefvsNew_bam <- bam[[1]]
  
  ##################################################################################################
  # 1 Get alignments data and Trim masked nucleotides from alignments and remove not aligned reads #
  ##################################################################################################
  
  RefvsRead <- data.frame(pos=RefvsRead_bam$pos,cigar=RefvsRead_bam$cigar, id=RefvsRead_bam$seq, qname=RefvsRead_bam$qname)
  
  corrected_cigar <- RefvsRead
  corrected_cigar$cigar <- lapply(RefvsRead$cigar,ignore_masked_p)
  corrected_cigar$cigar <- unlist(corrected_cigar$cigar, recursive = TRUE, use.names = TRUE)
  corrected_cigar <- corrected_cigar[-which(is.na(corrected_cigar$cigar)),]
  corrected_cigar$pos  <- replace(corrected_cigar$pos, is.na(corrected_cigar$pos),0)
  RefvsRead_bam <- corrected_cigar
  
  RefvsNew <- data.frame(pos=RefvsNew_bam$pos,cigar=RefvsNew_bam$cigar, id=RefvsNew_bam$seq)
  
  corrected_cigar <- RefvsNew
  corrected_cigar$cigar <- lapply(RefvsNew$cigar,ignore_masked_p)
  corrected_cigar$cigar <- unlist(corrected_cigar$cigar, recursive = TRUE, use.names = TRUE)
  corrected_cigar$cigar <- replace(corrected_cigar$cigar, is.na(corrected_cigar$cigar), "")
  corrected_cigar$pos  <- replace(corrected_cigar$pos, is.na(corrected_cigar$pos),0)
  RefvsNew_bam <- corrected_cigar

  
  #######################################################################
  # 3 Obtain the count for all the modifications that occur in the data #
  #######################################################################
  templates_qname <- template_ids$x
  
  indels_df <- data.frame(indels_ids$Ids, indels_ids$Modification)
  ins_df <- indels_df[which(indels_df$indels_ids.Modification=="ins"),]
  del_df <- indels_df[which(indels_df$indels_ids.Modification=="del"),]
  
  groups_df <- RefvsRead_bam
  groups_df$group <- NA
  cigar_lengths=explodeCigarOpLengths(groups_df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
  cigar_types=explodeCigarOps(groups_df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
  groups_df$modif <- NA
  
  template_based_qnames<- intersect(groups_df$qname,templates_qname)
  dels_qnames <- intersect(groups_df$qname,del_df$indels_ids.Ids)
  ins_qnames <- intersect(groups_df$qname,ins_df$indels_ids.Ids)
  
  #PARAMETERS FOR THE DATAFRAME
  DF <- data.frame(groups_df$cigar)
  colnames(DF) <- c("cigar")
  #Modification
  DF$modif <- NA
  #Group (A priori)
  DF$group <- NA
  #Distance to cut site
  DF$Dcutsite <- 0
  #Modification Size
  DF$size <- 0
  #Base change in case of substitution
  #DF$BaseModif <- "-"
  
  for (i in 1:length(groups_df$pos)){
    #GROUP
    if(groups_df$qname[i] %in% templates_qname){
      DF$group[i] <- "template"
    }else if(groups_df$qname[i] %in% ins_qnames){
      DF$group[i] <- "ins"
    }else if(groups_df$qname[i] %in% dels_qnames){
      DF$group[i] <- "del"
    }else{
      DF$group[i] <- "wt"
    }
    #MODIFICATION + Size + Distance to cut site
    if(length(cigar_types[[i]])==1){
      DF$modif[i] <- "wt"
      DF$Dcutsite[i] <- cigar_lengths[[i]][1]
      DF$size[i] <- 0
      
    }else if(length(cigar_types[[i]])<=3){
      if(cigar_types[[i]][2]=="I"){
        DF$modif[i] <- "ins"
        DF$Dcutsite[i] <- cigar_lengths[[i]][1] - cut_site
        DF$size[i] <- cigar_lengths[[i]][2]
        
      }else if(cigar_types[[i]][2]=="D"){
        DF$modif[i] <- "del"
        DF$Dcutsite[i] <- cigar_lengths[[i]][1] - cut_site
        DF$size[i] <- cigar_lengths[[i]][2]
      }
    }else{
      DF$modif[i] <- "delin"
      DF$Dcutsite[i] <- cigar_lengths[[i]][1] - cut_site
      ind <- which(cigar_types[[i]]=="M")
      DF$size[i] <- mean(cigar_lengths[[i]][-ind])
    }
    
    
  }
}
df_filtered <- as.data.frame(table(DF$cigar,DF$group,DF$modif,DF$Dcutsite, DF$size))
colnames(df_filtered) <- c("cigar", "Group", "Modif", "cut_site_dist", "size", "Freq")
df_filtered <- df_filtered[which(df_filtered$Freq>20),]
final_df <- df_filtered
colnames(final_df)
final_df$cigar <- NULL
final_df <- droplevels(final_df) #SOLUTION TO Error in randomForest.default(m, y, ...) : Can't have empty classes in y.

#data partition
d_1 <- which(final_df$group=="template")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.8 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_1 <- final_df[train_ind, ] #separate between train and test
test_1 <- final_df[setdiff(d_1,train_ind),]

d_2 <- which(final_df$group=="wt")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.8 * length(d_2)) #get the sample size desired in this case 0.7
train_ind <- sample(d_2, size = smp_size) #sample the amount specified before
train_2 <- final_df[train_ind, ] #separate between train and test
test_2 <- final_df[setdiff(d_2,train_ind),]

d_3 <- which(final_df$group=="del")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.8 * length(d_3)) #get the sample size desired in this case 0.7
train_ind <- sample(d_3, size = smp_size) #sample the amount specified before
train_3 <- final_df[train_ind, ] #separate between train and test
test_3 <- final_df[setdiff(d_3,train_ind),]

d_4 <- which(final_df$group=="ins")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.8 * length(d_4)) #get the sample size desired in this case 0.7
train_ind <- sample(d_4, size = smp_size) #sample the amount specified before
train_4 <- final_df[train_ind, ] #separate between train and test
test_4 <- final_df[setdiff(d_4,train_ind),]

d_6 <- which(final_df$group=="delin")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.8 * length(d_6)) #get the sample size desired in this case 0.7
train_ind <- sample(d_6, size = smp_size) #sample the amount specified before
train_6 <- final_df[train_ind, ] #separate between train and test
test_6 <- final_df[setdiff(d_6,train_ind),]

train <- rbind(train_1,train_2,train_3,train_4,train_6)
test <- rbind(test_1,test_2,test_3,test_4,train_6)



# DECISION TREES
colnames(train)

# tree <- rpart(  Group ~  Freq + cut_site_dist + wt_similarity + size + GC_content + Modif + BaseModif
#                 , data = train, method='class') #run rpart() function to build my tree
# 
# fancyRpartPlot(tree) #ploting the tree
# tree.pred <- predict(tree, test, type = "class") #building th predict
# table(tree.pred, test$Group)

Rforest <- randomForest(formula = group ~ Dcutsite + size + modif 
                        , data = train, method='class') #run rpart() function to build my tree
Rforest.pred <- predict(Rforest, test, type = "class") #building th predict

