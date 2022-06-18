library(Rsamtools)
library(plyr)
library(Rtsne)
library(dplyr)
library(stringr)
library(ShortRead)
library(seqinr)
library(GenomicAlignments)
library(data.table)
library(biomaRt) 
library(ggplot2)
### HDR
library(DECIPHER)
library(Biostrings)
library(parallel)
library(minimap)
library(MASS)
library(tidyverse)
library(caret)
library(ggfortify)
library(rpart)
library(rpart.plot)
library(caret)
library(readr)
library(rattle)
library(randomForest)

library(DescTools)
#################
### Functions ###
#################

#######
#### Reverse complement
#######
strComp=function(X){
  return(c2s(rev(comp(s2c(X)))))
}

#######
### Trim masked nucleotides at the begining and the end of the alignment
#######
ignore_masked_p <- function(df){
  # Ignore S at the beginning or end of the sequence
  new_cigar <- c()
  for (cig in 1:length(df)) {
    startS <- str_replace_all(df[cig], "^\\d+S", "")
    endS <- str_replace_all(startS, "\\d+S$", "")
    new_cigar <- c(new_cigar, endS)
  }
  df <- new_cigar
  
  # Ignore H at the beginning or end of the sequence
  new_cigar <- c()
  for (cig in 1:length(df)) {
    startH <- str_replace_all(df[cig], "^\\d+H", "")
    endH <- str_replace_all(startH, "\\d+H$", "")
    new_cigar <- c(new_cigar, endH)
  }
  df <- new_cigar
  
  return(df)
}
#########
### Get new reference for template-based edition 
#########
newReference <- function(template_sequence, reference_sequence){
  ##########
  ##### This function generates a new reference with the change of the target integrated on it 
  ##### Its based in the assumption that a template is composed by two homology arms and a modification between them
  ##########
  
  # Align the template in both directions to check which is the correct one
  fw_alignment <- pairwiseAlignment(reference_sequence, template_sequence, type="local")
  revComp_alignment <- pairwiseAlignment(reference_sequence, reverseComplement(template_sequence), type="local")
  
  # Check if the alignment is better in forward or in reverse complement direction and save alignment and template in right 
  if (score(fw_alignment) > score(revComp_alignment)){
    rigth_seq <- template_sequence
    aln <- fw_alignment
  } else {
    rigth_seq <- reverseComplement(template_sequence)
    aln <- revComp_alignment
  }
  
  # Get template start (to see if left, right or both arms have aligned) and last position of template in reference
  temp_aln_start <- start(subject(aln))
  ref_aln_end <- start(pattern(aln)) + nmatch(aln) + nmismatch(aln)  #nchar(pattern(aln))
  
  # Which is the longer aligned arm? Template alignment starts at position 1 (--> left arm longer; aligns first)?
  if (temp_aln_start != 1 && start(pattern(aln)) != 1){ # when left arm has not aligned. We also check if the template starts before 
    left_arm <- substr(rigth_seq, 1, temp_aln_start-1)
    # The alignment starts after the change, we have to re-align the first part (short left arm) to see where in the reference starts the template
    left_start <- start(pattern(pairwiseAlignment(reference_sequence, left_arm, type="local"))) 
    new_reference <- paste0(substr(reference_sequence,1,left_start-1), rigth_seq, substr(reference_sequence,ref_aln_end,nchar(reference_sequence)))
    # Start and end location of the modification in the reference sequence
    change_start <- end(subject(pairwiseAlignment(reference_sequence, left_arm, type="local")))
    change_end <- temp_aln_start
  } else {
    # Has the whole template aligned or we have to look for the end of the right template?
    if (nchar(aln) >= (nchar(template_sequence) - temp_aln_start + 1) || (ref_aln_end-1 == nchar(reference_sequence)) ){ # Whole template has aligned --> corrected in case template starts before reference. Or we have arrived to the end of reference
      # Check if there is a deletion in the template. Otherwise, we have to take into account the deletion to jump its length in the reference
      if (!identical(width(insertion(aln))[[1]], integer(0))) {  ## There is a deletion in the template
        ref_aln_end <- ref_aln_end + sum(width(insertion(aln))[[1]])
        change_start <- start(insertion(aln)[[1]][1])
        change_end <- change_start + sum(width(insertion(aln))[[1]])
      } else if (!identical(width(deletion(aln))[[1]], integer(0))){ ## There is an insertion in the template
        change_start <- start(deletion(aln)[[1]][1])
        change_end <- change_start + 1 #since this is a deletion the end should be the the same position, but this causes errors later...
      } else { ## It's just a substitution
        change_start <- mismatchTable(aln)$SubjectStart[1]
        change_end <- mismatchTable(aln)$SubjectEnd[dim(mismatchTable(aln))[1]]+1
        if (is.na(change_start)){ ### This means that there is no change in the template!!
          change_start <- 0
          change_end <- 0
        }
      }
      new_reference <- paste0(substr(reference_sequence,1,start(pattern(aln))-1), rigth_seq, substr(reference_sequence,ref_aln_end,nchar(reference_sequence)))
    } else { # We have to look for the end of the right template
      right_arm <- substr(rigth_seq, nchar(aln)+1, length(rigth_seq))
      right_end <- end(pattern(pairwiseAlignment(reference_sequence, right_arm, type="local"))) 
      new_reference <- paste0(substr(reference_sequence,1,start(pattern(aln))-1), rigth_seq, substr(reference_sequence,right_end+1,nchar(reference_sequence)))
      change_start <- nchar(subject(aln)) 
      change_end <- nchar(aln)+1
    }
  }
  
  NewRef <- paste(unlist(new_reference),collapse='')
  write.fasta(new_reference,"NewRef", file.out = "NewRef.fasta")
  
  return(c(new_reference, change_start, change_end))
}

######
### Cut site
######
get_cutSite <- function(gRNA_seq, reference){
  ### From gRNA sequence get the position of the cut site in teh reference sequence
  # Split by gRNA to get non-edited expected region 
  gRNA_seq <- toupper(gRNA_seq)
  seq_splited <- strsplit(as.character(sread(reference)[[1]]), split=gRNA_seq, fixed=TRUE)
  cut_site <- nchar(seq_splited[[1]][1]) + 17
  if (length(seq_splited[[1]]) == 1){
    seq_splited <- strsplit(as.character(sread(reference)[[1]]), split=as.character(reverseComplement(DNAString(gRNA_seq))), fixed=TRUE)
    cut_site <- nchar(seq_splited[[1]][1]) + 3
  }
  return(cut_site)
}

#############
### Template-based quantification
#############
templateCount <- function(ori_ref_file, new_ref_temp, readsAlnInfo_collaps, sampleID, alnCompleteInfo){ #ref_fasta, "newRef.fa", collapsed_df, sample_id
  ########
  ### Function to get number of reads supporting template-based edits. The return will be a list with two values: count and kind of substitution
  ########
  #
  #### Align new reference with the change led by template with original reference
  #system(paste0("minimap2 -d ", sampleID, "_reference.mmi ", ori_ref_file)) #QUITAR EL COMENTARIO!!!!!!!!!!!!!!!!!!!!!!!
  #system(paste0("minimap2 -a -A 29 -B 17 -O 25 -E 2 ", sampleID, "_reference.mmi ", new_ref_temp, "> template.sam;")) #QUITAR EL COMENTARIO!!!!!!!!!!!!!!!!!!!!!!!
  #### Get cigar
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isDuplicate = FALSE, isSupplementaryAlignment = FALSE), ## Get just primary alignments
    what=c("cigar"))
  sam_file_path = "template.sam"
  bam=scanBam("/Users/Asus/Desktop/BIOINFO/PRBB_2nd/template.bam", param=param) 
  temp_cigar=bam[[1]]$cigar[1]
  #### Cigar length and type
  cigar_lengths=explodeCigarOpLengths(temp_cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
  cigar_types=explodeCigarOps(temp_cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
  if ("I" %in% cigar_types[[1]] || "D" %in% cigar_types[[1]]){
    #### When it is not a simple substitution; not all nucleotides are matches of mismatches (or clipped nucleotides!!)
    rowReads <- readsAlnInfo_collaps %>% filter(cigar == temp_cigar)
    if ("I" %in% cigar_types[[1]] && "D" %in% cigar_types[[1]]){
      tempCount <- c(rowReads$count, "delin")
    } else if ("I" %in% cigar_types[[1]]) {
      if(cigar_lengths[[1]][which(cigar_types[[1]] == "I")] %% 3 != 0){
        tempCount <- c(rowReads$count, "ins-out") 
      } else {
        tempCount <- c(rowReads$count, "ins-in") 
      }
    } else if ("D" %in% cigar_types[[1]]) {
      if(cigar_lengths[[1]][which(cigar_types[[1]] == "D")] %% 3 != 0){
        tempCount <- c(rowReads$count, "dels-out") 
      } else {
        tempCount <- c(rowReads$count, "dels-in") 
      }
    }
  } else { # there is a substitution. Alignment an pattern match of the substitution in that certain region
    ref_tempSeq <- sread(readFasta(new_ref_temp))[[1]]
    ref_oriSeq <- sread(readFasta(ori_ref_file))[[1]]
    aln <- pairwiseAlignment(ref_tempSeq, ref_oriSeq)
    if (nmismatch(aln) == 0){
      tempCount <- c(0, "no-changes")
    } else {
      s <- mismatchTable(aln)$PatternStart[1]
      e <- mismatchTable(aln)$PatternStart[length(mismatchTable(aln)$PatternStart)]
      diffinTemp <- substr(ref_tempSeq, s, e)
      matchChange <- lapply(alnCompleteInfo$seq, function(x){ nchar(x) == nchar(ref_tempSeq) && substr(x,s,e) == as.character(diffinTemp) })
      tempCount <- c(sum(unlist(matchChange)), "subs")
    }
  }
  return(tempCount)
}


#################
### Files #######
#################


#################################################################################
########################### TOY DATA SETS #######################################
#################################################################################

setwd("/Users/ASUS/Desktop/BIOINFO/PRBB/Training_Dataset/PAPER_dataset_creation/TRUE_DATASET/Output/")


samples <-  c("SRR10286462",
                       "SRR10286463",
                       "SRR10286464",
                       "SRR10286465",
                       "SRR10286469",
                       "SRR10286470",
                       "SRR10286471",
                       "SRR10286473",
                       "SRR10286474",
                       "SRR10286475",
                       "SRR10286476",
                       "SRR10286477",
                       "SRR10286478",
                       "SRR10286479",
                       "SRR10286480",
                       "SRR10286481",
                       "SRR10286482",
                       "SRR10286483",
                       "SRR10286484",
                       "SRR10286485",
                       "SRR10286486",
                       "SRR10286487",
                       "SRR10286488",
                       "SRR10286489",
                       "SRR10286490",
                       "SRR10286491",
                       "SRR10286492",
                       "SRR10286493",
                       "SRR10286494",
                       "SRR10286498",
                       "SRR10286499",
                       "SRR10286500",
                       "SRR10286501",
                       "SRR10286502",
                       "SRR10286503",
                       "SRR10286504",
                       "SRR10286505",
                       "SRR10286506",
                       "SRR10286507",
                       "SRR10286508",
                       "SRR10286509",
                       "SRR10286510",
                       "SRR10286511",
                       "SRR10286512",
                       "SRR10286513",
                       "SRR10286514",
                       "SRR10286515",
                       "SRR10286516",
                       "SRR10286517",
                       "SRR10286518",
                       "SRR10286519",
                       "SRR10286520",
                       "SRR10286521",
                       "SRR10286522",
                       "SRR10286523",
                       "SRR10286524",
                       "SRR10286525",
                       "SRR10286526",
                       "SRR10286527",
                       "SRR10286528",
                       "SRR10286529",
                       "SRR10286530",
                       "SRR10286534",
                       "SRR10286535",
                       "SRR10286537",
                       "SRR10286538",
                       "SRR10286539",
                       "SRR10286540",
                       "SRR10286541",
                       "SRR10286542",
                       "SRR10286543",
                       "SRR10286544",
                       "SRR10286545",
                       "SRR10286546",
                       "SRR10286547",
                       "SRR10286548",
                       "SRR10286549",
                       "SRR10286550",
                       "SRR10286551",
                       "SRR10286552",
                       "SRR10286553",
                       "SRR10286554",
                       "SRR10286555",
                       "SRR10286556",
                       "SRR10286557",
                       "SRR10286558",
                       "SRR10286559",
                       "SRR10286560",
                       "SRR10286561",
                       "SRR10286562",
                       "SRR10286563",
                       "SRR10286564",
                       "SRR10286565",
                       "SRR10286566",
                       "SRR10286567",
                       "SRR10286568",
                       "SRR10286569",
                       "SRR10286570",
                       "SRR10286571",
                       "SRR10286572",
                       "SRR10286573",
                       "SRR10286574",
                       "SRR10286575",
                       "SRR10286576",
                       "SRR10286577",
                       "SRR10286578",
                       "SRR10286579",
                       "SRR10286580",
                       "SRR10286581",
                       "SRR10286582",
                       "SRR10286583",
                       "SRR10286584",
                       "SRR10286585",
                       "SRR10286586",
                       "SRR10286587",
                       "SRR10286588",
                       "SRR10286589",
                       "SRR10286590",
                       "SRR10286591",
                       "SRR10286592",
                       "SRR10286593",
                       "SRR10286594",
                       "SRR10286595",
                       "SRR10286596",
                       "SRR10286597",
                       "SRR10286598",
                       "SRR10286599",
                       "SRR10286603",
                       "SRR10286604",
                       "SRR10286605",
                       "SRR10286606",
                       "SRR10286607",
                       "SRR10286608",
                       "SRR10286609",
                       "SRR10286610",
                       "SRR10286611",
                       "SRR10286612",
                       "SRR10286613",
                       "SRR10286614",
                       "SRR10286615",
                       "SRR10286616",
                       "SRR10286617",
                       "SRR10286618",
                       "SRR10286619",
                       "SRR10286620",
                       "SRR10286621",
                       "SRR10286622",
                       "SRR10286626",
                       "SRR10286627",
                       "SRR10286628",
                       "SRR10286629",
                       "SRR10286630",
                       "SRR10286631",
                       "SRR10286632",
                       "SRR10286633",
                       "SRR10286634",
                       "SRR10286635",
                       "SRR10286636",
                       "SRR10286637",
                       "SRR10286638",
                       "SRR10286639",
                       "SRR10286640",
                       "SRR10286641",
                       "SRR10286642",
                       "SRR10286643",
                       "SRR10286644",
                       "SRR10286645",
                       "SRR10286649",
                       "SRR10286650",
                       "SRR10286651",
                       "SRR10286652",
                       "SRR10286654",
                       "SRR10286655",
                       "SRR10286656",
                       "SRR10286660",
                       "SRR10286661",
                       "SRR10286662",
                       "SRR10286663",
                       "SRR10286664",
                       "SRR10286666",
                       "SRR10286667",
                       "SRR10286668",
                       "SRR10286669",
                       "SRR10286670",
                       "SRR10286671",
                       "SRR10286672",
                       "SRR10286673",
                       "SRR10286674",
                       "SRR10286675",
                       "SRR10286677",
                       "SRR10286678",
                       "SRR10286679",
                       "SRR10286683",
                       "SRR10286684",
                       "SRR10286685",
                       "SRR10286686",
                       "SRR10286687",
                       "SRR10286688",
                       "SRR10286689",
                       "SRR10286690",
                       "SRR10286691",
                       "SRR10286692",
                       "SRR10286693",
                       "SRR10286694",
                       "SRR10286695",
                       "SRR10286696",
                       "SRR10286697",
                       "SRR10286698",
                       "SRR10286699",
                       "SRR10286700",
                       "SRR10286701",
                       "SRR10286702",
                       "SRR10286706",
                       "SRR10286707",
                       "SRR10286708",
                       "SRR10286709",
                       "SRR10286710",
                       "SRR10286711",
                       "SRR10286712",
                       "SRR10286716",
                       "SRR10286717",
                       "SRR10286718",
                       "SRR10286719",
                       "SRR10286720",
                       "SRR10286721",
                       "SRR10286722",
                       "SRR10286723",
                       "SRR10286724",
                       "SRR10286725",
                       "SRR10286726",
                       "SRR10286727",
                       "SRR10286728",
                       "SRR10286729",
                       "SRR10286730",
                       "SRR10286731",
                       "SRR10286732",
                       "SRR10286733",
                       "SRR10286734",
                       "SRR10286735",
                       "SRR10286739",
                       "SRR10286740",
                       "SRR10286741",
                       "SRR10286742",
                       "SRR10286743",
                       "SRR10286744",
                       "SRR10286745",
                       "SRR10286746",
                       "SRR10286747",
                       "SRR10286748",
                       "SRR10286749",
                       "SRR10286750",
                       "SRR10286751",
                       "SRR10286752",
                       "SRR10286754",
                       "SRR10286755",
                       "SRR10286756",
                       "SRR10286757",
                       "SRR10286758",
                       "SRR10286762",
                       "SRR10286763",
                       "SRR10286765",
                       "SRR10286766",
                       "SRR10287764",
                       "SRR10287765",
                       "SRR10287766",
                       "SRR10287767",
                       "SRR10287768",
                       "SRR10287769",
                       "SRR10287770",
                       "SRR10287771",
                       "SRR10287772",
                       "SRR10287773",
                       "SRR10287774",
                       "SRR10287775",
                       "SRR10287776",
                       "SRR10287777",
                       "SRR10287778",
                       "SRR10287779",
                       "SRR10287780",
                       "SRR10287781",
                       "SRR10287782",
                       "SRR10287783",
                       "SRR10287784",
                       "SRR10287785",
                       "SRR10287786",
                       "SRR10287787",
                       "SRR10287788",
                       "SRR10287789",
                       "SRR10287790",
                       "SRR10287791",
                       "SRR10287792",
                       "SRR10287793",
                       "SRR10287794",
                       "SRR10287795",
                       "SRR10287796",
                       "SRR10287797",
                       "SRR10287798",
                       "SRR10287799",
                       "SRR10287800",
                       "SRR10287801",
                       "SRR10287802",
                       "SRR10287803",
                       "SRR10287804",
                       "SRR10287805",
                       "SRR10287806",
                       "SRR10287807",
                       "SRR10287808",
                       "SRR10287809",
                       "SRR10287810",
                       "SRR10287811",
                       "SRR10287812",
                       "SRR10287813",
                       "SRR10287814",
                       "SRR10287815",
                       "SRR10287816",
                       "SRR10287817",
                       "SRR10287818",
                       "SRR10287819",
                       "SRR10287820",
                       "SRR10287821",
                       "SRR10287822",
                       "SRR10287823",
                       "SRR10287824",
                       "SRR10287825",
                       "SRR10287826",
                       "SRR10287827",
                       "SRR10287828",
                       "SRR10287829",
                       "SRR10287830",
                       "SRR10287831",
                       "SRR10287832",
                       "SRR10287834",
                       "SRR10287835",
                       "SRR10287836",
                       "SRR10287837",
                       "SRR10287838",
                       "SRR10287839",
                       "SRR10287840",
                       "SRR10287841",
                       "SRR10287842",
                       "SRR10287843",
                       "SRR10287844",
                       "SRR10287845",
                       "SRR10287846",
                       "SRR10287847",
                       "SRR10287848",
                       "SRR10287849",
                       "SRR10287850",
                       "SRR10287851",
                       "SRR10287852",
                       "SRR10287853",
                       "SRR10287854",
                       "SRR10287855",
                       "SRR10287856",
                       "SRR10287857",
                       "SRR10287858",
                       "SRR10287859",
                       "SRR10287860",
                       "SRR10287861",
                       "SRR10287862",
                       "SRR10287863",
                       "SRR10287864",
                       "SRR10287865",
                       "SRR10287866",
                       "SRR10287867",
                       "SRR10287868",
                       "SRR10287869",
                       "SRR10287870",
                       "SRR10287871",
                       "SRR10287872",
                       "SRR10287873",
                       "SRR10287874",
                       "SRR10287875",
                       "SRR10287876",
                       "SRR10287877",
                       "SRR10287878",
                       "SRR10287879",
                       "SRR10287880",
                       "SRR10287881",
                       "SRR10287882",
                       "SRR10287883",
                       "SRR10287884",
                       "SRR10287885",
                       "SRR10287886",
                       "SRR10287887",
                       "SRR10287888",
                       "SRR10287889",
                       "SRR10287890",
                       "SRR10287891",
                       "SRR10287892",
                       "SRR10287893",
                       "SRR10287894",
                       "SRR10287895",
                       "SRR10287896",
                       "SRR10287897",
                       "SRR10287898",
                       "SRR10287899",
                       "SRR10287900",
                       "SRR10287901",
                       "SRR10287902",
                       "SRR10287903",
                       "SRR10287904",
                       "SRR10287905",
                       "SRR10287906",
                       "SRR10287907",
                       "SRR10287908",
                       "SRR10287909",
                       "SRR10287910",
                       "SRR10287911",
                       "SRR10287912",
                       "SRR10287913",
                       "SRR10287914",
                       "SRR10287915",
                       "SRR10287916",
                       "SRR10287917",
                       "SRR10287918",
                       "SRR10287919",
                       "SRR10287920",
                       "SRR10287921",
                       "SRR10287922",
                       "SRR10287923",
                       "SRR10287924",
                       "SRR10287925",
                       "SRR10287926",
                       "SRR10287927",
                       "SRR10287928",
                       "SRR10287929",
                       "SRR10287930",
                       "SRR10287931",
                       "SRR10287932",
                       "SRR10287933",
                       "SRR10287934",
                       "SRR10287935",
                       "SRR10287936",
                       "SRR10287937",
                       "SRR10287938",
                       "SRR10287939",
                       "SRR10287940",
                       "SRR10287941",
                       "SRR10287942",
                       "SRR10287943",
                       "SRR10287944",
                       "SRR10287945",
                       "SRR10287946",
                       "SRR10287947",
                       "SRR10287948",
                       "SRR10287949",
                       "SRR10287950",
                       "SRR10287951",
                       "SRR10287952",
                       "SRR10287953",
                       "SRR10287954",
                       "SRR10287955",
                       "SRR10287956",
                       "SRR10287957",
                       "SRR10287958",
                       "SRR10287959",
                       "SRR10287960",
                       "SRR10287961",
                       "SRR10287962",
                       "SRR10287963",
                       "SRR10287964",
                       "SRR10287965",
                       "SRR10287966",
                       "SRR10287967",
                       "SRR10287968",
                       "SRR10287969",
                       "SRR10287970",
                       "SRR10287971",
                       "SRR10287972",
                       "SRR10287973",
                       "SRR10287974",
                       "SRR10287975",
                       "SRR10287986",
                       "SRR10287987",
                       "SRR10287998",
                       "SRR10288009",
                       "SRR10288020",
                       "SRR10288031",
                       "SRR10288036",
                       "SRR10288037",
                       "SRR10288038",
                       "SRR10288039",
                       "SRR10288040",
                       "SRR10288041",
                       "SRR10288042",
                       "SRR10288043",
                       "SRR10288044",
                       "SRR10288045",
                       "SRR10288046",
                       "SRR10288047",
                       "SRR10288048",
                       "SRR10288049",
                       "SRR10288050",
                       "SRR10288051",
                       "SRR10288052",
                       "SRR10288053",
                       "SRR10288054",
                       "SRR10288055",
                       "SRR10288056",
                       "SRR10288057",
                       "SRR10288058",
                       "SRR10288059",
                       "SRR10288060",
                       "SRR10288061",
                       "SRR10288062",
                       "SRR10288063",
                       "SRR10288064",
                       "SRR10288065",
                       "SRR10288067",
                       "SRR10288068",
                       "SRR10288069",
                       "SRR10288070",
                       "SRR10288071",
                       "SRR10288072",
                       "SRR10288073",
                       "SRR10288074",
                       "SRR10288075",
                       "SRR10288076",
                       "SRR10288077",
                       "SRR10288078",
                       "SRR10288079",
                       "SRR10288080",
                       "SRR10288081",
                       "SRR10288082",
                       "SRR10288083",
                       "SRR10288084",
                       "SRR10288085",
                       "SRR10288086",
                       "SRR10288087",
                       "SRR10288088",
                       "SRR10288089",
                       "SRR10288090",
                       "SRR10288091",
                       "SRR10288092",
                       "SRR10288093",
                       "SRR10288094",
                       "SRR10288095",
                       "SRR10288096",
                       "SRR10288097",
                       "SRR10288098",
                       "SRR10288099",
                       "SRR10288100",
                       "SRR10288101",
                       "SRR10288102",
                       "SRR10288103",
                       "SRR10288104",
                       "SRR10288105",
                       "SRR10288106",
                       "SRR10288107",
                       "SRR10288108",
                       "SRR10288109",
                       "SRR10288110",
                       "SRR10288111",
                       "SRR17171765",
                       "SRR17171766",
                       "SRR17171767",
                       "SRR17171768",
                       "SRR17171770",
                       "SRR17171771",
                       "SRR17171772",
                       "SRR17171773",
                       "SRR17171774",
                       "SRR17171775",
                       "SRR17171776",
                       "SRR17171777",
                       "SRR17171778",
                       "SRR17171779",
                       "SRR17171782",
                       "SRR17171783",
                       "SRR17171784",
                       "SRR17171785",
                       "SRR17171786",
                       "SRR17171787",
                       "SRR17171788",
                       "SRR17171789",
                       "SRR17171790",
                       "SRR17171791")



#################################################################################
########################### REAL DATA SETS ######################################
#################################################################################


fig2<-read.csv("../../fig2e_fastqs-info.csv")
fig4<-read.csv("../../fig4_fastqs-info.csv")
fig10<-read.csv("../../fig10_fastqs-info.csv")

Run<-c(fig2$Run,fig4$Run,fig10$Run)
gRNA<-c(toupper(fig2$spacer),toupper(fig4$gRNA),toupper(fig10$gRNA))
templates <-c(toupper(fig2$template),toupper(fig4$template),toupper(fig10$template))
gRNAS_PAM <-data.frame(Run,gRNA,templates)

colnames(gRNAS_PAM) <- c("seq_name", "gRNA","template")
#################
### HDR       ###
#################
final_df <- data.frame()
for(chr in 420:length(samples)){
  #################
  # Load the data #
  #################
  sample_id <- samples[chr]
  sequence_name <- sample_id
  #gRNA_sequence = toupper(gRNAS_PAM$gRNA[which(gRNAS_PAM$seq_name==sequence_name)]) #for toy data
  gRNA_sequence = toupper(gRNAS_PAM$gRNA[which(gRNAS_PAM$seq_name==sequence_name)]) #for real data
  # z <- paste("/Users/ASUS/Desktop/BIOINFO/PRBB/HDR/tests_clustering/",sample_id,"/",sep = "")
  # setwd(z)
  # ref_fasta = "reference.fa" #toy
  # temp = "template.fasta" #toy
  #FASTA
  ref_fasta = paste("Data/",sample_id,"_reference-correctOrient.fasta", sep = "") #real
  template = paste("Data/",sample_id,"_template.fasta", sep = "") #real

  #SEQUENCES
  temp_seq <- sread(readFasta(template))[[1]]
  ref_seq <- sread(readFasta(ref_fasta))[[1]]
  ref <- toupper(read.fasta(ref_fasta,as.string = T))
  
  #CUT_SITE
  cut_site = get_cutSite(gRNA_seq = gRNA_sequence,reference = readFasta(ref_fasta))
  
  #A PRIORI GROUPS
  indels_df <- read.csv(paste("Results/",sample_id,"_indels.csv", sep = ""))
  template_df <- read.csv(paste("../Templata_reads/",sample_id,"_template-reads.csv",sep=""))
  ### Generate new reference (reference with change done by the template)
  newRef_Results <- newReference(template_sequence = temp_seq, reference_sequence = ref_seq)
  NewRef <- newRef_Results[1]
  
  # Bam of the Reference vs reads
  bamPath <- paste("Data/",sample_id,".sorted.bam", sep = "")
  bamFile <- BamFile(bamPath)
  bam <- scanBam(bamFile)
  RefvsRead_bam <- bam[[1]]
  
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
  
  
  # ##########################################################################################
  # # 2 Build the Imperfect_matches dataframe and the count of the perfect-imperfect mathces #
  # ##########################################################################################
  # match_ind <- which(RefvsRead$cigar == RefvsNew$cigar)
  # perfect_score <- 0
  # y = 0
  # 
  # Imperfect <- 0
  # imperfect_df <- data.frame(c(1:nchar(NewRef)),c(rep(0,nchar(NewRef))),c(rep(0,nchar(NewRef))),c(rep(0,nchar(NewRef))),c(rep(0,nchar(NewRef))))
  # imperfect_df$sequence <- NA
  # colnames(imperfect_df) <- c("Position","A","C","T","G","sequence")
  # for (i in 1:length(match_ind)){
  #   
  #   RR <- RefvsRead_bam$id[match_ind[i]]
  #   if (RR==NewRef){
  #     perfect_score = perfect_score + 1
  #   }else if(!is.na(RR)){
  #     y <-  y+1
  #     dif_pos <- c()
  #     dif_char <- c()
  #     char <- 1
  #     N <- toupper(strsplit(NewRef,"")[[1]])
  #     I <- toupper(strsplit(RR,"")[[1]])
  #     for(j in 1:length(N)){
  #       imperfect_df$sequence[j] <- N[j]
  #       if(N[j]!=I[j] && (I[j]!="-" && N[j]!="-")){
  #         if(I[j]=="A"){
  #           imperfect_df$A[j] = imperfect_df$A[j] + 1
  #         }else if(I[j]=="C"){
  #           imperfect_df$C[j] = imperfect_df$C[j] + 1          
  #         }else if(I[j]=="G"){
  #           imperfect_df$G[j] = imperfect_df$G[j] + 1
  #         }else{
  #           imperfect_df$T[j] = imperfect_df$T[j] + 1
  #         }
  #       }
  #     }
  #     Imperfect <- Imperfect + 1 
  #   }
  # }
  # 
  # perfect_score
  # Imperfect
  # imperfect_df
  # 
  
  #######################
  ## 2.2 Modification detection #
  #######################

  #######################################################################
  # 3 Obtain the count for all the modifications that occur in the data #
  #######################################################################
  groups_df <- RefvsRead_bam
  groups_df$group <- "unknown"
  cigar_lengths=explodeCigarOpLengths(groups_df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
  cigar_types=explodeCigarOps(groups_df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
  
  ins_df <- indels_df[which(indels_df$Modification=="ins"),]
  del_df <- indels_df[which(indels_df$Modification=="del"),]
  template_based_qnames<- intersect(groups_df$qname,template_df[,2])
  dels_qnames <- intersect(groups_df$qname,del_df$Ids)
  ins_qnames <- intersect(groups_df$qname,ins_df$Ids)
  
  for (i in 1:length(groups_df$pos)){
    #GROUP
    if(groups_df$qname[i] %in% template_based_qnames){
      groups_df$group[i] <- "template"
    }else if(groups_df$qname[i] %in% ins_qnames){
      groups_df$group[i] <- "ins"
    }else if(groups_df$qname[i] %in% dels_qnames){
      groups_df$group[i] <- "del"
    }
  }

  seq_freq <- as.data.frame(Freq(groups_df$id))
  seq_freq$freq <- (seq_freq$freq/sum(seq_freq$freq))*100
  seq_freq <- seq_freq[which(seq_freq$freq > 0.01),]
  df_with_freq <- data.frame(seq_freq$level,seq_freq$freq)
  colnames(df_with_freq) <- c("ID","Freq")
  df_with_freq$group <- "-"
  df_with_freq$cigar <- "-"
  df_with_freq$pos <- 0
  for(i in 1:length(df_with_freq$ID)){
    df_with_freq$cigar[i] <- groups_df$cigar[which(groups_df$id == df_with_freq$ID[i])[1]]
    df_with_freq$group[i] <- groups_df$group[which(groups_df$id == df_with_freq$ID[i])[1]]
    df_with_freq$pos[i] <- groups_df$pos[which(groups_df$id == df_with_freq$ID[i])[1]]
  }
  df <- df_with_freq
  which(df$group == "template")
  ################################
  #OLD WAY OF OBTAINING THE DATAFRAME TOO MUCH TIME CONSUMING  
  # df <-as.data.frame(table(groups_df$pos,groups_df$cigar,groups_df$group, groups_df$id)) #, groups_df$id deberia incluir el id? o buscar el id en base al cigar pese a que ignore posibles substituciones?, este cambio lo haria mucho más rapido y menos pesado, pero daria problemas en las substituciones
  # df$Freq <- (df$Freq/sum(df$Freq))*100
  # df <- df[which(df$Freq>0.05),]
  ################################
  
  
  #tmp_ind <- which(df$group=="template")
  #maxim <- max(df$Freq[tmp_ind]) #which(df$Freq == max(df$Freq[which(df$Var3=="template")]))
  #maxim/sum(df$Freq[tmp_ind])
  #to_eliminate <- tmp_ind[which(df$Freq[tmp_ind] != maxim)]
  #if(length(to_eliminate)!=0)  df <- df[-to_eliminate,]
  

  df$modif <- NA
  df$cut_site_dist <- 0
  df$size <- 0
  df$baseModif <- "-"
  cigar_lengths=explodeCigarOpLengths(df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
  cigar_types=explodeCigarOps(df$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
  

  df[,'group'] <- as.character(df[,'group'])
  df[,'Freq'] <- as.numeric(as.character(df[,'Freq']))
  
  pos_subs <- c()
  
  for (i in 1:dim(df)[1]){

    if (length(cigar_types[[i]]) == 1 && cigar_types[[i]] == "M"){ #CHECK WT
      if(df$group[i]=="unknown"){ 
        df$group[i] <- "wt"

      }
      R <- strsplit(ref, "")[[1]][df$pos[i]:(cigar_lengths[[i]])] #REFERENCE
      I <- strsplit(as.character(df$ID[i]), "")[[1]][1:(cigar_lengths[[i]])] #ID
      if(length(R)==length(I)){
       pos_subs <- which(R!=I)
      }
      if(length(pos_subs)==0){
       df$group[i] <- "wt"
       df$modif[i] <- "wt" 
       df$cut_site_dist[i] = nchar(ref) #setting the cut site dist to maximum
       #NO NEED TO MODIFY SIZE AS FOR WT IT STAYS 0
      }else if(is.na(R[pos_subs])==F && is.na(I[pos_subs])==F){ #CHECK SUBSTITUTIONS
       sub <- pos_subs[which(abs(pos_subs - cut_site) == min(abs(pos_subs - cut_site)))] #GET THE CLOSER SUBS TO THE CUTSITE 
       change <- paste(R[sub[1]],">",I[sub[1]],sep = "")
       df$modif[i] <- "subs"
       df$cut_site_dist[i] = sub[1]-cut_site #setting the cut site dist 
       df$size[i] <- 1 #FOR THE MOMENT SUBS WILL BE 1 THIS WILL CHANGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       df$baseModif[i] <- change
      }else{
        df$modif[i] <- "wt"
        df$cut_site_dist[i] = nchar(ref)
      }
    
    }else if(length(cigar_types[[i]]) <= 3 && cigar_types[[i]][2] == "D"){
        df$modif[i] <- "del"
        if(df$group[i]=="unknown"){
          df$group[i] <- "del"
        }
    }
     else if(length(cigar_types[[i]]) <= 3 && cigar_types[[i]][2] == "I"){
        df$modif[i] <- "ins"
        if(df$group[i]=="unknown"){
          df$group[i] <- "ins"
        }
     }
     else if(length(cigar_types[i][[1]]) > 3){
       df$modif[i] <- "Delin"
       if(df$group[i]=="unknown"){
         df$group[i] <- "Delin"
       }
    }
  }
  
  df_IDS <- df
  
  #####
  # 4 Filter the dataframe #
  #####
  
  
  df_filtered <- as.data.frame(table(df$cigar,df$group,df$modif,df$baseModif,df$cut_site_dist, df$size))
  colnames(df_filtered) <- c("cigar", "group", "modif", "baseModif", "cut_site_dist", "size", "Freq")
  ref_numbers = which(df_filtered$Freq>5)
  df_filtered <- df_filtered[which(df_filtered$Freq>5),]
  
  df<-df_filtered
  
  #REMOVAL OF FACTORS IN NUMERIC COLUMNS
  df[,'size'] <- as.numeric(as.character(df[,'size']))
  df[,'cut_site_dist'] <- as.numeric(as.character(df[,'cut_site_dist']))
  
  df$wt_similarity <- NA
  df$GC_content <- NA
  
  #####
  # 5 Build the data frame #
  #####
  wt_seq <- groups_df$id[which(df$group=="wt")[1]]
  if(is.na(wt_seq)){
    wt_seq <-ref
  }  
  for(i in 1:length(df$cigar)){
    data_related <- df_IDS %>% filter(grepl(df_filtered$cigar[i],cigar) )  %>% filter(grepl(df_filtered$modif[i],modif)) %>% filter(grepl(df_filtered$group[i],group)) %>% filter(grepl(df_filtered$baseModif[i],baseModif)) %>% filter(grepl(df_filtered$cut_site_dist[i],cut_site_dist)) %>% filter(grepl(df_filtered$size[i],size)) 
    tmp_seq <- data_related$ID[1] #THIS IS WRONG I THINK
    alignment <- pairwiseAlignment(wt_seq, tmp_seq, type="local")
    df$wt_similarity[i] <- score(alignment)
    
    #Obtain the GC content per sequence 
    num_g <- str_count(toupper(tmp_seq), "G")
    num_c <- str_count(toupper(tmp_seq), "C")
    gc_content <- (num_g + num_c) / str_length(tmp_seq) * 100 
    df$GC_content[i] <- gc_content
    
    #Add the size and cut_site_distance to the data frame 
    cig <- df$cigar[i]
    lens=explodeCigarOpLengths(cig)[[1]] ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
    types=explodeCigarOps(cig)[[1]] ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
    
    #compute size
    if(df$modif[i]!="wt" && df$modif[i]!="subs"){
      #compute size
      ind <- which(types!="M")
      df$size[i] <- sum(lens[ind])
      #compute distance
      dist = 0
      vdist = c()
      for(j in 1:length(types)){
        if(types[j]=="M"){
           dist = dist + lens[j] 
          }
        vdist = c(vdist,dist)
       }
       df$cut_site_dist[i] <- mean(vdist) - cut_site  
    }
  }
  
  #GC content and wt similarity naturalization
  df$wt_similarity <- df$wt_similarity/max(df$wt_similarity)
  df$GC_content <- df$GC_content/max(df$GC_content)
  # #####
  # # 6 #
  # #####
  # 
  # for (i in 1:length(df_filtered$cigar)){
  #   cig <- df_filtered$cigar[i]
  #   lens=explodeCigarOpLengths(cig) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
  #   types=explodeCigarOps(cig) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
  #   #compute size
  #     ind <- which(types[[1]]!="M")
  #     df_filtered$size[i] <- sum(lens[[1]][ind])
  #     
  #     #compute distance
  #     dist = 0
  #     vdist = c()
  #     for(j in 1:length(types[[1]])){
  #       if(types[[1]][j]=="M"){
  #         dist = dist + lens[[1]][j] 
  #       }else{
  #         vdist = c(vdist,dist)
  #       }
  #     }
  #     df_filtered$cut_site_dist[i] <- mean(vdist) - cut_site
  # }
  
  #####
  # 7 treat the data #
  #####
  
  #As factors
  df$cigar <- as.factor(df$cigar)
  df$group <- as.factor(df$group)
  df$baseModif <- as.factor(df$baseModif)
  df$modif <- as.factor(df$modif)
  final_df <- rbind(final_df,df)
  
  
}
colnames(final_df)
###
#USE SUBSTITUTIONS FROM SAMPLE 121 
###

final_df <- rbind(final_df,read.csv("../Real_Dataframe.csv"))
final_df <- final_df[which(final_df$Freq>100),]
#write.csv(final_df,"../Toy_Dataframe.csv", row.names = FALSE)

final_df <- droplevels(final_df) #SOLUTION TO Error in randomForest.default(m, y, ...) : Can't have empty classes in y.
#final_df$group  <- factor(final_df$group)

# final_df$simp_cig <- NA
# for (i in 1:length(cigar_types_final)){
#   if (length(cigar_types_final[[i]])>1 && length(cigar_types_final[[i]])<=3){  
#     ci <- cigar_types_final[[i]][cigar_types_final[[i]] != "M" ]
#     ci <- paste(ci,sep = "", collapse = "")
#     final_df$simp_cig[i] <- ci
#     
#   }else if(length(cigar_types_final[[i]])==1){
#     final_df$simp_cig[i] <- "M"
#   }else{
#     final_df$simp_cig <- "InDel"
#   }
# }
# final_df <- final_df[rowSums(is.na(final_df)) != ncol(final_df),]
# df_filtered <- df_filtered[-which(is.na(df_filtered)==T),]

#data partition
table(final_df$modif,final_df$group)



###########################
# DELETIONS DECISION TREE #
###########################

deletions_df <- final_df[which(final_df$modif=="del"),]#separate the DEATH_EVENT equal to 1
d_1 <- which(deletions_df$group=="template")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_del <- deletions_df[train_ind, ] #separate between train and test
test_del <- deletions_df[setdiff(d_1,train_ind),]

d_1 <- which(deletions_df$group=="del")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_del2 <- deletions_df[train_ind, ] #separate between train and test
test_del2 <- deletions_df[setdiff(d_1,train_ind),]

train <- rbind(train_del,train_del2)
test <- rbind(test_del,test_del2)

train <- droplevels(train) 
test <- droplevels(test)


train$modif <- as.factor(train$modif)
train$baseModif <- as.factor(train$baseModif)
train$cut_site_dist <- as.numeric(as.character(train$cut_site_dist))
train$size <- as.numeric(as.character(train$size))
train$Freq <- as.numeric(as.character(train$Freq))
train$wt_similarity <- as.numeric(as.character(train$wt_similarity))
train$GC_content <- as.numeric(as.character(train$GC_content))
train$group <- as.factor(train$group)

test$modif <- as.factor(test$modif)
test$baseModif <- as.factor(test$baseModif)
test$cut_site_dist <- as.numeric(as.character(test$cut_site_dist))
test$size <- as.numeric(as.character(test$size))
test$Freq <- as.numeric(as.character(test$Freq))
test$wt_similarity <- as.numeric(as.character(test$wt_similarity))
test$GC_content <- as.numeric(as.character(test$GC_content))
test$group <- as.factor(test$group)

# DECISION TREES
train$baseModif <- as.factor(train$baseModif)
test$baseModif <- as.factor(test$baseModif)
tree <- rpart(formula = group ~ Freq + cut_site_dist + size + wt_similarity + GC_content
              , data = train, method='class') #run rpart() function to build my tree

fancyRpartPlot(tree) #ploting the tree
tree.pred <- predict(tree, test, type = "class") #building th predict
table(tree.pred, test$group)

# RANDOM FOREST
Rforest <- randomForest(formula = group ~ Freq + cut_site_dist + size + wt_similarity + GC_content
                        , data = train, method='class') #run rpart() function to build my tree

Rforest.pred <- predict(Rforest, test, type = "class") #building th predict
table(Rforest.pred, test$group)

###########################
# Insertion DECISION TREE #
###########################

insertions_df <- final_df[which(final_df$modif=="ins"),]#separate the DEATH_EVENT equal to 1
d_1 <- which(insertions_df$group=="template")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_ins <- insertions_df[train_ind, ] #separate between train and test
test_ins <- insertions_df[setdiff(d_1,train_ind),]

d_1 <- which(insertions_df$group=="ins")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_ins2 <- insertions_df[train_ind, ] #separate between train and test
test_ins2 <- insertions_df[setdiff(d_1,train_ind),]

train <- rbind(train_ins,train_ins2)
test <- rbind(test_ins,test_ins2)

train <- droplevels(train) 
test <- droplevels(test)


train$modif <- as.factor(train$modif)
train$baseModif <- as.factor(train$baseModif)
train$cut_site_dist <- as.numeric(as.character(train$cut_site_dist))
train$size <- as.numeric(as.character(train$size))
train$Freq <- as.numeric(as.character(train$Freq))
train$wt_similarity <- as.numeric(as.character(train$wt_similarity))
train$GC_content <- as.numeric(as.character(train$GC_content))
train$group <- as.factor(train$group)

test$modif <- as.factor(test$modif)
test$baseModif <- as.factor(test$baseModif)
test$cut_site_dist <- as.numeric(as.character(test$cut_site_dist))
test$size <- as.numeric(as.character(test$size))
test$Freq <- as.numeric(as.character(test$Freq))
test$wt_similarity <- as.numeric(as.character(test$wt_similarity))
test$GC_content <- as.numeric(as.character(test$GC_content))
test$group <- as.factor(test$group)

# DECISION TREES
train$baseModif <- as.factor(train$baseModif)
test$baseModif <- as.factor(test$baseModif)
tree <- rpart(formula = group ~ Freq + cut_site_dist + size + wt_similarity + GC_content
              , data = train, method='class') #run rpart() function to build my tree

fancyRpartPlot(tree) #ploting the tree
tree.pred <- predict(tree, test, type = "class") #building th predict
table(tree.pred, test$group)

# RANDOM FOREST
Rforest <- randomForest(formula = group ~ Freq + cut_site_dist + size + wt_similarity + GC_content
                        , data = train, method='class') #run rpart() function to build my tree

Rforest.pred <- predict(Rforest, test, type = "class") #building th predict
table(Rforest.pred, test$group)

##############################
# SUBSTITUTION DECISION TREE #
##############################
final_df <- final_df[which(final_df$Freq>100),]
substitutions_df <- final_df[which(final_df$modif=="subs"),]#separate the DEATH_EVENT equal to 1
d_1 <- which(substitutions_df$group=="template")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_subs <- substitutions_df[train_ind, ] #separate between train and test
test_subs <- substitutions_df[setdiff(d_1,train_ind),]

d_1 <- which(substitutions_df$group=="wt")#separate the DEATH_EVENT equal to 1
smp_size <- floor(0.7 * length(d_1)) #get the sample size desired in this case 0.7
train_ind <- sample(d_1, size = smp_size) #sample the amount specified before
train_subs2 <- substitutions_df[train_ind, ] #separate between train and test
test_subs2 <- substitutions_df[setdiff(d_1,train_ind),]

train <- rbind(train_subs,train_subs2)
test <- rbind(test_subs,test_subs2)

train <- droplevels(train) 
test <- droplevels(test)


train$modif <- as.factor(train$modif)
train$baseModif <- as.factor(train$baseModif)
train$cut_site_dist <- as.numeric(as.character(train$cut_site_dist))
train$size <- as.numeric(as.character(train$size))
train$Freq <- as.numeric(as.character(train$Freq))
train$wt_similarity <- as.numeric(as.character(train$wt_similarity))
train$GC_content <- as.numeric(as.character(train$GC_content))
train$group <- as.factor(train$group)

test$modif <- as.factor(test$modif)
test$baseModif <- as.factor(test$baseModif)
test$cut_site_dist <- as.numeric(as.character(test$cut_site_dist))
test$size <- as.numeric(as.character(test$size))
test$Freq <- as.numeric(as.character(test$Freq))
test$wt_similarity <- as.numeric(as.character(test$wt_similarity))
test$GC_content <- as.numeric(as.character(test$GC_content))
test$group <- as.factor(test$group)

# DECISION TREES
train$baseModif <- as.factor(train$baseModif)
test$baseModif <- as.factor(test$baseModif)
tree <- rpart(formula = group ~ Freq + cut_site_dist + size + baseModif + wt_similarity + GC_content
              , data = train, method='class') #run rpart() function to build my tree

fancyRpartPlot(tree) #ploting the tree
tree.pred <- predict(tree, test, type = "class") #building th predict
table(tree.pred, test$group)

# RANDOM FOREST
Rforest <- randomForest(formula = group ~ Freq + cut_site_dist + size + wt_similarity + GC_content
                        , data = train, method='class') #run rpart() function to build my tree

Rforest.pred <- predict(Rforest, test, type = "class") #building th predict
table(Rforest.pred, test$group)



#write.csv(final_df,"FINAL_DF4.csv", row.names = FALSE)

f1 <- read.csv("FINAL_DF1.csv")
f2 <- read.csv("FINAL_DF2.csv")
f3 <- read.csv("FINAL_DF3.csv")
f4 <- read.csv("FINAL_DF4.csv")
final_df <- rbind(f4,f1,f2,f3)
