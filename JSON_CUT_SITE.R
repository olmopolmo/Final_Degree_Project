#CUT_SITE
cut_site = get_cutSite(gRNA_seq = gRNA_sequence,reference = readFasta(ref_fasta))
cut_site = NULL
library(jsonlite)
exportJson <- toJSON(cut_site)
write(exportJson, "cutSite.json")

#ADD THIS TO THE PHP:º