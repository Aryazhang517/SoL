#########################
## hmmserach output ###
########################
library(data.table)
library(stringr)
library(tibble)
library(stringr)
library(tidyverse)
library(dplyr)
library(plyr)

CYS_1 <- read.delim('CYS_hmmsearch_PF00291',sep="\t", check.names = FALSE)
CFAT_1 <- read.delim('CFAT_hmmsearch_PF00155',sep="\t", check.names = FALSE)
SDR_1 <- read.delim('SDR_hmmsearch_PF00106',sep="\t", check.names = FALSE)

####SDR#####
if(T){
  SDR_1$id <- sapply(strsplit(as.character(SDR_1$query_name), " "), function(x) x[1])
  SDR_2 <- SDR_1[ ,c(4,3)]
  ###Split delimited strings in a column and insert as new rows 
  SDR_3 <- SDR_2 %>% 
    mutate(pfams = strsplit(as.character(pfams), " ")) %>% 
    unnest(pfams)
  
  SDR_3$pfams_1 <- gsub("[()]", " ", SDR_3$pfams)
  setDT(SDR_3)[, paste0("TM", 1:2) := tstrsplit(SDR_3$pfams_1, " ")]
  SDR_4 <- SDR_3[ ,c(1,4,5)]
  colnames(SDR_4) <- c("ID","pfam","scores")
  SDR_5 <- SDR_4[SDR_4$pfam=="PF00106", ]
  SDR_5$scores <- as.numeric(as.character(SDR_5$scores))
  SDR_6 <- SDR_5[SDR_5$scores > 50, ]
  SDR_7 <- SDR_6[!duplicated(SDR_6$ID), ]
}
####CYS#####
if(T){
CYS_1$id <- sapply(strsplit(as.character(CYS_1$query_name), " "), function(x) x[1])
CYS_2 <- CYS_1[ ,c(4,3)]
###Split delimited strings in a column and insert as new rows 
CYS_3 <- CYS_2 %>% 
  mutate(pfams = strsplit(as.character(pfams), " ")) %>% 
  unnest(pfams)

CYS_3$pfams_1 <- gsub("[()]", " ", CYS_3$pfams)
setDT(CYS_3)[, paste0("TM", 1:2) := tstrsplit(CYS_3$pfams_1, " ")]
CYS_4 <- CYS_3[ ,c(1,4,5)]
colnames(CYS_4) <- c("ID","pfam","scores")
CYS_5 <- CYS_4[CYS_4$pfam=="PF00291", ]
CYS_5$scores <- as.numeric(as.character(CYS_5$scores))
CYS_6 <- CYS_5[CYS_5$scores > 50, ]
CYS_7 <- CYS_6[!duplicated(CYS_6$ID), ]
}
####CFAT#####
if(T){
CFAT_1$id <- sapply(strsplit(as.character(CFAT_1$query_name), " "), function(x) x[1])
CFAT_2 <- CFAT_1[ ,c(4,3)]
###Split delimited strings in a column and insert as new rows 
CFAT_3 <- CFAT_2 %>% 
  mutate(pfams = strsplit(as.character(pfams), " ")) %>% 
  unnest(pfams)

CFAT_3$pfams_1 <- gsub("[()]", " ", CFAT_3$pfams)
setDT(CFAT_3)[, paste0("TM", 1:2) := tstrsplit(CFAT_3$pfams_1, " ")]
CFAT_4 <- CFAT_3[ ,c(1,4,5)]
colnames(CFAT_4) <- c("ID","pfam","scores")
CFAT_5 <- CFAT_4[CFAT_4$pfam=="PF00155", ]
CFAT_5$scores <- as.numeric(as.character(CFAT_5$scores))
CFAT_6 <- CFAT_5[CFAT_5$scores > 50, ]
CFAT_7 <- CFAT_6[!duplicated(CFAT_6$ID), ]
}
###
##add genome info
Sol_7 <- rbind(CYS_7,CFAT_7,SDR_7)
Sol_7$genome <-sapply(strsplit(as.character(Sol_7$ID), '_'), function(x) x[2])
Sol_7$genome <- str_replace(Sol_7$genome, "GENOME", "GUT_GENOME")
genome28 <- read.delim("Genome28info",sep="\t")
genome28_1 <- genome28[,c(1,2,6)]
setDT(genome28_1)[, paste0("TM", 1:7) := tstrsplit(genome28_1$Genome_Taxonomy_Database_lineage, ";")]
colnames(genome28_1)[4:10] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
genome_2 <- genome28_1[,c(1,4,5)]
bacteria_genome <- genome28_1[!genome28_1$Kingdom == "d__Archaea",]

Sol_b <- Sol_7[Sol_7$genome %in% bacteria_genome$genome, ]
Sol_genome_stat <- table(Sol_b$genome,Sol_b$pfam)
write.table(Sol_genome_stat,file="Sol_genome_stat.tsv",sep ="\t", row.names=T,quote = F)
Sol_genome_stat <- read.delim("Sol_genome_stat.tsv",sep="\t")
#####genome contain 3Es#####
Sol_3E_genome_pfam <- Sol_genome_stat[apply(Sol_genome_stat,1, function(x) all(x!=0)),]
Sol_b_f3e <- Sol_b[Sol_b$genome %in% row.names(Sol_3E_genome_pfam), ]
table(Sol_b_f3e$pfam)
Sol_b_f3e_diamond <- Sol_3eh_diamond[Sol_3eh_diamond$qseqid %in% Sol_b_f3e$ID, ]

###sequence####
sequence <- read.delim2("Sol_pfam50_sequence.tsv")
colnames(sequence)[1] <- "qseqid"
pair <- read.delim2("Sol_pfam_cluster_id_pair.tsv",header = F)
colnames(pair) <- c("family_id","qseqid")
merged_1 <- merge(Sol_b_f3e_diamond,sequence,by="qseqid")
Sol_seq <- merge(pair,merged_1,by="qseqid")

#########################
## taxonomy###
######################

taxa_3eh <- table(Sol_b_f3e_diamond$Genome_Taxonomy_Database_lineage,Sol_b_f3e_diamond$type)
write.table(taxa_3eh,file="taxa_3eh",sep="\t",quote = F)
Sol_3ehpfam_genome_1 <- read.delim("taxa_3eh",sep="\t")
Sol_3ehpfam_genome_1 <- tibble::rownames_to_column(Sol_3ehpfam_genome_1,var="taxa")
###represent_genome
represent_genome <- read.delim("representative_genome",sep="\t")
Sol_3ehpfam_genome_2 <- merge(Sol_3ehpfam_genome_1,represent_genome,by = "taxa")
##genome_number
genome_f3eh <- Sol_b_f3e_diamond[!duplicated(Sol_b_f3e_diamond$genome),]
genome_3ehno <- as.data.frame(table(genome_f3eh$Genome_Taxonomy_Database_lineage))
colnames(genome_3ehno) <- c("taxa","genome_no")

Sol_3ehpfam_genome_3 <- merge(Sol_3ehpfam_genome_2,genome_3ehno,by="taxa")
Sol_3ehpfam_genome_4 <- Sol_3ehpfam_genome_3[!duplicated(Sol_3ehpfam_genome_3$taxa), ]
setDT(Sol_3ehpfam_genome_4)[, paste0("TM", 1:7) := tstrsplit(Sol_3ehpfam_genome_4$taxa, ";")]
colnames(Sol_3ehpfam_genome_4)[7:13] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
Sol_3ehpfam_genome_4$Species[Sol_3ehpfam_genome_4$Species == 's__'] <- 'species_unclassified'
Sol_3ehpfam_genome_4$Species <- str_replace(Sol_3ehpfam_genome_4$Species, "s__", "")

write.table(Sol_3ehpfam_genome_4,file="Sol_3ehpfam_genome_annotation.tsv",sep="\t",quote = F,row.names = F)
save(Sol_3eh,Sol_b_f3e_diamond,Sol_3ehpfam_genome_4,file="Sol_enzyme_3E_genome_higher_similarity_pfam50.Rdata")





