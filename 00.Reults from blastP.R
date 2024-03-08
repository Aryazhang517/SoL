######################################
## 01 process the raw data from BlastP
######################################

library(data.table)
require(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggsci)
library(reshape2)
library(ggplot2) 
library(ggsci)

Sol_0 <- read.delim('Sulfonolipid_BE_output', sep = '\t', check.names = FALSE)
Sol_2 <- as.data.table(Sol_0)
##select the row with the maximum value in each group
Sol_3 <- Sol_2[Sol_2[, .I[bitscore == max(bitscore)], by=qseqid]$V1]  # Top N highest values by group
Sol_40 <- Sol_3[!duplicated(Sol_3$qseqid),]
##add genome
Sol_40$genome <-sapply(strsplit(as.character(Sol_40$qseqid), '_'), function(x) x[2])
Sol_40$genome <- str_replace(Sol_40$genome, "GENOME", "GUT_GENOME")
##add genome annotation
genome28 <- read.delim("Genome28info",sep="\t")
genome28_1 <- genome28[,c(1,2,6)]

setDT(genome28_1)[, paste0("TM", 1:7) := tstrsplit(genome28_1$Genome_Taxonomy_Database_lineage, ";")]
colnames(genome28_1)[4:10] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
genome_2 <- genome28_1[,c(1,4,5)]
bacteria_genome <- genome28_1[!genome28_1$Kingdom == "d__Archaea",] ## remove genomes from archaea

Sol_4 <- merge(Sol_40, bacteria_genome,by="genome")

CFAT <- Sol_4[Sol_4$sseqid %in% c('EFK34533.1','SEA34915.1','EAZ83176.1','ABQ05446.1','MBM7419832.1','AFL77572.1'), ]
CFAT$type <- "CFAT"
CYS <- Sol_4[Sol_4$sseqid %in% c('EFK34534.1','SEA02653.1','EAZ80562.1','MBM7419833.1'), ]
CYS$type <- "CYS"
SDR <- Sol_4[Sol_4$sseqid %in% c('QAR30703.1'), ]
SDR$type <- "SDR"

colnames(CFAT)==colnames(CYS)
colnames(CFAT)==colnames(SDR)
Sol_m <- rbind(CYS,CFAT,SDR)

Sol_genome_stat <- table(Sol_m$genome,Sol_m$type)
write.table(Sol_genome_stat,file="Sol_genome_stat.tsv",sep ="\t", row.names=T,quote = F)
Sol_genome_stat <- read.delim("Sol_genome_stat.tsv",sep="\t")

genome_f <- Sol_m[!duplicated(Sol_m$genome),]
genome_no <- as.data.frame(table(genome_f$Genome_Taxonomy_Database_lineage))
colnames(genome_no) <- c("taxa","genome_no")

df <- Sol_genome_stat
df$no <- apply(df,1,function(x) sum(x>0))
table(df$no)

data <- df %>% 
  group_by(no) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

#####Supplementary Figure 1a #####
pdf(file="Proportion of enzyme types per genomes.pdf",width = 4,height = 4)
ggplot(data, aes(x = "", y = perc, fill = labels)) +
  geom_col(color = "black") +
  geom_label(aes(label = labels), color = c(1, "white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  guides(fill = guide_legend(title = "types")) +
  scale_fill_lancet() +
  coord_polar(theta = "y") + 
  theme_void()
dev.off()

#####Supplementary Figure 1b #####
##input normalized data (Total Enzyme number/ Number of Genomes in each phylum)
dfn <- read.table('Phylum_normalized_data', header = T)
data <- read.table('enzyme_phylum_stat',sep="\t",header = T)

##Bubblematrix
data_melt<-melt (data)
names(data_melt) = c('Biosynthetic_enzyme','Phylum','Proportion')
data_melt <- data.frame(data_melt)
data_melt$Biosynthetic_enzyme <- factor(data_melt$Biosynthetic_enzyme,levels=c("CYS","CFAT","SDR"), ordered = TRUE)

pdf(file="BE distribution in phylum.pdf",width = 6,height = 5)
ggplot(data_melt, aes(x =Biosynthetic_enzyme , y = Phylum , size = Proportion, color= Phylum)) + geom_point()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA)) +
  scale_size_continuous(range = c(2, 15)) + scale_color_lancet()
#scale_color_d3()
#scale_color_aaas()
#scale_color_lancet()
dev.off()

