library(dplyr)
library(tidyr)
library(tidyverse)
library(ggsci)
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(vegan)
library(ggsci)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)

######### TPM ###########
Sulfonolipid_BE <-read.table("IBD_MT_SolBE_raw_counts_f.csv",sep="\t",row.names = 1,header=T,encoding="UTF-8")
##calculate tpm
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

ftr.tpm <- Sulfonolipid_BE %>%
  gather(sample, cnt, 7:ncol(Sulfonolipid_BE)) %>%
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)

ftr.tpm[is.na(ftr.tpm)]<- 0
ftr.tpm_out <- select(ftr.tpm,c(1,7:ncol(ftr.tpm)))
write.table(ftr.tpm_out,file = "IBD_p1p2_MT_TPM_all.tsv", sep="\t",quote = F, row.names = F)

### Relative abundance of clusters ###############
countdata_1 <-read.csv("IBD_TPM_rev_cluster90_pfam50_MT.csv",sep="\t",row.names = 1)
countdata_2 <-countdata_1[rowSums(countdata_1) != 0,]
countdata_2 <-countdata_2[,colSums(countdata_2) != 0]
countdata_30 <- countdata_2[rowSums(countdata_2>0)>ncol(countdata_2)*0.05,]
###metadata
metadata <- read.table("IBD_metadata",sep="\t",header=T,encoding="UTF-8")
colnames(metadata) <- c("id","condition")

IBD <- metadata %>% filter(metadata$condition %in% c("IBD"))
Ctr <- metadata %>% filter(metadata$condition %in% c("nonIBD"))

countdata_IBD <- countdata_30[,as.vector(IBD$id)]
countdata_IBD <-countdata_IBD[,colSums(countdata_IBD) != 0]
countdata_Ctr <- countdata_30[,as.vector(Ctr$id)]
countdata_Ctr <-countdata_Ctr[,colSums(countdata_Ctr) != 0]

countdata_IBD$ID <- rownames(countdata_IBD)
countdata_Ctr$ID <- rownames(countdata_Ctr)

####### Prevalence #######
countdata_IBD$percent <- (rowSums(countdata_IBD>0)/ncol(countdata_IBD))*100
countdata_Ctr$percent <- (rowSums(countdata_Ctr>0)/ncol(countdata_Ctr))*100

Sol_BE_perc_IBD <- countdata_IBD %>% select(c(percent)) %>% setNames(c("percent_IBD"))
Sol_BE_perc_Ctr <- countdata_Ctr %>% select(c(percent)) %>% setNames(c( "percent_Ctr"))

Sol_BE_cluster_IBD <- merge(Sol_BE_perc_Ctr,Sol_BE_perc_IBD,by="row.names")
colnames(Sol_BE_cluster_IBD)[1] <-c("SOl_BE")

#####Fisher.exact test######
##count proportion
df_ctr <- countdata_Ctr[ ,c(1:2)]
df_ctr$ctr <- rowSums(countdata_Ctr[ ,1:(ncol(countdata_Ctr))]>0)
df_ctr$nonctr <- (ncol(countdata_Ctr))- df_ctr$ctr
df_ctr <- df_ctr[ ,c(3,4)]

df_case <- countdata_IBD[ ,c(1:2)]
df_case$case <- rowSums(countdata_IBD[ ,1:(ncol(countdata_IBD))]>0)
df_case$noncase <-  (ncol(countdata_IBD))- df_case$case
df_case <- df_case[ ,c(3,4)]  

df <- merge(df_ctr,df_case,by="row.names")
df <- df[ ,c(1,2,4,3,5)]
colnames(df)[1] <- "cluster"
##Fisher.exact test###
analysis  <- df %>% 
  setNames(c("cluster","control","case","nonctrol","noncase")) %>% 
  nest(-cluster) %>% 
  mutate(matrix = map(data, ~matrix(unlist(.x), nrow = 2))) %>% 
  mutate(fisher = map(matrix, ~fisher.test(.x, alternative = "two.sided"))) %>% 
  mutate(stats = map(fisher, ~broom::glance(.x)))

analysis_1 <- analysis %>% 
  unnest(stats) %>%
  select(cluster, p.value, odds = estimate)

Fisher <- merge(df,analysis_1,by = "cluster")
Fisher_1 <- Fisher[Fisher$p.value < 0.05, ]
Fisher_2 <- Fisher[ ,c(1,6,7)]
colnames(Fisher_2)[2] <- "p.value_fisher"
###### Abundance #####
group <- tibble::column_to_rownames(metadata,var="id")
countdata_3T <- t(countdata_30)
for_stat <- merge (countdata_3T,group,by="row.names",all.x=TRUE )
colnames(for_stat)[c(1)] <- "ID"
colnames(for_stat)[ncol(for_stat)] <- "Group"
write.table(for_stat, file = 'Sol_TPM_IBD_MT_cluster90_pfam50stat.tsv',sep ="\t", row.names=F,quote = F)

##two_sided###
stat_IBD <- read.table('Sol_in_IBD_only_MT_sta_tpm_p12_pfam50.tsv',header=T,sep="\t")
colnames(stat_IBD)[1] <-c("SOl_BE")
stat_IBD$adjust_Pvalue <- p.adjust(stat_IBD$pvalue,method="BH",n=length(stat_IBD$pvalue))

Sol_BE_cluster_IBD_1 <- merge(Sol_BE_cluster_IBD,stat_IBD, by="SOl_BE")
##higher###
stat_IBD_h <- read.table('Sol_in_IBD_only_MT_sta_tpm_p12_pfam50_greater.tsv',header=T,sep="\t")
colnames(stat_IBD_h)[1] <-c("SOl_BE")
stat_IBD_h$adjust_Pvalue <- p.adjust(stat_IBD_h$pvalue,method="BH",n=length(stat_IBD_h$pvalue))
colnames(stat_IBD_h) <- c("SOl_BE","statistic_h","pvalue_h","adjust_Pvalue_h")
Sol_BE_cluster_IBD_2 <- merge(Sol_BE_cluster_IBD_1,stat_IBD_h, by="SOl_BE")
###lower##
stat_IBD_l <- read.table('Sol_in_IBD_only_MT_sta_tpm_p12_pfam50_less.tsv',header=T,sep="\t")
colnames(stat_IBD_l)[1] <-c("SOl_BE")
stat_IBD_l$adjust_Pvalue <- p.adjust(stat_IBD_l$pvalue,method="BH",n=length(stat_IBD_l$pvalue))
colnames(stat_IBD_l) <- c("SOl_BE","statistic_l","pvalue_l","adjust_Pvalue_l")
Sol_BE_cluster_IBD_3 <- merge(Sol_BE_cluster_IBD_2,stat_IBD_l, by="SOl_BE")

###add group###
Sol_BE_cluster_IBD_3$group <- "other"
Sol_BE_cluster_IBD_3[which(Sol_BE_cluster_IBD_3$adjust_Pvalue < 0.05),"group"] = "differential SoL biosynthetic family"

all_sta <- merge(Fisher_2,Sol_BE_cluster_IBD_3,by.x="cluster",by.y="SOl_BE")
all_stat_1 <- all_sta[ ,c(1:6,8,15,11,14)]
colnames(all_stat_1)[1] <- "family_id"
write.table(all_stat_1, file = 'Sol_enzyme_prevalence_abundance_all_test.tsv',sep ="\t", row.names=F,quote = F)

save(metadata,countdata_2,countdata_30,Sol_BE_cluster_IBD_3,countdata_Ctr,countdata_IBD,metadata,for_stat,file="Sol_3EH_counts_IBD_cluster90_pfam50_all.Rdata")
########################################
##### plot for differential family #####
########################################
rm(list = ls())
load("Sol_3EH_counts_IBDonly_cluster90_pfam50.Rdata")
##average abundance##
countdata_Ctr_1 <- countdata_Ctr[ ,c(1:(ncol(countdata_Ctr)-2))]
countdata_Ctr_1$rowMean <- rowMeans(countdata_Ctr_1)
countdata_Ctr_1 <- tibble::rownames_to_column(countdata_Ctr_1,var = "SOl_BE")
countdata_Ctr_2 <- countdata_Ctr_1[,c(1,ncol(countdata_Ctr_1))]
colnames(countdata_Ctr_2)[2] <- "nonIBD"

countdata_IBD_1 <- countdata_IBD[ ,c(1:(ncol(countdata_IBD)-2))]
countdata_IBD_1$rowMean <- rowMeans(countdata_IBD_1)
countdata_IBD_1 <- tibble::rownames_to_column(countdata_IBD_1,var = "SOl_BE")
countdata_IBD_2 <- countdata_IBD_1[,c(1,ncol(countdata_IBD_1))]
colnames(countdata_IBD_2)[2] <- "IBD"

Sol_abundance <- merge(countdata_Ctr_2,countdata_IBD_2,by="SOl_BE")
SOl_all <- merge(Sol_BE_cluster_IBD,Sol_abundance,by="SOl_BE")

##for differential families##
differ <- Sol_BE_cluster_IBD[Sol_BE_cluster_IBD$class=="differential Sol_BE", ]
row.names(for_stat) <- for_stat[,1]
differ_family <- for_stat[, colnames(for_stat) %in% differ$SOl_BE]
differ_family1 <- merge(differ_family,metadata,by.x = "row.names",by.y = "id")
df <- differ_family1

my_comparison <- list(c("nonIBD","IBD"))
df$condition <- factor(df$condition,
                       levels = c("IBD",
                                  "nonIBD"))

df_all <- df %>% pivot_longer(cols = colnames(df)[2:(ncol(df)-1)],
                              names_to = "family",
                              values_to = "abundance")

colnames(df_all)[1] <- "id"
df_all$condition <- factor(df_all$condition,
                           levels = c("IBD",
                                      "nonIBD"))
df_all$family <- factor(df_all$family, levels = c("SDR_family73",
                                                  "CFAT_family3",
                                                  "CYS_family1",
                                                  "CYS_family7",
                                                  "CYS_family28",
                                                  "CYS_family45", 
                                                  "CYS_family53",
                                                  "CYS_family24"))
p1 <- ggplot(df_all,aes(family,log10(abundance),fill = condition)) +  
  geom_boxplot(na.rm = F)+
  theme(aspect.ratio = 4/8)+
  theme(axis.text.x = element_blank())+
  labs(x = "",
       y = "log10(abundance)") +
  theme_bw() +
  scale_fill_lancet()+
  #theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 45)) +
  #stat_summary(fun = "mean",geom="point",shape=23,size=3,fill="white")+
  theme(text=element_text(size = 12))

p1
dev.off()

##prevalence###
differ_pre <- Sol_BE_cluster_IBD[Sol_BE_cluster_IBD$class == "differential Sol_BE", ]
differ_pre_1 <- differ_pre[, c(1,2,3)]
colnames(differ_pre_1) <- c("SOl_BE","nonIBD","IBD")
differ_pre_1$SOl_BE <- as.factor(differ_pre_1$SOl_BE)
df_pre <- differ_pre_1 %>% pivot_longer(cols = colnames(differ_pre_1)[2:3],
                                        names_to = "condition",
                                        values_to = "prevalence")

df_pre$condition <- factor(df_pre$condition,
                           levels = c("IBD",
                                      "nonIBD"))

df_pre$SOl_BE <- factor(df_pre$SOl_BE,levels = c("SDR_family73",
                                                 "CFAT_family3",
                                                 "CYS_family1",
                                                 "CYS_family7",
                                                 "CYS_family28",
                                                 "CYS_family45", 
                                                 "CYS_family53",
                                                 "CYS_family24"))
#pdf(file="prevalence.pdf",width = 10, height = 10)
pre <- ggplot(df_pre, aes(x=SOl_BE, y=prevalence, fill=condition)) + 
  geom_bar(stat="identity",
           width = 0.75,
           position=position_dodge()) +
  theme(aspect.ratio = 4/8)+
  labs(x = "",
       y = "Prevalence (%)") +
  theme_bw() +
  scale_fill_lancet()+
  theme(text=element_text(size = 12))
pre

library(patchwork)
pdf(file="prevalence_abundance_DF.pdf",width = 9, height = 5.5)  ## figure 1f
pre/p1 
dev.off()

#### Beta diversity ####
countdata_10 <- countdata_30[,colSums(countdata_30)!=0]
countdata <- as.data.frame(countdata_10)
ct <- as.data.frame(t(countdata))
library(tibble)
ct_1 <- tibble::rownames_to_column(ct,var="id")
ct_2 <-select(ct_1,c(1))
group_1<-metadata
group <- merge(group_1,ct_2,by="id")
write.table(group, file="mymeta_group.csv",sep = "\t",quote=F,row.names = F)

data <- countdata[,as.vector(group$id)]
data <- data[,colSums(data)>0]
distance <- vegdist(t(data), method = 'bray')
distance
pcoa <- cmdscale(distance, k = (nrow(t(data)) - 1), eig = TRUE)

# save the distance info
write.csv(as.matrix(distance), 'distance.csv', quote = F)
write.csv(as.matrix(pcoa$points), 'pcoa_value.csv', quote = F)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
# check info
summary(pcoa)
# check the coordinate value
pcoa$eig
point <- data.frame(pcoa$point)
# save the coordinate info
write.csv(point, 'pcoa.sample.csv')
# compute first two coordinate value
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
# extract first two coordinate value
sample_condition <- data.frame({pcoa$point})[1:2]
sample_condition$id <- rownames(sample_condition)
names(sample_condition)[1:2] <- c('PCoA1', 'PCoA2')

sample_condition <- merge(sample_condition, group, by = 'id', all.x = FALSE)
write.csv(sample_condition, 'sample_condition.csv', quote = F)
### PERMANOVA for conditions
dt <- data.frame(t(data))
write.table(data,file="data_f.tsv",quote = F)
group_1 <-read.table('mymeta_group.csv',header=T,sep="\t")

adonis_result  <- adonis(dt ~ condition, group_1, distance = 'bray', permutations = 999)
#adonis_result[is.na(adonis_result)] <- 0
# save result
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all_condition.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
# pairwise
if(F){
  group_name <- unique(group_1$condition)
  
  adonis_result_two <- NULL
  for (i in 1:(length(group_name) - 1)) {
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group_1, condition %in% c(group_name[i], group_name[j]))
      otu_ij <- dt[group_ij$id, ]
      adonis_result_otu_ij <- adonis(otu_ij~condition, group_ij, permutations = 999, distance = 'bray')     #随机置换检验 999 次
      adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
    }
  }
  adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
  names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
  
  # marker significance
  for (i in 1:nrow(adonis_result_two)) {
    if (adonis_result_two[i, 'Pr (>F)'] <= 0.001) adonis_result_two[i, 'Sig'] <- '***'
    else if (adonis_result_two[i, 'Pr (>F)'] <= 0.01) adonis_result_two[i, 'Sig'] <- '**'
    else if (adonis_result_two[i, 'Pr (>F)'] <= 0.05) adonis_result_two[i, 'Sig'] <- '*'
  }
  
  # save result
  write.table(adonis_result_two, 'PERMANOVA.result_pairwise_9999_all_conditions.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
}
pdf(file="Pcoa.pdf",width = 6,height = 6)
sample_condition$condition <- factor(sample_condition$condition, levels = c("IBD","nonIBD"))
ggscatter(sample_condition, x= "PCoA1", y = "PCoA2",
          color = "condition",
          #shape = "condition",
          ellipse = TRUE, 
          mean.point = TRUE, star.plot = T, 
          ellipse.level = 0.95, 
          ggtheme = theme_minimal()) +
  labs(title=" ",subtitle = "P-value= 0.001",
       x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'),
       y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme_classic()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_point(size = 4)+ 
  theme(panel.grid = element_line(color = 'black', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
  theme(axis.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 12),legend.position="bottom")+scale_color_npg()

dev.off()






