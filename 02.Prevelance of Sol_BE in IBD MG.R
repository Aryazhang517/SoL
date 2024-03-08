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
Sulfonolipid_BE <-read.table("IBD_MG_SolBE_raw_counts_f.csv",sep="\t",row.names = 1,header=T,encoding="UTF-8")
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
write.table(ftr.tpm_out,file = "IBD_p1p2_MG_TPM_all.tsv", sep="\t",quote = F, row.names = F)

### Relative abundance of clusters ###############
countdata_1 <-read.csv("IBD_TPM_rev_cluster90_pfam50_MG.csv",sep="\t",row.names = 1)
countdata_2 <-countdata_1[rowSums(countdata_1) != 0,]
countdata_2 <-countdata_2[,colSums(countdata_2) != 0]
countdata_3 <- countdata_2[rowSums(countdata_2>0)>ncol(countdata_2)*0.05,]
###metadata
metadata <- read.table("IBD_metadata",sep="\t",header=T,encoding="UTF-8")
colnames(metadata) <- c("id","condition")

###Fisher.exact test for prevalence###
Case <- metadata %>% filter(metadata$condition %in% c("IBD"))
Ctr <- metadata %>% filter(metadata$condition %in% c("nonIBD"))
countdata_Case <- countdata_3[,colnames(countdata_3) %in% Case$id]
countdata_Case <-countdata_Case[,colSums(countdata_Case) != 0]
countdata_Ctr <- countdata_3[,colnames(countdata_3) %in% Ctr$id]
countdata_Ctr <-countdata_Ctr[,colSums(countdata_Ctr) != 0]
### prevalence
countdata_Case$percent <- (rowSums(countdata_Case>0)/ncol(countdata_Case))*100
countdata_Ctr$percent <- (rowSums(countdata_Ctr>0)/ncol(countdata_Ctr))*100
Sol_BE_perc_IBD <- countdata_Case %>% select(c(percent)) %>% setNames(c("percent_IBD"))
Sol_BE_perc_Ctr <- countdata_Ctr %>% select(c(percent)) %>% setNames(c( "percent_Ctr"))
Sol_BE_cluster_IBD <- merge(Sol_BE_perc_Ctr,Sol_BE_perc_IBD,by="row.names")
colnames(Sol_BE_cluster_IBD)[1] <-c("SOl_BE")
#####Fisher.exact test######
df_ctr <- countdata_Ctr[ ,c(1:2)]
df_ctr$ctr <- rowSums(countdata_Ctr[ ,1:(ncol(countdata_Ctr))]>0)
df_ctr$nonctr <- (ncol(countdata_Ctr))- df_ctr$ctr
df_ctr <- df_ctr[ ,c(3,4)]

df_case <- countdata_Case[ ,c(1:2)]
df_case$case <- rowSums(countdata_Case[ ,1:(ncol(countdata_Case))]>0)
df_case$noncase <-  (ncol(countdata_Case))- df_case$case
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

all_stat <- merge(Fisher_2,Sol_BE_cluster_IBD,by.x="cluster",by.y="SOl_BE")
colnames(all_stat)[1] <- "family_id"
write.table(all_stat, file = 'Sol_enzyme_prevalence_abundance_all_test.tsv',sep ="\t", row.names=F,quote = F)

######### Beta-diversity ########
countdata_10 <- countdata_3[,colSums(countdata_3)!=0]
countdata <- as.data.frame(countdata_10)
ct <- as.data.frame(t(countdata))
ct_1 <- tibble::rownames_to_column(ct,var="id")
ct_2 <-select(ct_1,c(1))
group_1<-metadata
group <- merge(group_1,ct_2,by="id")
#write.table(group, file="mymeta_group.csv",sep = "\t",quote=F,row.names = F)
data <- countdata[,as.vector(group$id)]
data <- data[,colSums(data)>0]
distance <- vegdist(t(data), method = 'jaccard')
distance
pcoa <- cmdscale(distance, k = (nrow(t(data)) - 1), eig = TRUE)
# save the distance info
write.csv(as.matrix(distance), 'distance.csv', quote = F)
write.csv(as.matrix(pcoa$points), 'pcoa_value.csv', quote = F)

ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
# check info
summary(pcoa)
pcoa$eig
point <- data.frame(pcoa$point)
write.csv(point, 'pcoa.sample.csv')

pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_condition <- data.frame({pcoa$point})[1:2]
sample_condition$id <- rownames(sample_condition)
names(sample_condition)[1:2] <- c('PCoA1', 'PCoA2')

sample_condition <- merge(sample_condition, group, by = 'id', all.x = FALSE)
write.csv(sample_condition, 'sample_condition.csv', quote = F)
### PERMANOVA for conditions
dt <- data.frame(t(data))
write.table(data,file="data_f.tsv",quote = F)
group_1 <-read.table('mymeta_group.csv',header=T,sep="\t")

adonis_result  <- adonis(dt ~ condition, group_1, distance = 'jaccard', permutations = 999)
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all_condition.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
##plot
pdf(file="Pcoa_IBD_MG.pdf",width = 6,height = 6)
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
  theme(panel.grid = element_line(color = 'black', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
  theme(axis.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 12),legend.position="bottom")+scale_color_npg()
dev.off()






