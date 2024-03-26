library(stringr)
library(tidyr)
library(openxlsx)
library(readxl)
library(ggplot2)
library(reshape2)
library(filenamer)
library(plyr)
library(tibble)
library(Hmisc)

rm(list = ls())
if(T){
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  flat_cor_mat <- function(cor_r, cor_p){
    cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
    cor_r <- gather(cor_r, column, cor, -1)
    cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
    cor_p <- gather(cor_p, column, p, -1)
    cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
    cor_p_matrix
  }
  
  Cor_heatmap <- function(df,metabolites){
    df_t <-t(df)
    ## check correlation and related p-values ##
    cormat_p <- rcorr(as.matrix(df_t),type = c("pearson"))
    melted_cormat_p <- flat_cor_mat(cormat_p$r, cormat_p$P)
    
    melted_cormat_p$stars <- " "
    melted_cormat_p[which(0.05 > melted_cormat_p$p), "stars"]="*"
    melted_cormat_p[which(0.01 > melted_cormat_p$p), "stars"]="**"
    melted_cormat_p[which(0.001 > melted_cormat_p$p), "stars"]="***"
    write.xlsx(melted_cormat_p,file = paste0("pearson_cor_pvalue",metabolites,".xlsx"))

    ## visualization 
    cormat <- round(cor(df_t,method="pearson"),2)
    cormat <- reorder_cormat(cormat)
    lower_tri <- get_lower_tri(cormat)
    # Melt the correlation matrix
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0.8, limit = c(0.7,1), space = "Lab", 
                           name="Pearson\nCorrelation") + 
      theme(axis.text = element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 90)) +
      #sprintf("%0.2f", round(a, digits = 2))))
      geom_text(aes(Var2, Var1, label = sprintf("%0.1f", round(value, digits = 2))), color = "black", size = 2) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.8, 0.3),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
    ggsave(ggheatmap,filename=paste("pearson_",metabolites,".pdf",sep=""))
  }
  Cor_heatmap_L <- function(df,metabolites){
    df_t <-t(df)
    
    ## check correlation and related p-values ##
    cormat_p <- rcorr(as.matrix(df_t),type = c("pearson"))
    melted_cormat_p <- flat_cor_mat(cormat_p$r, cormat_p$P)
    
    melted_cormat_p$stars <- " "
    melted_cormat_p[which(0.05 > melted_cormat_p$p), "stars"]="*"
    melted_cormat_p[which(0.01 > melted_cormat_p$p), "stars"]="**"
    melted_cormat_p[which(0.001 > melted_cormat_p$p), "stars"]="***"
    write.xlsx(melted_cormat_p,file = paste0("pearson_cor_pvalue",metabolites,".xlsx"))

    ## visualization 
    cormat <- round(cor(df_t,method="pearson"),2)
    cormat <- reorder_cormat(cormat)
    lower_tri <- get_lower_tri(cormat)
    # Melt the correlation matrix
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0.8, limit = c(0.7,1), space = "Lab", 
                           name="Pearson\nCorrelation") + 
      theme(axis.text = element_text(size=2.5)) +
      theme(axis.text.x = element_text(angle = 90)) +
      #sprintf("%0.2f", round(a, digits = 2))))
      #geom_text(aes(Var2, Var1, label = sprintf("%0.1f", round(value, digits = 2))), color = "black", size = 2) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.8, 0.3),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
    ggsave(ggheatmap,filename=paste("pearson_",metabolites,".pdf",sep=""),height = 16,width = 16)
  }
  
}
load("NM_HMP2_metatable_feature_table.Rdata")
excel_sheets("20230708_SoL_potential_candidates_non_strictly.xlsx")
############################################################################
### 01 correlation with co-elutions in nature microbiology 2019 IBD ########
############################################################################
SoL_NM_neg <- read_excel("20230708_SoL_potential_candidates_non_strictly.xlsx", sheet =2)
SoL_NM_pos <- read_excel("20230708_SoL_potential_candidates_non_strictly.xlsx",sheet = 3)

SoL_NM_0 <- rbind(SoL_NM_neg,SoL_NM_pos)
SoL_NM_1 <- SoL_NM_0[ ,c(1:3,8:12,4:7)]
###ppm =< 5
SoL_NM_2 <- SoL_NM_1[(5 >= SoL_NM_1$ppm & SoL_NM_1$ppm >= -5), ]
###featurecounts
FT_meta$Retention.Time <- as.numeric(FT_meta$Retention.Time)
HILIC_neg_meta <- FT_meta[str_detect(FT_meta$Metabolomic.Feature,"HILIC-neg"), ]
HILIC_pos_meta <- FT_meta[str_detect(FT_meta$Metabolomic.Feature,"HILIC-pos"), ]
C8_pos_meta <- FT_meta[str_detect(FT_meta$Metabolomic.Feature,"C8-pos"), ]
###correlation
FT_FTC <- FT_table[-c(1:7), ]
rownames(FT_FTC) <- FT_FTC$X..Feature...Sample
FT_FTC <- FT_FTC[-1]
#FT_FTC_1 <- FT_FTC[, str_detect(colnames(FT_FTC),"PRISM")]
write.table(FT_FTC,file="FT_FTC_data.tsv",sep="\t",quote=F)
FT_FTC_0 <- read.csv("FT_FTC_data.tsv",sep="\t")

HILIC_neg_1793_1821_1835 <- HILIC_neg_meta[(HILIC_neg_meta$Retention.Time >= 2.71 &  HILIC_neg_meta$Retention.Time <= 2.85),]
HILIC_pos_2155_2115 <- HILIC_pos_meta[(HILIC_pos_meta$Retention.Time >= 1.54 & HILIC_pos_meta$Retention.Time <= 1.65), ]
C8_pos_1239 <- C8_pos_meta[(C8_pos_meta$Retention.Time >= 6.88 & C8_pos_meta$Retention.Time <= 6.98), ]
C8_pos_1177 <- C8_pos_meta[(C8_pos_meta$Retention.Time >= 7.44 & C8_pos_meta$Retention.Time <= 7.54), ]
C8_pos_1188 <- C8_pos_meta[(C8_pos_meta$Retention.Time >= 6.22 & C8_pos_meta$Retention.Time <= 6.32), ]
C8_pos_1296 <- C8_pos_meta[(C8_pos_meta$Retention.Time >= 7.16 & C8_pos_meta$Retention.Time <= 7.26), ]

HILIC_neg_1793_1821_1835_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% HILIC_neg_1793_1821_1835$Metabolomic.Feature, ]
HILIC_pos_2155_2115_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% HILIC_pos_2155_2115$Metabolomic.Feature, ]
C8_pos_1239_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% C8_pos_1239$Metabolomic.Feature, ]
C8_pos_1177_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% C8_pos_1177$Metabolomic.Feature, ]
C8_pos_1188_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% C8_pos_1188$Metabolomic.Feature, ]
C8_pos_1296_pre <- FT_FTC_0[rownames(FT_FTC_0) %in% C8_pos_1296$Metabolomic.Feature, ]

Cor_heatmap(HILIC_neg_1793_1821_1835_pre,"HILIC_neg_1793_1821_1835")
Cor_heatmap_L(HILIC_pos_2155_2115_pre,"HILIC_pos_2155_2115")
Cor_heatmap(C8_pos_1239_pre,"C8_pos_1239")
Cor_heatmap(C8_pos_1177_pre,"C8_pos_1177")
Cor_heatmap(C8_pos_1188_pre,"C8_pos_1188")
Cor_heatmap(C8_pos_1296_pre,"C8_pos_1296")

############################################################################
### 02. correlation with co-elutions in HMP2 ###############################
############################################################################
Sol_hmp_neg <- read_excel("20230708_SoL_potential_candidates_non_strictly.xlsx",sheet = 4)
Sol_hmp_pos <- read_excel("20230708_SoL_potential_candidates_non_strictly.xlsx",sheet = 5)
Sol_hmp_0 <- rbind(Sol_hmp_neg,Sol_hmp_pos)
Sol_hmp_1 <- Sol_hmp_0[ ,c(7,3,4,1,8:11,5,6)]
SoL_hmp_2 <- Sol_hmp_1[(Sol_hmp_1$ppm <= 5 & Sol_hmp_1$ppm >= -5), ]
## counts table ##
hmp_FTC <- HMP2[,-c(1:6)]
rownames(hmp_FTC) <- hmp_FTC$Compound
hmp_FTC <- hmp_FTC[-1]
write.table(hmp_FTC,file="hmp_FTC_data.tsv",sep="\t",quote=F)
hmp_FTC_0 <- read.csv("hmp_FTC_data.tsv",sep="\t")
hmp_FTC_0[is.na(hmp_FTC_0)] <- 0

HMP2_meta$RT <- as.numeric(HMP2_meta$RT)
hmp_HILIC_neg_meta <- HMP2_meta[str_detect(HMP2_meta$Method,"HILIC-neg"), ]
hmp_C18_neg_meta <- HMP2_meta[str_detect(HMP2_meta$Method,"C18-neg"), ]
hmp_C8_pos_meta <- HMP2_meta[str_detect(HMP2_meta$Method,"C8-pos"), ]
hmp_HILIC_pos_meta <- HMP2_meta[str_detect(HMP2_meta$Method,"HILIC-pos"), ]

HILn_QI18296_QI28227 <- hmp_HILIC_neg_meta[(hmp_HILIC_neg_meta$RT >= 2.56 & hmp_HILIC_neg_meta$RT <= 2.68), ]
HILn_QI9934 <- hmp_HILIC_neg_meta[(hmp_HILIC_neg_meta$RT >= 2.92 & hmp_HILIC_neg_meta$RT <= 3.02), ]
C18n_QI12075 <- hmp_C18_neg_meta[(hmp_C18_neg_meta$RT >= 18.23 & hmp_C18_neg_meta$RT <= 18.32), ]
HILp_QI20990_QI21780_QI19947 <- hmp_HILIC_pos_meta[(hmp_HILIC_pos_meta$RT >= 1.38 & hmp_HILIC_pos_meta$RT <= 1.49), ]
HILp_QI19946 <- hmp_HILIC_pos_meta[(hmp_HILIC_pos_meta$RT >= 1.48 & hmp_HILIC_pos_meta$RT <= 1.58), ]
HILp_QI21779_QI20989 <- hmp_HILIC_pos_meta[(hmp_HILIC_pos_meta$RT >= 1.65 & hmp_HILIC_pos_meta$RT <= 1.76), ]
C8p_QI18677 <- hmp_C8_pos_meta[(hmp_C8_pos_meta$RT >= 6.68 & hmp_C8_pos_meta$RT <= 6.78), ]
C8p_QI17936_QI15640 <- hmp_C8_pos_meta[(hmp_C8_pos_meta$RT >= 6.98 & hmp_C8_pos_meta$RT <= 7.09), ]
C8p_QI16404 <- hmp_C8_pos_meta[(hmp_C8_pos_meta$RT >= 7.28 & hmp_C8_pos_meta$RT <= 7.38), ]
C8p_QI17160_QI18372 <- hmp_C8_pos_meta[(hmp_C8_pos_meta$RT >= 7.52 & hmp_C8_pos_meta$RT <= 7.62), ]

HILn_QI18296_QI28227_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% HILn_QI18296_QI28227$Compound, ]
HILn_QI9934_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% HILn_QI9934$Compound, ]
C18n_QI12075_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% C18n_QI12075$Compound, ]
HILp_QI20990_QI21780_QI19947_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% HILp_QI20990_QI21780_QI19947$Compound, ]
HILp_QI19946_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% HILp_QI19946$Compound, ]
HILp_QI21779_QI20989_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% HILp_QI21779_QI20989$Compound, ]
C8p_QI18677_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% C8p_QI18677$Compound, ]
C8p_QI17936_QI15640_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% C8p_QI17936_QI15640$Compound, ]
C8p_QI16404_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% C8p_QI16404$Compound, ]
C8p_QI17160_QI18372_pre <- hmp_FTC_0[rownames(hmp_FTC_0) %in% C8p_QI17160_QI18372$Compound, ]

Cor_heatmap_L(HILn_QI18296_QI28227_pre,"HILn_QI18296_QI28227")
Cor_heatmap_L(HILn_QI9934_pre,"HILn_QI9934")
Cor_heatmap(C18n_QI12075_pre,"C18n_QI12075")
Cor_heatmap_L(HILp_QI20990_QI21780_QI19947_pre,"HILp_QI20990_QI21780_QI19947")
Cor_heatmap_L(HILp_QI19946_pre,"HILp_QI19946")
Cor_heatmap_L(HILp_QI21779_QI20989_pre,"HILp_QI21779_QI20989") #
Cor_heatmap_L(C8p_QI18677_pre,"C8p_QI18677")
Cor_heatmap_L(C8p_QI17936_QI15640_pre,"C8p_QI17936_QI15640")
Cor_heatmap_L(C8p_QI16404_pre,"C8p_QI16404")
Cor_heatmap_L(C8p_QI17160_QI18372_pre,"C8p_QI17160_QI18372")

save(SoL_NM_2,SoL_hmp_2,FT_meta,FT_FTC_0,HMP2_meta,hmp_FTC_0,
     file="001_Correlation_related_raw_data_pearson.Rdata")

############################################################################
######### Plot for specific candidates ######
############################################################################
###HILIC_neg_1793_1821_1835
HILIC_neg_1793_1821_1835_f <- FT_FTC_0[rownames(FT_FTC_0) %in% c("HILIC-neg_Cluster_1835","HILIC-neg_Cluster_1821","HILIC-neg_Cluster_1793"), ]
###HILIC_pos_2155_2115
HILIC_pos_2155_2115_f <- FT_FTC_0[rownames(FT_FTC_0) %in% c("HILIC-pos_Cluster_2115","HILIC-pos_Cluster_2155"), ]
#HILn_QI18296_QI28227
HILn_QI8871_QI11564_f <- hmp_FTC_0[rownames(hmp_FTC_0) %in% c("HILn_QI8871","HILn_QI9663","HILn_QI11564","HILn_QI15355","HILn_QI28227","HILn_QI21386"), ]
#HILp_QI20990_QI21780_QI19947
HILp_QI20990_QI21780_QI19947_f <- hmp_FTC_0[rownames(hmp_FTC_0) %in% c("HILp_QI21780","HILp_QI20990","HILp_QI22371","HILp_QI19947"), ]
#C8p_QI17936_QI15640
C8p_QI17936_QI15640_f <- hmp_FTC_0[rownames(hmp_FTC_0) %in% c("C8p_QI15640","C8p_QI17936","C8p_QI19098"), ]
#C8p_QI17160_QI18372
C8p_QI17160_QI18372_f <- hmp_FTC_0[rownames(hmp_FTC_0) %in% c("C8p_QI18372","C8p_QI17160","C8p_QI22700"),]

Cor_heatmap(HILIC_neg_1793_1821_1835_f,"HILIC_neg_1793_1821_1835_f")
Cor_heatmap(HILIC_pos_2155_2115_f,"HILIC_pos_2155_2115_f")
Cor_heatmap(HILn_QI8871_QI11564_f,"HILn_QI8871_QI11564_f")
Cor_heatmap(HILp_QI20990_QI21780_QI19947_f,"HILp_QI20990_QI21780_QI19947_f")
Cor_heatmap(C8p_QI17936_QI15640_f,"C8p_QI17936_QI15640_f")
Cor_heatmap(C8p_QI17160_QI18372_f,"C8p_QI17160_QI18372_f")





