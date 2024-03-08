library(stringr)
library(tidyr)
library(openxlsx)
library(readxl)
library(ggplot2)
library(reshape2)
library(filenamer)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggradar)
library(tidyverse)
library(dplyr)
library(gg.gap)
library(ggsci)
library(plyr)
library(tibble)
rm(list=ls())

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
  
  Cor_heatmap <- function(df,metabolites){
    df_t <-t(df)
    cormat <- round(cor(df_t,method="spearman"),2)
    cormat <- reorder_cormat(cormat)
    melted_cormat <- melt(cormat, na.rm = TRUE)
    #lower_tri <- get_lower_tri(cormat)
    # Melt the correlation matrix
    #melted_cormat <- melt(lower_tri, na.rm = TRUE)
    
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0.8, limit = c(0.5,1), space = "Lab", 
                           name="Spearman\nCorrelation") + 
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
        #legend.position = c(0.8, 0.3),
        legend.position = "bottom",
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
    ggsave(ggheatmap,filename=paste("spearman_",metabolites,".pdf",sep=""),width = 6,height = 6)
  }
  Cor_heatmap_L <- function(df,metabolites){
    df_t <-t(df)
    cormat <- round(cor(df_t,method="spearman"),2)
    cormat <- reorder_cormat(cormat)
    lower_tri <- get_lower_tri(cormat)
    # Melt the correlation matrix
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0.8, limit = c(0.5,1), space = "Lab", 
                           name="spearman\nCorrelation") + 
      theme(axis.text = element_text(size=4)) +
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
    ggsave(ggheatmap,filename=paste("spearman_",metabolites,".pdf",sep=""),height = 16,width = 16)
  }
  
}
############################################################################
### 01 correlation with co-elutions in nature microbiology 2019 IBD ########
############################################################################
###correlation
metabolites_features <- read.csv("../metabolomics_ISF_v2/FT_FTC_data.tsv",sep="\t")
species_19nm_1 <- read_xlsx("./nature_microbiology_19/41564_2018_306_MOESM6_ESM.xlsx")

FT_FTC <- species_19nm_1[-c(1:8), ]
write.table(FT_FTC,file="species_FTC_data.tsv",sep="\t",quote=F,row.names = F)

#### differential metabolites ######
metabolite_1 <- metabolites_features[rownames(metabolites_features) %in% c("HILIC-neg_Cluster_1835","HILIC-neg_Cluster_1821","HILIC-neg_Cluster_1793","HILIC-pos_Cluster_2115","HILIC-pos_Cluster_2155"), ]

metabolite_t1 <- as.data.frame(t(metabolite_1))
FT_meta <- read.delim2("nature_microbiology_19/Supplementary Dataset 2- Per-subject metabolite relative abundance profiles")

FT_meta <- as.data.frame(t(FT_meta[c(1:7), ]))
colnames(FT_meta) <- FT_meta[1,]
FT_meta <- FT_meta[-1, ]
FT_meta <- FT_meta[2]

FT_mb_meta <- merge(FT_meta,metabolite_t1,by="row.names")
FT_mb_meta$Diagnosis <- str_replace_all(FT_mb_meta$Diagnosis,"Control","nonIBD")
if(T){
FT_mb_meta$Diagnosis <- str_replace_all(FT_mb_meta$Diagnosis,"UC","IBD")
FT_mb_meta$Diagnosis <- str_replace_all(FT_mb_meta$Diagnosis,"CD","IBD")
}
colnames(FT_mb_meta) <- str_replace_all(colnames(FT_mb_meta),"HILIC-","HILIC_")

FT_mb_meta$HILIC_pos_Cluster_2155

FT_mb_meta_w <- pivot_longer(FT_mb_meta,
                             cols=HILIC_neg_Cluster_1793:HILIC_pos_Cluster_2155,
                             names_to = "metabolites",
                             values_to = "abundance")

my_comparison <- list(c("nonIBD","IBD"))

if(T){
  p1 <- ggplot(FT_mb_meta_w,aes(Diagnosis,log2(abundance),fill = Diagnosis)) +  
    geom_boxplot(na.rm = F)+
    theme(axis.text.x = element_blank())+
    stat_compare_means(comparisons = my_comparison,
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       paired = FALSE)+
    labs(x = "",
         y = "log2(SoLs relative abundance)") +
    theme_bw(base_size = 12) +
    scale_fill_lancet()+
    #theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 45)) +
    theme(text=element_text(size = 12),
          legend.position = "right")+
    facet_wrap(~metabolites,scales="free",nrow = 1)
  
  p1
}


############################################################################
### 02. correlation with co-elutions in HMP2 ###############################
############################################################################
##metadata
metadata <- read.delim2("./2018-01-01078D-s6/hmp2_metadata.csv",sep=",")
metadata_0 <- metadata[,c(2:5)]
table(metadata_0$data_type)
MG <- metadata_0[metadata_0$data_type=="metagenomics", ]
MB <- metadata_0[metadata_0$data_type=="metabolomics", ]
colnames(MG)[1] <- c("MG_ID")
colnames(MB)[1] <- c("MB_ID")

MG_MB_pair <- merge(MB,MG,by="site_sub_coll")
MG_MB_pair_1 <- MG_MB_pair[!str_detect(MG_MB_pair$MG_ID,"_TR"), ]  ##remove technology duplicates
MG_MB_pair_2 <- MG_MB_pair_1[!str_detect(MG_MB_pair_1$MG_ID,"_P"), ] 
MG_MB_pair_3 <- MG_MB_pair_1[!MG_MB_pair_1$site_sub_coll %in% MG_MB_pair_2$site_sub_coll,  ]

MG_MB_pair_f <- rbind(MG_MB_pair_2,MG_MB_pair_3)
length(unique(MG_MB_pair_f$site_sub_coll))

###metabolomics
meta_hmp <- read.delim2("./nature_microbiology_19/HMP2_metabolomics.csv",sep=",")

###relative abundance
features_hmp_0 <- meta_hmp[,c(7:ncol(meta_hmp))]
rownames(features_hmp_0) <- features_hmp_0$Compound
features_hmp_0 <- features_hmp_0[-1]
features_hmp_0[] <- lapply(features_hmp_0,as.numeric)
features_hmp_0[is.na(features_hmp_0)] <- 0
HILn <- features_hmp_0[str_detect(row.names(features_hmp_0),"HILn"), ]
HILp <- features_hmp_0[str_detect(row.names(features_hmp_0),"HILp"), ]
C18n <- features_hmp_0[str_detect(row.names(features_hmp_0),"C18"), ]
C8p <- features_hmp_0[str_detect(row.names(features_hmp_0),"C8p"), ]

HILn_r <- as.data.frame(prop.table(as.matrix(HILn),2) * 1e6)
HILp_r <- as.data.frame(prop.table(as.matrix(HILp),2) * 1e6)
C18n_r <- as.data.frame(prop.table(as.matrix(C18n),2) * 1e6)
C8p_r <- as.data.frame(prop.table(as.matrix(C8p),2) * 1e6)

features_hmp_re_ab <- rbind(HILn_r,HILp_r,C18n_r,C8p_r)

features_hmp <- features_hmp_re_ab[ ,colnames(features_hmp_re_ab) %in% MG_MB_pair_f$MB_ID]


metabolite_hmp <- features_hmp[rownames(features_hmp) %in% c("HILp_QI21780","HILp_QI20990","HILn_QI8871","HILn_QI9663","HILn_QI11564"), ]
metabolite_hmp[is.na(metabolite_hmp)] <- 0

## metagenome_taxa ##
species_hmp <- read.delim2("./2018-01-01078D-s6/IBDMDB_taxonomic_profiles_3.tsv")
colnames(species_hmp) <- str_replace_all(colnames(species_hmp), "_profile","")
species_hmp_s <- species_hmp[str_detect(species_hmp$Feature.Sample,"s__"), ]

Species <- function(x){
  res <- strsplit(x,"s__")[[1]]
  paste0(res[2])
}

species_hmp_s$species <- unlist(lapply(species_hmp_s$Feature.Sample,Species))
rownames(species_hmp_s) <- species_hmp_s$species
species_hmp_s <- species_hmp_s[,c(2:1638)]
species_hmp_s_1 <- species_hmp_s[ ,colnames(species_hmp_s) %in% MG_MB_pair_f$MG_ID]
species_hmp_s_1[] <- lapply(species_hmp_s_1,as.numeric)

species_hmp_s_1 <- species_hmp_s_1[rowSums(species_hmp_s_1)!=0, ]
species_hmp_s_1 <- species_hmp_s_1[rowSums(species_hmp_s_1 >0) > ncol(species_hmp_s_1)*0.05, ]
species_hmp_s_2 <- as.data.frame(t(species_hmp_s_1))
MG_MB_0 <- MG_MB_pair_f[,c(2,5)]
species_hmp_s_3 <- merge(species_hmp_s_2,MG_MB_0,by.x = "row.names",by.y="MG_ID")
rownames(species_hmp_s_3) <- species_hmp_s_3$MB_ID
species_hmp_s_4 <- species_hmp_s_3[,c(2:166)]
species_hmp_sf <- as.data.frame(t(species_hmp_s_4))

metabolite_hmp_1 <- metabolite_hmp[,match(colnames(species_hmp_sf),colnames(metabolite_hmp))]

metab_sp_hmp <- rbind(metabolite_hmp_1,species_hmp_sf)

Cor_heatmap_L(metab_sp_hmp,"metab_sp_hmp_relative")

metab_sp_hmp_f <- metab_sp_hmp[rownames(metab_sp_hmp) %in% c("HILp_QI21780","HILp_QI20990","HILn_QI8871","HILn_QI9663","HILn_QI11564","Alistipes_putredinis","Alistipes_finegoldii"), ]

Cor_heatmap(metab_sp_hmp_f,"metab_sp_hmp_f_relative")

### differential analysis
##metabolites
metabolite_hmp_t <- as.data.frame(t(metabolite_hmp))

metadata_1 <- metadata[,c(2,5,71)]
metadata_1 <- metadata_1[metadata_1$data_type=="metabolomics", ]
metadata_1 <- metadata_1[metadata_1$External.ID %in% MG_MB_pair_f$MB_ID, ]
metadata_1 <- metadata_1[-2]

###merged
metabolite_hmp_meta <- merge(metabolite_hmp_t,metadata_1,by.x = "row.names",by.y = "External.ID")

metabolite_hmp_meta_w <- pivot_longer(metabolite_hmp_meta,
                                      cols = HILn_QI9663:HILp_QI21780,
                                      names_to = "metabolites",
                                      values_to = "abundance")
metabolite_hmp_meta_w$log_abundance <- log10(metabolite_hmp_meta_w$abundance)

metabolite_hmp_meta_w$diagnosis <- str_replace_all(metabolite_hmp_meta_w$diagnosis,"UC","IBD")
metabolite_hmp_meta_w$diagnosis <- str_replace_all(metabolite_hmp_meta_w$diagnosis,"CD","IBD")
table(metabolite_hmp_meta_w$diagnosis)
my_comparison <- list(c("nonIBD","IBD"))

if(T){
  p2 <- ggplot(metabolite_hmp_meta_w,aes(diagnosis,log2(abundance),fill = diagnosis)) +  
    geom_boxplot(na.rm = F)+
    theme(axis.text.x = element_blank())+
    stat_compare_means(comparisons = my_comparison,
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       paired = FALSE)+
    labs(x = "",
         y = "log2(SoLs relative abundance)") +
    theme_bw(base_size = 12) +
    scale_fill_lancet()+
    #theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 45)) +
    theme(text=element_text(size = 12),
          legend.position = "right")+
    facet_wrap(~metabolites,scales="free",nrow = 1)
  
  p2
}

library(patchwork)
pdf(file = "Differential analysis_SoLs relative abundance_nrow5.pdf",width = 11,height = 5)
p2/p1
dev.off()


### supplementary data 3 ####
features_hmp_re_ab$Sample <- rownames(features_hmp_re_ab)
metab_sp_hmp$features <- rownames(metab_sp_hmp)
meta_spec$features <- row.names(meta_spec)


data_frames_1 <- list("sheet1"= metabolite_hmp_meta,
                      "sheet2"= features_hmp_re_ab,
                      "sheet3"= metab_sp_hmp)

write.xlsx(data_frames_1,file="Supplementary Data 3 SoL metabolites collected from matabolomics dataset 1.xlsx")

data_frames_2 <- list("sheet1"= FT_mb_meta,
                      "sheet2"= meta_spec)

write.xlsx(data_frames_2,file="Supplementary Data 4 SoL metabolites collected from matabolomics dataset 2.xlsx")











