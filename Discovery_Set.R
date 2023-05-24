###############################################################
# Code for analysis of MIS-C transcriptomic data
# Identification of diagnostic signature
###############################################################

library(stringr)
library(DESeq2)
library(biomaRt)
library(lubridate)
library(tidyverse)
library(sva)
library(gtools)
library(pROC)
library(nnet)
library(tidyverse)
library(cowplot)


###############################################################
# functions
###############################################################

deseq.results <- function(res, gene.names, outfile){
  res <- res[order(res$padj),]
  res.df <- data.frame(res)
  rownames(res.df) <- rownames(res)
  
  gene.names <- gene.names[match(rownames(res.df), gene.names$ensembl_gene_id),]
  res.df$gene <- gene.names$external_gene_name
  res.df$Significant <- 0
  res.df$Significant[res.df$padj < 0.05] <- 1
  
  
  # create volcano plot column
  res.df$Color <- 'black'
  res.df$Color[abs(res.df$log2FoldChange) > 1] <- 'limegreen'
  res.df$Color[res.df$padj < 0.05] <- 'gold'
  res.df$Color[res.df$padj < 0.05 & 
                 abs(res.df$log2FoldChange) > 1] <- 'red'
  
  res.df$Color <- factor(res.df$Color, levels = c("black", 'limegreen', "gold", 'red'))
  
  res.df <- res.df[!(is.na(res.df$log2FoldChange)),]
  res.df <- res.df[!(is.na(res.df$padj)),]
  print(range(res.df$log2FoldChange))
  write.csv(file = paste(outfile, '.csv', sep = ""), res.df)
  return(res.df)
}

roc_plot <- function(roc, title, colour){
  ci <- ci.se(roc, specificities=seq(0, 1, l=25))
  ci <- data.frame(x = as.numeric(rownames(ci)),
                   lower = ci[, 1],
                   upper = ci[, 3])
  p <- ggroc(roc) + theme_bw() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")+
    labs(x = 'Specificity', y = "Sensitivity")+
    geom_ribbon(
      data = ci,
      aes(x = x, ymin = lower, ymax = upper),
      alpha = 0.3,
      fill = colour,
      inherit.aes = F) +
    annotate("text", x=0.48, y=0.1, label= paste("AUC: ", 
                                                 round(auc(roc), 3)*100, 
                                                 "% (95% CI: ",
                                                 round(ci(roc)[1], 3)*100,
                                                 '%-', 
                                                 round(ci(roc)[3], 3)*100, "%)",
                                                 sep = ''))+
    theme(plot.title = element_text(size=12))+
    ggtitle(title)
}

###############################################################
# Differential expression analysis
###############################################################

# Load normalised data and name dds_samples (DESeq2 object)

# update design
design(dds_samples) <- model.matrix(~ 0 + age + phenotype + sex + rna_batch, data = data.frame(colData(dds_samples)))

# run DESeq2
dds_samples <- DESeq(dds_samples, 
                     betaPrior=FALSE, 
                     minReplicatesForReplace = Inf)

dds_samples <- dds_samples[which(mcols(dds_samples)$betaConv),]

norm <- counts(dds_samples, normalized=T)

# Perform differential expression analysis
# MIS-C vs DB
misc_db <- results(dds_samples, 
                   independentFiltering=TRUE,
                   alpha=0.05, 
                   pAdjustMethod="BH", 
                   parallel=FALSE, 
                   contrast=list(c("phenotypeMIS_C"), c("phenotypeDB")))
misc_db <- deseq.results(misc_db, gene.names = gene_names, 'SDE_Genes_MIS_C_vs_DB')

table(misc_db$Significant)
dim(misc_db[misc_db$Significant==1 & misc_db$log2FoldChange > 0,])
dim(misc_db[misc_db$Significant==1 & misc_db$log2FoldChange < 0,])

# MIS_C vs DV
misc_dv <- results(dds_samples, 
                   independentFiltering=TRUE,
                   alpha=0.05, 
                   pAdjustMethod="BH", 
                   parallel=FALSE, 
                   contrast=list(c("phenotypeMIS_C"), c("phenotypeDV")))
misc_dv <- deseq.results(misc_dv, gene.names = gene_names, 'SDE_Genes_MIS_C_vs_DV')

table(misc_dv$Significant)
dim(misc_dv[misc_dv$Significant==1 & misc_dv$log2FoldChange > 0,])
dim(misc_dv[misc_dv$Significant==1 & misc_dv$log2FoldChange < 0,])

# MIS_C vs KD
misc_kd <- results(dds_samples, 
                   independentFiltering=TRUE,
                   alpha=0.05, 
                   pAdjustMethod="BH", 
                   parallel=FALSE, 
                   contrast=list(c("phenotypeMIS_C"), c("phenotypeKD")))
misc_kd <- deseq.results(misc_kd, gene.names = gene_names, 'SDE_Genes_MIS_C_vs_KD')

table(misc_kd$Significant)
dim(misc_kd[misc_kd$Significant==1 & misc_kd$log2FoldChange > 0,])
dim(misc_kd[misc_kd$Significant==1 & misc_kd$log2FoldChange < 0,])

# MIS_C vs rest
misc_rest <- results(dds_samples,
                     independentFiltering=TRUE,
                     #cooksCutoff = 29,
                     alpha=0.05,
                     pAdjustMethod="BH",
                     parallel=FALSE,
                     contrast=list(c("phenotypeMIS_C"), c("phenotypeKD", "phenotypeDB", "phenotypeDV")),
                     listValues=c(1, -1/3))
misc_rest <- deseq.results(misc_rest, gene.names = gene_names, 'SDE_Genes_MIS_C_vs_Rest')

table(misc_rest$Significant)
dim(misc_rest[misc_rest$Significant==1 & misc_rest$log2FoldChange > 0,])
dim(misc_rest[misc_rest$Significant==1 & misc_rest$log2FoldChange < 0,])

###############################################################
# Signature discovery
###############################################################

# run FS-PLS to identify a signature that can distinguish MIS-C from DB, DV and KD but using genes SDE in MIS-C vs DB, MIS-C vs DV, MIS-C vs KD together 
misc_db_shortlist <- misc_db[abs(misc_db$log2FoldChange) > 0.5 & misc_db$padj < 0.001,]
misc_dv_shortlist <- misc_dv[abs(misc_dv$log2FoldChange) > 0.5 & misc_dv$padj < 0.001,]
misc_kd_shortlist <- misc_kd[abs(misc_kd$log2FoldChange) > 0.5 & misc_kd$padj < 0.001,]

# subset gene lists so significant genes only but also subset based on 'base mean' 
norm_samples <- counts(dds_samples, normalized = T)
norm_samples <- norm_samples[,!colnames(norm_samples) %in% dds_samples$sample_ID[dds_samples$phenotype=="HC"]]

misc_mean <- rowMeans(norm_samples[,colnames(norm_samples) %in% dds_samples$sample_ID[dds_samples$phenotype=="MIS_C"]])
misc_mean_50 <- misc_mean[misc_mean > 50]
misc_mean <- misc_mean[misc_mean > 100]

db_mean <- rowMeans(norm_samples[,colnames(norm_samples) %in% dds_samples$sample_ID[dds_samples$phenotype=="DB"]])
db_mean_50 <- db_mean[db_mean > 50]
db_mean <- db_mean[db_mean > 100]

dv_mean <- rowMeans(norm_samples[,colnames(norm_samples) %in% dds_samples$sample_ID[dds_samples$phenotype=="DV"]])
dv_mean_50 <- dv_mean[dv_mean > 50]
dv_mean <- dv_mean[dv_mean > 100]

kd_mean <- rowMeans(norm_samples[,colnames(norm_samples) %in% dds_samples$sample_ID[dds_samples$phenotype=="KD"]])
kd_mean_50 <- kd_mean[kd_mean > 50]
kd_mean <- kd_mean[kd_mean > 100]

genes_to_include <- unique(c(names(misc_mean), names(db_mean), names(dv_mean), names(kd_mean)))

genes_to_include <- genes_to_include[genes_to_include %in% intersect(intersect(intersect(names(db_mean_50), names(misc_mean_50)),
                                                                               names(dv_mean_50)), 
                                                                     names(kd_mean_50))]

genes_to_include <- genes_to_include[genes_to_include %in% 
                                       c(rownames(misc_db_shortlist), rownames(misc_dv_shortlist), rownames(misc_kd_shortlist))]

length(genes_to_include)

# set path to FS-PLS
options("fspls.path"="User/Path/to/FSPLS")

###STANDARD FSLPS OPTIONS
options("fspls.outfile"="./outfile.txt")
options("fspls.install"=FALSE)
options("fspls.numreps"=1)
options("method" = "fspls")
options("fspls.lambda" = "0") # #lambda.1se ##aa$lambda.min       ## if 0, no shrinkage. IF specified, then uses that shrinkage, if NULL, then uses CV
options("fspls.lambda1" = 1)      ## for pvalue adjustment
options("fspls.debug" = "off")     ## a lot of debugging information printed
options("fspls.log" = NULL)
options("fspls.applySidakCorrection"=FALSE)
options("fspls.inclLast"=FALSE)
options("fspls.useBIC"=F)
options("fspls.max"=10)
options("fspls.lambda" = 0) 
options("fspls.tikhonov"=F) ## set to F, using T is experimental
options("fspls.lars"=T)
options("fspls.alpha"=0)
options("fspls.coeff_se_reduce_percentile"=NA)  ##set to NA for no effect, or 0.5. Anything less than 0.5 shrinks
options("fspls.stopping"="auc")

#SOURCE CODE REQUIRED FROM FSPLS FOLDER
modules = c("fspls_lars_multi_looc1.R","extra_functions.R")
path.fspls = paste(getOption("fspls.path","."),modules,sep="/")
for(jk1 in path.fspls) source(jk1) 

##READ PHENOTYPES
# remove HC from dds_samples
dds_samples_fspls <- dds_samples[,dds_samples$phenotype != "HC"]
phen <- data.frame(colData(dds_samples_fspls))

# set up the phen object with phenotypic information of interest
phen$MIS_C_Batch <- as.character(phen$phenotype)
phen$MIS_C_Batch[phen$rna_batch=="New" & phen$phenotype != "MIS_C"] <- "Rest_New"
phen$MIS_C_Batch[phen$rna_batch=="A" & phen$phenotype != "MIS_C"] <- "Rest_A"
phen$rest <- 'rest'
phen$rest[phen$phenotype=="MIS_C"] <- 'MIS_C'

levs = c(levels(factor(phen$phenotype)), 
         levels(factor(phen$MIS_C_Batch)), 
         levels(factor(phen$rest)))

colnames(phen)[colnames(phen) == "phenotype"] <- 'disease_class'

all_vars1 = list(
  "MISCvsRestAll" = list("rest"=c("rest", "MIS_C")),
  "MISCvsRestOld" = list('MIS_C_Batch'=c("Rest_A", "MIS_C")), 
  "MISCvsRestNew" = list('MIS_C_Batch'=c( "Rest_New", "MIS_C")),
  "MISCvsKD"= list("disease_class" =  c("KD", "MIS_C")),
  "MISCvsDB"=  list("disease_class"=c("DB", "MIS_C")),
  "MISCvsDV" = list("disease_class"=c('DV', "MIS_C")))

normalise = NULL #"log"
dat_bm_50 <- norm_samples[genes_to_include,]

genes <- gene_names
phen = phen[,c(3, 7, 24:25)]

# log data
dat_bm_50 <- log2(data.frame(t(dat_bm_50))+1)

colnames(phen) <- c('Sequencing_Sample_ID','disease_class', 'MIS_C_Batch', "rest")

rownames(phen) <- phen$Sequencing_Sample_ID
phen <- phen[match(rownames(dat_bm_50), phen$Sequencing_Sample_ID),]

set.seed(101)
# create object that FSPLS can read
res_bm_50 = readRapidsData(all_vars1, dat_bm_50,   phen=phen, rand_inds =NULL);
# shuffle data
rand_inds_50 = attr(res_bm_50[[1]],"rand_inds")  ## HOW THE DATA WAS SHUFFLED

# run FSPLS
fspls_101 = .runFSPLS(res_bm_50, max=10, beam=20, batch=10, misclass=F)

# extract results
fspls_101$ggps$MISCvsKD$annotM$value <- str_sub(fspls_101$ggps$MISCvsKD$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsDB$annotM$value <- str_sub(fspls_101$ggps$MISCvsDB$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsDV$annotM$value <- str_sub(fspls_101$ggps$MISCvsDV$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsRestNew$annotM$value <- str_sub(fspls_101$ggps$MISCvsRestNew$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsRestOld$annotM$value <- str_sub(fspls_101$ggps$MISCvsRestOld$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsRestAll$annotM$value <- str_sub(fspls_101$ggps$MISCvsRestAll$annotM$value, - 5, - 1)
fspls_101$ggps$MISCvsTSS$annotM$value <- str_sub(fspls_101$ggps$MISCvsTSS$annotM$value, - 5, - 1)

# identify top model
top_models_101 <- rbind(DB = fspls_101$ggps$MISCvsDB$annotM[which(fspls_101$ggps$MISCvsDB$annotM$value == max(fspls_101$ggps$MISCvsDB$annotM$value)),],
                        DV = fspls_101$ggps$MISCvsDV$annotM[which(fspls_101$ggps$MISCvsDV$annotM$value == max(fspls_101$ggps$MISCvsDV$annotM$value)),],
                        KD = fspls_101$ggps$MISCvsKD$annotM[which(fspls_101$ggps$MISCvsKD$annotM$value == max(fspls_101$ggps$MISCvsKD$annotM$value)),],
                        New = fspls_101$ggps$MISCvsRestNew$annotM[which(fspls_101$ggps$MISCvsRestNew$annotM$value == max(fspls_101$ggps$MISCvsRestNew$annotM$value)),],
                        All = fspls_101$ggps$MISCvsRestOld$annotM[which(fspls_101$ggps$MISCvsRestOld$annotM$value == max(fspls_101$ggps$MISCvsRestOld$annotM$value)),], 
                        REST_ALL = fspls_101$ggps$MISCvsRestAll$annotM[which(fspls_101$ggps$MISCvsRestAll$annotM$value == max(fspls_101$ggps$MISCvsRestAll$annotM$value)),])

beta_df_101 <- data.frame(ID = str_extract_all(top_models_101$gene_nme1, 'ENSG...........')[which(top_models_101$size==max(as.numeric(top_models_101$size)))[1]])
colnames(beta_df_101) <- 'ID'
beta_df_101$gene <- genes$external_gene_name[match(beta_df_101$ID, genes$ensembl_gene_id)]

###############################################################
# evaluate performance of the 5 genes
###############################################################

# combine the signature with the top gene
genes_signature <- c(beta_df_101$ID, 
                     rownames(misc_rest)[1])

# extract counts
norm_samples_sig <- counts(dds_samples, normalized = T)
norm_samples_sig <- data.frame(t(norm_samples_sig[rownames(norm_samples_sig) %in% genes_signature,]))

norm_samples_sig$group <- dds_samples$phenotype

norm_samples_sig_hc <- norm_samples_sig

norm_samples_sig <- norm_samples_sig[norm_samples_sig$group != "HC",]
norm_samples_sig[,1:5] <- log2(norm_samples_sig[,1:5]+1)

norm_samples_sig$MISC <- 0
norm_samples_sig$MISC[norm_samples_sig$group=="MIS_C"] <- 1

# extract gene coefficients for each gene
mod <- summary(glm(MISC ~., data = norm_samples_sig[,c(1:5,7)]))$coefficients[,1]
mod <- mod[2:length(mod)]

# calculate weighted DRS
norm_samples_sig$DRS <- as.vector(as.matrix(norm_samples_sig[,match(names(mod), 
                                                                    colnames(norm_samples_sig))]) %*% as.matrix(mod))

# ROC curve analysis
roc_5_all <- roc(norm_samples_sig$MISC, norm_samples_sig$DRS, plot = T, ci  = T)

roc_5_db <- roc(norm_samples_sig$MISC[norm_samples_sig$group %in% c("MIS_C", "DB")], 
                norm_samples_sig$DRS[norm_samples_sig$group %in% c("MIS_C", "DB")],
                plot = T, ci  = T)

roc_5_dv <- roc(norm_samples_sig$MISC[norm_samples_sig$group %in% c("MIS_C", "DV")], 
                norm_samples_sig$DRS[norm_samples_sig$group %in% c("MIS_C", "DV")],
                plot = T, ci  = T)

roc_5_kd <- roc(norm_samples_sig$MISC[norm_samples_sig$group %in% c("MIS_C", "KD")], 
                norm_samples_sig$DRS[norm_samples_sig$group %in% c("MIS_C", "KD")],
                plot = T, ci  = T)

roc_plot_all_5 <- roc_plot(roc_5_all, title = "MIS-C vs. all disease groups", 'darkgrey')
roc_plot_db_5 <- roc_plot(roc_5_db, title = "MIS-C vs. bacterial", 'red2')
roc_plot_dv_5 <- roc_plot(roc_5_dv, title = "MIS-C vs. viral", 'navy')
roc_plot_kd_5 <- roc_plot(roc_5_kd, title = "MIS-C vs. KD", 'gold')

# calculate threshold sensitiviy and specificity
coords <- coords(roc=roc_5_all, input = 'threshold', transpose = FALSE, x = 'best')

drs_plt <- ggplot(norm_samples_sig, 
                  aes(x = group,  
                      y = DRS, 
                      fill = group))+
  theme_bw()+
  geom_boxplot(alpha = 0.6, 
               outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3), size = 0.2)+
  theme(legend.position = 'none')+
  labs(x = 'Disease group', y = "Disease Risk Score")+
  scale_fill_manual(values = c('#9310B9', '#FDE725', '#46ACFF', '#FF022A'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = coords$threshold, linetype = 'dashed')+
  ggtitle("Disease risk score - RNA-seq discovery set")

ggdraw() +
  draw_plot(drs_plt, x = 0, y = .4, width = 0.5, height = 0.6) +
  draw_plot(roc_plot_all_5, x = 0.5, y = .45, width = .5, height = .55) +
  draw_plot(roc_plot_kd_5, x = 0, y = 0, width = .333, height = .4) +
  draw_plot(roc_plot_dv_5, x = 0.3333, y = 0, width = .333, height = .4) +
  draw_plot(roc_plot_db_5, x = 0.6666, y = 0, width = .333, height = .4) +
  draw_plot_label(label = c("A", "B", "C", "D", 'E'), size = 15,
                  x = c(0, 0.5, 0, 0.33, 0.66), 
                  y = c(1, 1, 0.4, 0.4, 0.4))





