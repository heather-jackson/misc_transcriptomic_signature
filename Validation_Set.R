###############################################################
# Code for analysis of RT-qPCR validation data for MIS-C signature
###############################################################

library(stringr)
library(clipr)
library(tidyverse)
library(lubridate)
library(stringi)
library(caret)
library(pROC)
library(caret)
library(ggpubr)
library(tidyr)
library(sp)
library(cowplot)

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
# read in CT values averaged across replicates and GAPDH normalised

assays_averaged <- read.csv(file = 'CT_Values_RTqPCR_MISC_Signature.csv')
rownames(assays_averaged) <- assays_averaged$X
assays_averaged <- assays_averaged[,2:ncol(assays_averaged)]

###############################################################
# test performance
# perform logistic regression with cross validation to determine the coefficients used for signature evaluation
# define training control
set.seed(10)
train_control <- trainControl(method = "cv", 
                              number = 10, 
                              savePredictions = TRUE, 
                              classProbs = F)

five_gene_seen <- assays_averaged[!assays_averaged$group %in% c("DV_COVID", "CONTROL"),]

# train the model
model <- train(as.factor(MISC) ~ .,
               data =  five_gene_seen[,c(1:5, 8)],
               trControl = train_control,
               method = "glm",
               family=binomial())

# extract out the coefficients
mod_coefficients <- model$finalModel$coefficients
mod_coefficients <- mod_coefficients[2:length(mod_coefficients)]

# calculate weighted DRS
assays_averaged$DRS <- as.vector(as.matrix(assays_averaged[,match(names(mod_coefficients), 
                                                                          colnames(assays_averaged))]) %*% as.matrix(mod_coefficients))

roc_seen_cv_coef <- roc(assays_averaged$MISC[!assays_averaged$group %in% c("DV_COVID",  "CONTROL")], 
                        assays_averaged$DRS[!assays_averaged$group %in% c("DV_COVID",  "CONTROL")], 
                        ci = T)

roc_db_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("Full Mis-C Criteria", "DEFINITE BACTERIAL", "Non-Sterile DB")], 
                      assays_averaged$DRS[assays_averaged$group %in% c("Full Mis-C Criteria", "DEFINITE BACTERIAL", "Non-Sterile DB")], 
                      ci = T)

roc_dv_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("Full Mis-C Criteria", "DEFINITE VIRAL")], 
                      assays_averaged$DRS[assays_averaged$group %in% c("Full Mis-C Criteria", "DEFINITE VIRAL")], 
                      ci = T)

roc_kd_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("Full Mis-C Criteria", "KAWASAKI DISEASE")], 
                      assays_averaged$DRS[assays_averaged$group %in% c("Full Mis-C Criteria", "KAWASAKI DISEASE")], 
                      ci = T)

roc_all_seen_unseen_cv_coef <- roc(assays_averaged$MISC[!assays_averaged$group == "CONTROL"], 
                                   assays_averaged$DRS[!assays_averaged$group == "CONTROL"], 
                                   ci = T)

roc_all_unseen_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("DV_COVID",  "Full Mis-C Criteria")], 
                              assays_averaged$DRS[assays_averaged$group %in% c("DV_COVID",  "Full Mis-C Criteria")], 
                              ci = T)

roc_covid_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("DV_COVID", "Full Mis-C Criteria")], 
                         assays_averaged$DRS[assays_averaged$group %in% c("DV_COVID", "Full Mis-C Criteria")], 
                         ci = T)

roc_db_s_cv_coef <- roc(assays_averaged$MISC[assays_averaged$group %in% c("Non-Sterile DB", "Full Mis-C Criteria")], 
                        assays_averaged$DRS[assays_averaged$group %in% c("Non-Sterile DB", "Full Mis-C Criteria")], 
                        ci = T)

# create plots 
roc_5_coefficient_seen_p <- roc_plot(roc_seen_cv_coef, "MIS-C vs all seen", 'darkgrey')

roc_5_coefficient_db_p <- roc_plot(roc_db_cv_coef, "MIS-C vs bacterial", 'red2')

roc_5_coefficient_dv_p <- roc_plot(roc_dv_cv_coef, "MIS-C vs viral", 'navy')

roc_5_coefficient_kd_p <- roc_plot(roc_kd_cv_coef, "MIS-C vs KD", 'gold')

roc_5_coefficient_seen_unseen_p <- roc_plot(roc_all_seen_unseen_cv_coef, "MIS-C vs all disease groups", 'darkgrey')

roc_5_coefficient_unseen_p <- roc_plot(roc_all_unseen_cv_coef, "MIS-C vs all unseen", 'darkgrey')

roc_5_coefficient_covid_p <- roc_plot(roc_covid_cv_coef, "MIS-C vs COVID-19", 'lightblue')

roc_5_coefficient_db_s_p <- roc_plot(roc_db_s_cv_coef, "MIS-C vs non-sterile bacterial", 'red4')

# weighted DRS plot
# calculate threshold sensitiviy and specificity
coords <- coords(roc=roc_seen_cv_coef, input = 'threshold', transpose = FALSE, x = 'best')

assays_averaged$group[assays_averaged$group == 'DEFINITE BACTERIAL'] <- 'Bacterial'
assays_averaged$group[assays_averaged$group == 'DEFINITE VIRAL'] <- 'Viral'
assays_averaged$group[assays_averaged$group == 'KAWASAKI DISEASE'] <- 'KD'
assays_averaged$group[assays_averaged$group == 'DV_COVID'] <- 'COVID-19'
assays_averaged$group[assays_averaged$group == 'Full Mis-C Criteria'] <- 'MIS-C'

assays_averaged$group <- fct_relevel(assays_averaged$group, 
                                         "MIS-C", 
                                         "KD", 
                                         "Viral", 
                                         "Bacterial", 
                                         "COVID-19")

colkey = c("MIS-C" = '#9310B9', 
           "Viral" = '#46ACFF', 
           "KD" = '#FDE725', 
           "Bacterial" = '#FF022A', 
           "COVID-19" = '#00008b')

DRS_plot <- ggplot(assays_averaged[assays_averaged$group != "CONTROL",], aes(x = group, 
                                                                                     y = DRS, 
                                                                                     fill = group))+
  theme_bw()+
  geom_boxplot(alpha = 0.6, 
               outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3), size = 0.2)+
  theme(legend.position = 'none')+
  labs(x = 'Disease group', y = "Disease Risk Score")+
  scale_fill_manual(values = colkey)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = coords$threshold, linetype = 'dashed')+
  ggtitle("Disease risk score - RT-qPCR validation set")


# plot together 
ggdraw() +
  draw_plot(DRS_plot, x = 0, y = .4, width = 0.3, height = 0.6) +
  draw_plot(roc_5_coefficient_seen_unseen_p, x = 0.3, y = .5, width = .23, height = .5) +
  draw_plot(roc_5_coefficient_seen_p, x = 0.53, y = .5, width = .23, height = .5) +
  draw_plot(roc_5_coefficient_unseen_p, x = 0.76, y = .5, width = .23, height = .5) +
  draw_plot(roc_5_coefficient_kd_p, x = 0, y = 0, width = .2, height = .4) +
  draw_plot(roc_5_coefficient_dv_p, x = .2, y = 0, width = .2, height = .4) +
  draw_plot(roc_5_coefficient_db_p, x = 0.4, y = 0, width = .2, height = .4) +
  draw_plot(roc_5_coefficient_covid_p, x = 0.6, y = 0, width = .2, height = .4) +
  draw_plot(roc_5_coefficient_db_s_p, x = 0.8, y = 0, width = .2, height = .4) +
  draw_plot_label(label = c("A", "B", "C", "D", 'E', "F", "G", "H", "I"), size = 15,
                  x = c(0, 0.3, 0.53, 0.76,
                        0, 0.2, 0.4, 0.6, 0.8), 
                  y = c(1, 1, 1, 1, 0.4, 0.4, 0.4, 0.4, 0.4))


