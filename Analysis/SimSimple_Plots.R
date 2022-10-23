## ***************************************************************************
## SimSimple_Plots.R #####
##
## Script name: SimSimple_Plots
##
## Purpose of script: Plot the results from the simple simulation
##
## Author: Annika Swenne
##
## Date Created: 2022-07-28
##
## Version/ Date: V01 / 2022-07-28
##                V02 / 2022-09-02
##
## ***************************
##
## Modification Notes:
## V02/ 2022-09-02: Different plots 
##   
## ***************************
##
## Notes/ to do:
##
## ***************************************************************************

#load the required packages:  (uncomment as required)
library(batchtools)
library(ggplot2)
library(dplyr)
library(rpart)
library(rpart.plot) #for pretty plots

set.seed(42)

## ***************************************************************************
load(file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Batchtools/Res_simple.RData")

npred <- 10
predictors <- paste(rep("x", times=npred), seq(1,npred), sep="")

#change the structure of the results for the vim to long format to visualize them with ggplot2
reslong <- res[which(!is.na(res$x1)),]
if ("error" %in% colnames(res)){
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors, drop=c("z1","error","accuracy"))
} else {
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors)
}
reslong <- reslong[order(reslong$job.id),]
reslong <- subset(reslong, select=-id)
colnames(reslong)[names(reslong)=="time"] <- "predictor"

reslong$predictor <- factor(reslong$predictor, levels=unique(reslong$predictor))

#create names for the algorithms & leaf models
reslong$rf <- NA
reslong$rf[which(reslong$algorithm=="rf_var" | reslong$algorithm=="freq_rf_var")] <- "Variance splitting"
reslong$rf[which(reslong$algorithm=="rf_zhao" | reslong$algorithm=="freq_rf_zhao")] <- "Zhao et al. (2012)"
reslong$rf[which(reslong$algorithm=="rf_res" | reslong$algorithm=="freq_rf_res")] <- "Residual variance splitting"
reslong$rf[which(reslong$algorithm=="rf_maxstat" | reslong$algorithm=="freq_rf_maxstat")] <- "Maximally selected residual rank statistic"
reslong$rf <- factor(reslong$rf, levels=unique(reslong$rf))

reslong$leaf <- NA
reslong$leaf[which(reslong$predleaf=="Meanout")] <- "Leaf model: Meanout"
reslong$leaf[which(reslong$predleaf=="Meanresid")] <- "Leaf model: Meanresid"
reslong$leaf[which(reslong$predleaf=="GLM")] <- "Leaf model: GLM"
reslong$leaf <- factor(reslong$leaf, levels=unique(reslong$leaf))

#create data set for the prediction accuracy
respred <- res[which(!is.na(res$error)),-which(colnames(res) %in% predictors),with=FALSE]
respred$rf <- NA
respred$rf[which(respred$algorithm=="pred_rf_var")] <- "Variance \nsplitting"
respred$rf[which(respred$algorithm=="pred_rf_zhao")] <- "Zhao et al. \n(2012)"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="Meanout")] <- "Residual \nsplitting \nMeanout"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="Meanresid")] <- "Residual \nsplitting \nMeanresid"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="GLM")] <- "Residual \nsplitting \nGLM"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.05 & respred$predleaf=="Meanout")] <- "Maxstat \nsplitting \nMeanout"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.05 & respred$predleaf=="Meanresid")] <- "Maxstat \nsplitting \nMeanresid"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.05 & respred$predleaf=="GLM")] <- "Maxstat \nsplitting \nGLM"

respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.3 & respred$predleaf=="Meanout")] <- "Meanout \nalpha: 0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.3 & respred$predleaf=="Meanresid")] <- "Meanresid \nalpha: 0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.3 & respred$predleaf=="GLM")] <- "GLM \nalpha: 0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.5 & respred$predleaf=="Meanout")] <- "Meanout \nalpha: 0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.5 & respred$predleaf=="Meanresid")] <- "Meanresid \nalpha: 0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.5 & respred$predleaf=="GLM")] <- "GLM \nalpha: 0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.7 & respred$predleaf=="Meanout")] <- "Meanout \nalpha: 0.7"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.7 & respred$predleaf=="Meanresid")] <- "Meanresid \nalpha: 0.7"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.7 & respred$predleaf=="GLM")] <- "GLM \nalpha: 0.7"


respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.05 & respred$problem=="pred_simsimple_det")] <- "Maxstat \nsplitting \nalpha: 0.05"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.3 & respred$problem=="pred_simsimple_det")] <- "Maxstat \nsplitting \nalpha: 0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.5 & respred$problem=="pred_simsimple_det")] <- "Maxstat \nsplitting \nalpha: 0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.7 & respred$problem=="pred_simsimple_det")] <- "Maxstat \nsplitting \nalpha: 0.7"

respred$rf <- factor(respred$rf, levels=unique(respred$rf))

setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA")


#*0) Example tree ----
data("PimaIndiansDiabetes2", package = "mlbench")
tree <- rpart(diabetes ~ age + glucose + insulin, data=PimaIndiansDiabetes2, method="class")
#plot
png("img/ExampleTree.png", width=862, height=515)
rpart.plot(tree, box.palette = "Grays")
dev.off()


#*1) VIM Null case: Comparison RFs ----
idvar <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_var")
idzhao <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_zhao")
idres <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_res"
               & reslong$predleaf=="Meanout")
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.05 & reslong$predleaf=="Meanout")
ggplot(reslong[c(idvar,idzhao,idres,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure1_null_comp.png", height=8, width=7)


#*2) VIM Null case: Alpha ----
#predleaf="Meanout"
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$predleaf=="Meanout" & reslong$alpha>0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure2_null_alpha_meanout.png", height=6, width=7)


#*3) VIM Null case: Leaves Residual splitting ----
idres <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_res")
ggplot(reslong[idres,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure3_null_leaves_res.png", height=6, width=7)


#*4) VIM Null case: Leaves Maxstat splitting ----
#**a) alpha=0.05 ----
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure4_null_leaves_max_alpha005.png", height=6, width=7)

#**b) alpha=0.3 ----
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.3)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure4_null_leaves_max_alpha03.png", height=6, width=7)

#**c) alpha=0.5 ----
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.5)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure4_null_leaves_max_alpha05.png", height=6, width=7)

#**d) alpha=0.7 ----
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.7)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure4_null_leaves_max_alpha07.png", height=6, width=7)

#close-up
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.7 & reslong$predleaf=="Meanresid")
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller = label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure4_null_leaves_max_closeup.png", height=4, width=7)


#*5) Selection frequency Null case ----
idvar <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="freq_rf_var")
idzhao <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="freq_rf_zhao")
idres <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="freq_rf_res")
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="freq_rf_maxstat"
               & reslong$alpha==0.05)
#add the percentages
#prop.table(table(reslong$job.id[c(idvar,idzhao,idres,idmax)],reslong$vim[c(idvar,idzhao,idres,idmax)]))

ggplot(reslong[c(idvar,idzhao,idres,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Selection frequency")
ggsave("img/Figure5_null_select.png", height=8, width=7)


#*6) Selection frequency Null case: Alpha ----
idmax <- which(reslong$problem=="simsimple_null" & reslong$algorithm=="freq_rf_maxstat"
               & reslong$alpha>0.05)
#add the percentages
#prop.table(table(reslong$job.id[c(idvar,idzhao,idres,idmax)],reslong$vim[c(idvar,idzhao,idres,idmax)]))

ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  facet_wrap(~alpha, ncol=1, scales="free", labeller = label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Selection frequency")
ggsave("img/Figure6_null_select_alpha.png", height=6, width=7)


#*7) Effect estimates Null case: Residual splitting ----
#get data set for the effect estimates
load(file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Batchtools/Res_simple_effect.RData")
res$leaf <- NA
res$leaf[which(res$predleaf=="Meanout")] <- "Leaf model: Meanout"
res$leaf[which(res$predleaf=="Meanresid")] <- "Leaf model: Meanresid"
res$leaf[which(res$predleaf=="GLM")] <- "Leaf model: GLM"
res$leaf <- factor(res$leaf, levels=unique(reslong$leaf))

res$stop <- NA
res$stop[which(res$minbucket==10)] <- "Minbucket=10"
res$stop[which(res$minbucket==20)] <- "Minbucket=20"
res$stop <- factor(res$stop, levels=unique(res$stop))

#**a) minsplit & minbucket ----
idres1 <- which(res$problem=="simsimple_null" & res$algorithm=="beta_rf_res" 
                & res$predleaf=="Meanresid" & res$minsplit==10)
idres2 <- which(res$problem=="simsimple_null" & res$algorithm=="beta_rf_res" 
                & res$predleaf=="GLM" & res$minbucket==10)
#get the number of leaves
nleaves <- table(res[c(idres1,idres2),]$leaf, res[c(idres1,idres2),]$x1split)[-1,]
dat_text <- data.frame(
  label = paste(c(nleaves[1,1],nleaves[1,2],nleaves[2,1],nleaves[2,2]),"leaves", sep="\n"),
  x1split = c("no","yes","no","yes"),
  leaf = c(rep("Leaf model: GLM",times=2),rep("Leaf model: Meanresid",times=2)))
dat_text$leaf <- factor(dat_text$leaf, levels=unique(reslong$leaf))

ggplot(res[c(idres1,idres2),], aes(x1split, y=betaz1)) + 
  geom_violin(fill="lightgrey") + 
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1, linetype="dotted") +
  stat_summary(geom = "point", fun = "mean", shape = 4) +
  facet_wrap(~leaf, ncol=2) +
  geom_text(data=dat_text[which(dat_text$x1split=="no"),], mapping=aes(x=1.2,y=2,label=label), 
            hjust=0,vjust=0, size=3)+
  geom_text(data=dat_text[which(dat_text$x1split=="yes"),], mapping=aes(x=2.2,y=2,label=label), 
            hjust=0,vjust=0, size=3)+
  theme_bw() + 
  coord_flip() +
  labs(x="Affected by split at x1", y="Estimated regression coefficient for the confounder in the leaves")
ggsave("img/Figure7_null_effect.png", height=4, width=7)

#**b) minbucket ----
idres1 <- which(res$problem=="simsimple_null" & res$algorithm=="beta_rf_res" 
                & res$predleaf=="Meanresid" & res$minbucket==10)
idres2 <- which(res$problem=="simsimple_null" & res$algorithm=="beta_rf_res" 
                & res$predleaf=="GLM")
#get the number of leaves
nleaves <- table(res[c(idres1,idres2),]$leaf, res[c(idres1,idres2),]$x1split,
                 res[c(idres1,idres2),]$minbucket)
dat_text <- data.frame(
  label = paste(c(nleaves[2,1,1],nleaves[2,2,1],nleaves[2,1,2],nleaves[2,2,2],
                  nleaves[3,1,1],nleaves[3,2,1]),"leaves", sep="\n"),
  x1split = c("no","yes","no","yes","no","yes"),
  leaf = c(rep("Leaf model: GLM",times=4),rep("Leaf model: Meanresid",times=2)),
  stop = c(rep("Minbucket=10",times=2),rep("Minbucket=20",times=2),rep("Minbucket=10",times=2)))
dat_text$leaf <- factor(dat_text$leaf, levels=unique(reslong$leaf))

ggplot(res[c(idres1,idres2),], aes(x1split, y=betaz1)) + 
  geom_violin(fill="lightgrey") + 
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1, linetype="dotted") +
  stat_summary(geom = "point", fun = "mean", shape = 4) +
  facet_wrap(~leaf+stop, ncol=3) +
  geom_text(data=dat_text[which(dat_text$x1split=="no"),], mapping=aes(x=1.2,y=2,label=label), 
            hjust=0,vjust=0, size=3)+
  geom_text(data=dat_text[which(dat_text$x1split=="yes"),], mapping=aes(x=2.2,y=2,label=label), 
            hjust=0,vjust=0, size=3)+
  theme_bw() + 
  coord_flip() +
  labs(x="Affected by split at x1", y="Estimated regression coefficient for the confounder in the leaves")
ggsave("img/Figure7_null_effect_minbucket.png", height=4, width=7)


#*8) MSE Null case ----
id <- which(respred$problem=="pred_simsimple_null" & respred$algorithm!="pred_rf_maxstat")
idmax <- which(respred$problem=="pred_simsimple_null" & respred$algorithm=="pred_rf_maxstat"
            & respred$alpha==0.05)
ggplot(respred[c(id,idmax),], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="Mean squared error")
ggsave("img/Figure8_null_mse.png", height=3, width=7)


#*9) MSE Null case: Alpha ----
idmax <- which(respred$problem=="pred_simsimple_null" & respred$algorithm=="pred_rf_maxstat"
               & respred$alpha>0.05)
ggplot(respred[c(idmax),], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="Mean squared error")
ggsave("img/Figure9_null_mse_alpha.png", height=3, width=7)


#*10) VIM Power case: Comparison RFs ----
idvar <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_var")
idzhao <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_zhao")
idres <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_res"
               & reslong$predleaf=="Meanout")
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.05 & reslong$predleaf=="Meanout")
ggplot(reslong[c(idvar,idzhao,idres,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure10_power_comp.png", height=8, width=7)


#*11) VIM Power case: Alpha ----
#predleaf="Meanout"
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$predleaf=="Meanout" & reslong$alpha>0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller = label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure11_power_alpha_meanout.png", height=6, width=7)


#*12) VIM Power case: Leaves Residual splitting ----
idres <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_res")
ggplot(reslong[idres,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure12_power_leaves_res.png", height=6, width=7)


#*13) VIM Power case: Leaves Maxstat splitting ----
#**a) alpha=0.05 ----
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure13_power_leaves_max_alpha005.png", height=6, width=7)

#**b) alpha=0.3 ----
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.3)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure13_power_leaves_max_alpha03.png", height=6, width=7)

#**c) alpha=0.5 ----
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.5)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure13_power_leaves_max_alpha05.png", height=6, width=7)

#**d) alpha=0.7 ----
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.7)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~leaf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure13_power_leaves_max_alpha07.png", height=6, width=7)


#*14) VIM Power case: Sensitivity ----
idvar <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_var" 
               & reslong$betax2==0.15)
idzhao <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_zhao"
                & reslong$betax2==0.15)
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha==0.05 & reslong$betax2==0.15)
ggplot(reslong[c(idvar,idzhao,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure14_power_sensitivity.png", height=6, width=7)


#*15) VIM Power case: Sensitivity Alpha----
idmax <- which(reslong$problem=="simsimple_power" & reslong$algorithm=="rf_maxstat"
               & reslong$alpha<1 & reslong$alpha>0.05 & reslong$betax2==0.15)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller = label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure15_power_sensitivity_alpha.png", height=6, width=7)


#*16) MSE Power case ----
id <- which(respred$problem=="pred_simsimple_power" & respred$algorithm!="pred_rf_maxstat")
idmax <- which(respred$problem=="pred_simsimple_power" & respred$algorithm=="pred_rf_maxstat"
               & respred$alpha==0.05)
ggplot(respred[c(id,idmax),], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="Mean squared error")
ggsave("img/Figure16_power_mse.png", height=3, width=7)


#*17) MSE Power case: Alpha ----
idmax <- which(respred$problem=="pred_simsimple_power" & respred$algorithm=="pred_rf_maxstat"
               & respred$alpha>0.05)
ggplot(respred[c(idmax),], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="Mean squared error")
ggsave("img/Figure17_power_mse_alpha.png", height=3, width=7)


#*18) VIM Local ----
#**a) mtry=3 ----
id <- which(reslong$problem=="simsimple_local" & reslong$mtry==3 
            & reslong$algorithm != "rf_maxstat" & is.na(reslong$addx2))
idmax <- which(reslong$problem=="simsimple_local" & reslong$mtry==3 
            & reslong$algorithm == "rf_maxstat" & reslong$alpha==0.05 & is.na(reslong$addx2))
ggplot(reslong[c(id,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure18_local_mtry3.png", height=4, width=7)

#**b) mtry=10 ----
id <- which(reslong$problem=="simsimple_local" & reslong$mtry==10 
            & reslong$algorithm != "rf_maxstat" & is.na(reslong$addx2))
idmax <- which(reslong$problem=="simsimple_local" & reslong$mtry==10 
               & reslong$algorithm == "rf_maxstat" & reslong$alpha==0.05 & is.na(reslong$addx2))
ggplot(reslong[c(id,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure18_local_mtry10.png", height=4, width=7)

#**c) addx2=TRUE (mtry=10) ----
id <- which(reslong$problem=="simsimple_local" & reslong$mtry==10 
            & reslong$algorithm != "rf_maxstat" & reslong$addx2==TRUE)
idmax <- which(reslong$problem=="simsimple_local" & reslong$mtry==10 
               & reslong$algorithm == "rf_maxstat" & reslong$alpha==0.05 & reslong$addx2==TRUE)
ggplot(reslong[c(id,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure18_local_addx2.png", height=4, width=7)


#*19:) VIM Local: Alpha ----
#**a) mtry=3 ----
idmax <- which(reslong$problem=="simsimple_local" & reslong$mtry==3 
               & reslong$algorithm == "rf_maxstat" & reslong$alpha>0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure19_local_alpha_mtry3.png", height=4, width=7)

#**b) mtry=10 ----
idmax <- which(reslong$problem=="simsimple_local" & reslong$mtry==10 
               & reslong$algorithm == "rf_maxstat" & reslong$alpha>0.05)
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~alpha, ncol=1, scales="free", labeller=label_both) +
  theme_bw() + 
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure19_local_alpha_mtry10.png", height=4, width=7)


#*20) MSE Local ----
id <- which(respred$problem=="pred_simsimple_local")
ggplot(respred[id,], aes(x=rf, y=error)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  facet_wrap(~mtry, ncol=1, scales="fixed", labeller=label_both) +
  theme_bw() + 
  labs(x="Algorithm", y="Mean squared error")
ggsave("img/Figure20_local_mse.png", height=5, width=7)
