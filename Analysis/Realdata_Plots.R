## ***************************************************************************
## Realdata_Plots.R #####
##
## Script name: Realdata_Plots
##
## Purpose of script: Plot the results from the real data analysis
##
## Author: Annika Swenne
##
## Date Created: 2022-08-29
##
## Version/ Date: V01 / 2022-08-29
##
## ***************************
##
## Modification Notes:
## V01/ 2022-08-29: 
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
library(data.table)
library(xtable)


set.seed(42)

## ***************************************************************************
load(file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Realdata/Res_real.RData")

npred <- 100
predictors <- colnames(res)[16:ncol(res)]

#read the data
data <- readRDS("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Realdata/selected_snps_cc.rds")
data$isced <- as.factor(data$isced)

#change the structure of the results for the vim to long format to visualize them with ggplot2
reslong <- res[which(!is.na(res$rs1620977)),]
if ("result.1" %in% colnames(res)){
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors, drop=c("result.1","sex","age","country","isced"))
} else {
  reslong <- reshape(reslong, varying=predictors, direction="long", v.names="vim", times=predictors)
}
reslong <- reslong[order(reslong$job.id),]
reslong <- subset(reslong, select=-id)
colnames(reslong)[names(reslong)=="time"] <- "predictor"

reslong$predictor <- factor(reslong$predictor, levels=unique(reslong$predictor))

#create names for the algorithms
reslong$rf <- NA
reslong$rf[which(reslong$algorithm=="rf_var" | reslong$algorithm=="freq_rf_var")] <- "Variance splitting"
reslong$rf[which(reslong$algorithm=="rf_zhao" | reslong$algorithm=="freq_rf_zhao")] <- "Zhao et al. (2012)"
reslong$rf[which(reslong$algorithm=="rf_res" | reslong$algorithm=="freq_rf_res")] <- "Residual variance splitting"
reslong$rf[which(reslong$algorithm=="rf_maxstat")] <- "Maximally selected residual rank statistic"
reslong$rf <- factor(reslong$rf, levels=unique(reslong$rf))

#extract top 10 median vim
medianvim <- setDT(reslong)[,list(Median=as.numeric(median(vim))), by=.(algorithm, predictor)]
medianvim <- medianvim[order(medianvim$Median, decreasing = TRUE),]
medianvimvar <- medianvim[which(medianvim$algorithm=="rf_var"),]
medianvimzhao <- medianvim[which(medianvim$algorithm=="rf_zhao"),]
medianvimmaxstat <- medianvim[which(medianvim$algorithm=="rf_maxstat"),]
topvimvar <- medianvimvar$predictor[c(1:10)]
topvimzhao <- medianvimzhao$predictor[c(1:10)]
topvimmaxstat <- medianvimmaxstat$predictor[c(1:10)]

#create dataset for the prediction accuracy
respred <- res[which(!is.na(res$result.1)),-which(colnames(res) %in% predictors),with=FALSE]
#change the job ids
respred$job.id[which(respred$job.id==71)] <- 1
respred$job.id[which(respred$job.id>=5)] <- respred$job.id[which(respred$job.id>=5)]+1
respred$job.id[which(respred$job.id==71)] <- 5
respred <- respred[order(respred$job.id),]
#add names
respred$rf <- NA
respred$rf[which(respred$algorithm=="pred_rf_base")] <- "Variance \nsplitting \nwithout SNPs"
respred$rf[which(respred$algorithm=="pred_rf_var")] <- "Variance \nsplitting"
respred$rf[which(respred$algorithm=="pred_rf_zhao")] <- "Zhao et al. \n(2012)"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="Meanout")] <- "Residual \nsplitting \nMeanout"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="Meanresid")] <- "Residual \nsplitting \nMeanresid"
respred$rf[which(respred$algorithm=="pred_rf_res" & respred$predleaf=="GLM")] <- "Residual \nsplitting \nGLM"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.05)] <- "Maxstat \nsplitting \nalpha: 0.05"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.1)] <- "Maxstat \nsplitting \nalpha: 0.1"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.2)] <- "Maxstat \nsplitting \nalpha: 0.2"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.3)] <- "Maxstat \nsplitting \nalpha: 0.3"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.4)] <- "Maxstat \nsplitting \nalpha: 0.4"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.5)] <- "Maxstat \nsplitting \nalpha: 0.5"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" & respred$alpha==0.6)] <- "Maxstat \nsplitting \nalpha: 0.6"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" &  respred$alpha==0.7)] <- "Maxstat \nsplitting \nalpha: 0.7"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" &  respred$alpha==0.8)] <- "Maxstat \nsplitting \nalpha: 0.8"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" &  respred$alpha==0.9)] <- "Maxstat \nsplitting \nalpha: 0.9"
respred$rf[which(respred$algorithm=="pred_rf_maxstat" &  respred$alpha==1)] <- "Maxstat \nsplitting \nalpha: 1"
respred$rf <- factor(respred$rf, levels=c("Variance \nsplitting \nwithout SNPs",
                                          "Variance \nsplitting",
                                          "Zhao et al. \n(2012)",
                                          "Maxstat \nsplitting \nalpha: 0.05",
                                          "Maxstat \nsplitting \nalpha: 0.1",
                                          "Maxstat \nsplitting \nalpha: 0.2",
                                          "Maxstat \nsplitting \nalpha: 0.3", 
                                          "Maxstat \nsplitting \nalpha: 0.4",
                                          "Maxstat \nsplitting \nalpha: 0.5", 
                                          "Maxstat \nsplitting \nalpha: 0.6",
                                          "Maxstat \nsplitting \nalpha: 0.7",
                                          "Maxstat \nsplitting \nalpha: 0.8",
                                          "Maxstat \nsplitting \nalpha: 0.9", 
                                          "Maxstat \nsplitting \nalpha: 1"))
setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA")
#setwd("C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Batchtools")


#*A) Study description ----
table(data$country)
#ITA EST CYP BEL SWE GER HUN ESP 
#626 296   0 210 424 622 445 374 
prop.table(table(data$country))
#ITA        EST        CYP        BEL        SWE        GER        HUN        ESP 
#0.20887554 0.09876543 0.00000000 0.07007007 0.14147481 0.20754087 0.14848182 0.12479146 

table(data$sex)
#male female 
#1484   1513
prop.table(table(data$sex))
#male    female 
#0.4951618 0.5048382

table(data$isced)
#1    2    3 
#185 1330 1482 
prop.table(table(data$isced))
#1         2         3 
#0.0617284 0.4437771 0.4944945 

mean(data$age) #11.41211
sd(data$age) #1.920415
median(data$age) #11.2
IQR(data$age) #3.5
min(data$age) #4
max(data$age) #16.2


#*29) Tuning alpha ----
ggplot(respred[1:11,], aes(x=alpha, y=result.1)) + 
  geom_point() +
  theme_bw() + 
  geom_hline(yintercept = min(respred$result.1[1:11]), linetype="dotted") +
  labs(x="Alpha", y="OOB Error")
ggsave("img/Figure29_alpha.png", height=4, width=7)


#*30) VIM ----
#**a) single plot ----
idvar <- which(reslong$algorithm=="rf_var")
idzhao <- which(reslong$algorithm=="rf_zhao")
idmax <- which(reslong$algorithm=="rf_maxstat")
ggplot(reslong[c(idvar,idzhao,idmax),], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  theme(text = element_text(size=7), axis.text.x = element_text(angle = 90)) +
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure30_vim.png", height=9, width=7.5)

#**b) individual plots ----
idvar <- which(reslong$algorithm=="rf_var")
ggplot(reslong[idvar,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  theme(text = element_text(size=7), axis.text.x = element_text(angle = 90)) +
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure30b_vim_var.png", height=4, width=7.5)

idzhao <- which(reslong$algorithm=="rf_zhao")
ggplot(reslong[idzhao,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  theme(text = element_text(size=7), axis.text.x = element_text(angle = 90)) +
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure30b_vim_zhao.png", height=4, width=7.5)

idmax <- which(reslong$algorithm=="rf_maxstat")
ggplot(reslong[idmax,], aes(x=predictor, y=vim)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~rf, ncol=1, scales="free") +
  theme_bw() + 
  theme(text = element_text(size=7), axis.text.x = element_text(angle = 90)) +
  labs(x="Predictor", y="Permutation importance")
ggsave("img/Figure30b_vim_max.png", height=4, width=7.5)

#**c) top median VIM ----
#create table for export
medianvim2 <- data.frame(rf_var=medianvimvar$predictor,
                         rf_zhao=medianvimzhao$predictor,
                         rf_maxstat=medianvimmaxstat$predictor)
print(xtable(medianvim2[c(1:10),], type="latex"),
      file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA/table/MedianVIM_topten.tex")
write.csv(medianvim2, quote=FALSE,
          file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA/table/MedianVIM.csv")
medianvim2 <- data.frame(rf_var1=medianvimvar$predictor[1:50],
                         rf_zhao1=medianvimzhao$predictor[1:50],
                         rf_maxstat1=medianvimmaxstat$predictor[1:50],
                         nr=seq(51,100,by=1),
                         rf_var2=medianvimvar$predictor[51:100],
                         rf_zhao2=medianvimzhao$predictor[51:100],
                         rf_maxstat2=medianvimmaxstat$predictor[51:100])
print(xtable(medianvim2, type="latex"),
      file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA/table/MedianVIM.tex")


#*31) MSE ----
id <- which(respred$alpha!=0.2)
ggplot(respred[-id,], aes(x=rf, y=result.1)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="OOB Error")
ggsave("img/Figure31_mse.png", height=3, width=7)


#*32) Runtime ----
load(file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Code/Realdata/Res_runtime.RData")
#create names for the algorithms
res$rf <- NA
res$rf[which(res$algorithm=="run_rf_var")] <- "Variance splitting"
res$rf[which(res$algorithm=="run_rf_zhao")] <- "Zhao et al. (2012)"
res$rf[which(res$algorithm=="run_rf_maxstat")] <- "Maximally selected residual rank statistic"
res$rf <- factor(res$rf, levels=unique(res$rf))
#calculate time in minutes
res$time_m <- res$time*10^-9/60

ggplot(res, aes(x=rf, y=time_m)) + 
  stat_boxplot(geom="errorbar", width = 0.3) +
  geom_boxplot(fill="lightgrey", notch=FALSE) + 
  theme_bw() + 
  labs(x="Algorithm", y="Runtime (m)")
ggsave("img/Figure33_runtime.png", height=3, width=7)

#create table for export
runtime <- setDT(res)[,list(Minimum=round(as.numeric(min(time_m)),digits=1),
                            Q1=round(as.numeric(quantile(time_m,probs=0.25)),digits=1),
                            Mean=round(as.numeric(mean(time_m)),digits=1),
                            Median=round(as.numeric(median(time_m)),digits=1),
                            Q3=round(as.numeric(quantile(time_m,probs=0.75)),digits=1),
                            Maximum=round(as.numeric(max(time_m)),digits=1)), 
                      by=.(rf)]
print(xtable(runtime, type="latex"),
      file="C:/Users/Annika Swenne/Documents/Uni_Bremen/MedicalBiometry/Masterarbeit/Latex_MA/table/Runtime.tex")

